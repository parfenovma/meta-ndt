module MetaNDTGridapExt

using MetaNDT, Gridap, GridapGmsh, Gmsh, SparseArrays, LinearAlgebra
using Gridap.Geometry: get_grid, get_node_coordinates, get_cell_node_ids
using Gridap.FESpaces: get_cell_dof_ids
using Gridap.ReferenceFEs: get_node_and_comp_to_dof

const MESH_LOCK = ReentrantLock()

function mesh(s,path)
    lock(MESH_LOCK) do
        gmsh.isInitialized() == 0 || throw(ArgumentError(
            "MetaNDT cannot mesh while a caller-owned Gmsh session is active. " *
            "Finish that session first, or run the experiment in a separate process."))
        gmsh.initialize()
        try
            gmsh.option.setNumber("General.Terminal",0)
            gmsh.model.add("water_cylinder")
            o = gmsh.model.occ
            rect = o.addRectangle(-s.width/2,-s.height/2,0,s.width,s.height)
            if s.scenario == :empty
                fluid = [(2,rect)]
                solid = Tuple{Int,Int}[]
            else
                if isnothing(s.shape)
                    solid = [(2,o.addDisk(p[1],p[2],0,s.radius,s.radius)) for p in s.centers]
                else
                    solid = [(2,brick_surface(o,outline,s.mesh_size)) for outline in specimen_outlines(s)]
                end
                fluid,_ = o.cut([(2,rect)],solid,-1,true,false)
            end
            o.synchronize()
            gmsh.model.addPhysicalGroup(2,[t for (_,t) in fluid],-1,"water")
            !isempty(solid) && gmsh.model.addPhysicalGroup(2,[t for (_,t) in solid],-1,"solid")
            boundary = gmsh.model.getBoundary(fluid,true,false,false)
            outer = Int[]
            for (_,tag) in boundary
                bb = gmsh.model.getBoundingBox(1,tag)
                # A translated circle near an edge is still an interface, not an absorber.
                tolerance=1e-6
                on_x_edge=any(x -> abs(bb[1]-x)<tolerance && abs(bb[4]-x)<tolerance,
                    (-s.width/2,s.width/2))
                on_y_edge=any(y -> abs(bb[2]-y)<tolerance && abs(bb[5]-y)<tolerance,
                    (-s.height/2,s.height/2))
                if on_x_edge || on_y_edge
                    push!(outer,tag)
                end
            end
            length(outer) == 4 || error("unexpected outer rectangle boundary")
            gmsh.model.addPhysicalGroup(1,outer,-1,"outer")
            !isempty(solid) && gmsh.model.addPhysicalGroup(1,
                [t for (_,t) in gmsh.model.getBoundary(solid,true,false,false)],-1,"interface")
            gmsh.option.setNumber("Mesh.MeshSizeMin",s.mesh_size)
            gmsh.option.setNumber("Mesh.MeshSizeMax",s.mesh_size)
            gmsh.option.setNumber("Mesh.ElementOrder",1)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",0)
            gmsh.option.setNumber("Mesh.RandomSeed",1)
            gmsh.model.mesh.generate(2)
            gmsh.write(path)
            return GmshDiscreteModel(gmsh)
        finally
            gmsh.finalize()
        end
    end
end

function brick_surface(o,outline,mesh_size)
    points=[o.addPoint(p[1],p[2],0,mesh_size) for p in outline]
    curves=[o.addLine(points[i],points[mod1(i+1,length(points))]) for i in eachindex(points)]
    o.addPlaneSurface([o.addCurveLoop(curves)])
end

function nodal_data(V,Ω,reffe,offset=0)
    grid = get_grid(Ω)
    coords = get_node_coordinates(grid)
    cells = get_cell_node_ids(grid)
    used = sort!(unique!(reduce(vcat,[collect(c) for c in cells])))
    index = Dict(n=>i for (i,n) in enumerate(used))
    xy = hcat([Float64[p[1],p[2]] for p in coords[used]]...)
    tris = hcat([Int[index[n] for n in c] for c in cells]...)
    nc = get_node_and_comp_to_dof(reffe)
    components = nc[1] isa Integer ? 1 : length(nc[1])
    dofs = zeros(Int,components,length(used))
    cdofs = get_cell_dof_ids(V)
    for (cell,nodes) in enumerate(cells), j in eachindex(nodes), k in 1:components
        localid = components == 1 ? nc[j] : nc[j][k]
        dofs[k,index[nodes[j]]] = cdofs[cell][localid]+offset
    end
    xy,tris,dofs
end

function sensor_sampler(xy,tris,dofs,positions,n)
    # Spatial bins keep section/field sampling local rather than scanning the
    # complete FE mesh for each point. Weights are still exact P1 barycentrics.
    xmin,ymin=minimum(xy;dims=2)[:,1]
    xmax,ymax=maximum(xy;dims=2)[:,1]
    spacing=2sqrt((xmax-xmin)*(ymax-ymin)/size(tris,2))
    key(p)=(floor(Int,(p[1]-xmin)/spacing),floor(Int,(p[2]-ymin)/spacing))
    bins=Dict{Tuple{Int,Int},Vector{Int}}()
    for c in axes(tris,2)
        ids=tris[:,c]; px=xy[1,ids]; py=xy[2,ids]
        lo=key((minimum(px),minimum(py))); hi=key((maximum(px),maximum(py)))
        for ix in lo[1]:hi[1], iy in lo[2]:hi[2]
            push!(get!(bins,(ix,iy),Int[]),c)
        end
    end
    I=Int[]; J=Int[]; W=Float64[]
    for (i,p) in enumerate(positions)
        found=false
        for c in get(bins,key(p),Int[])
            ids=tris[:,c]; ax,ay=xy[:,ids[1]]
            bx,by=xy[:,ids[2]].-(ax,ay); cx,cy=xy[:,ids[3]].-(ax,ay)
            qx,qy=p.-(ax,ay); det=bx*cy-by*cx
            z1=(qx*cy-qy*cx)/det; z2=(bx*qy-by*qx)/det
            w = (1-z1-z2,z1,z2)
            if minimum(w) >= -1e-9
                for k in 1:3
                    push!(I,i); push!(J,dofs[1,ids[k]]); push!(W,w[k])
                end
                found=true; break
            end
        end
        found || throw(ArgumentError("sample $(i) lies outside the meshed domain"))
    end
    sparse(I,J,W,length(positions),n)
end

function assemble(s,path)
    model = mesh(s,path)
    Ωf=Triangulation(model; tags="water"); df=Measure(Ωf,4)
    Γo=BoundaryTriangulation(model; tags="outer"); do_=Measure(Γo,4)
    rf=ReferenceFE(TRI,lagrangian,Float64,1)
    Vf=TestFESpace(Ωf,rf; conformity=:H1); Uf=TrialFESpace(Vf)
    ρ=s.water_density; c=s.water_speed; L=s.width; F=c*L
    g(x)=exp(-((x[1]-s.source[1])^2+(x[2]-s.source[2])^2)/(2s.source_width^2))
    source_profile=CellField(g,Ωf)
    normalization=sum(∫(source_profile)*df)
    if s.scenario == :elastic
        Ωs=Triangulation(model; tags="solid"); ds=Measure(Ωs,4)
        rs=ReferenceFE(TRI,lagrangian,VectorValue{2,Float64},1)
        Vs=TestFESpace(Ωs,rs; conformity=:H1); Us=TrialFESpace(Vs)
        Y=MultiFieldFESpace([Vf,Vs]); X=MultiFieldFESpace([Uf,Us])
        Γ=InterfaceTriangulation(Ωs,Ωf); di=Measure(Γ,4)
        n=get_normal_vector(Γ).⁺
        λ=s.young_modulus*s.poisson_ratio/((1+s.poisson_ratio)*(1-2s.poisson_ratio))
        μ=s.young_modulus/(2*(1+s.poisson_ratio))
        me((φ,u),(q,v))=∫((ρ*F^2/c^2)*φ*q)*df + ∫((s.solid_density*L^2)*(u⋅v))*ds
        ke((φ,u),(q,v))=∫((ρ*F^2)*(∇(φ)⋅∇(q)))*df +
            ∫(L^2*(λ*tr(ε(u))*tr(ε(v))+2μ*(ε(u)⊙ε(v))))*ds
        cae((φ,u),(q,v))=∫((ρ*F^2/c)*φ*q)*do_
        ci((φ,u),(q,v))=∫((ρ*F*L)*(q.⁻*(u.⁺⋅n)-φ.⁻*(v.⁺⋅n)))*di
        le((q,v))=∫((-ρ*F/normalization)*source_profile*q)*df
        M=assemble_matrix(me,X,Y); K=assemble_matrix(ke,X,Y)
        Ca=assemble_matrix(cae,X,Y); Ci=assemble_matrix(ci,X,Y)
        b=assemble_vector(le,Y)
        sxy,stri,sdofs=nodal_data(Vs,Ωs,rs,num_free_dofs(Vf))
    else
        m(φ,q)=∫((ρ*F^2/c^2)*φ*q)*df
        k(φ,q)=∫((ρ*F^2)*(∇(φ)⋅∇(q)))*df
        ca(φ,q)=∫((ρ*F^2/c)*φ*q)*do_
        l(q)=∫((-ρ*F/normalization)*source_profile*q)*df
        M=assemble_matrix(m,Uf,Vf); K=assemble_matrix(k,Uf,Vf)
        Ca=assemble_matrix(ca,Uf,Vf); Ci=spzeros(size(M)...)
        b=assemble_vector(l,Vf)
        sxy=zeros(2,0); stri=zeros(Int,3,0); sdofs=zeros(Int,2,0)
    end
    !s.absorbing_boundary && (Ca=spzeros(size(M)...))
    fxy,ftri,fdofs=nodal_data(Vf,Ωf,rf)
    P=sensor_sampler(fxy,ftri,fdofs,s.sensors,size(M,1))
    (;M,K,Ca,Ci,b,fxy,ftri,fdofs,sxy,stri,sdofs,P,F,L,nf=num_free_dofs(Vf))
end

ε(u)=symmetric_gradient(u)

function pulse(s,t)
    T=s.pulse_cycles/s.frequency
    0 <= t <= T ? s.source_strength*sin(2π*s.frequency*t)*sinpi(t/T)^2 : 0.0
end

function MetaNDT.simulate(s::CylinderExperiment; mesh_path=nothing)
    # Snapshot mutable inputs so later edits cannot change a result's provenance.
    s=remake_experiment(s)
    started=time_ns()
    tmp=mktempdir()
    try
        path=isnothing(mesh_path) ? joinpath(tmp,"experiment.msh") : mesh_path
        a=assemble(s,path)
        assembled=(time_ns()-started)/1e9
        ratio=s.duration/s.time_step
        N=ceil(Int,ratio-8eps(ratio)); dt=s.duration/N
        times=collect(range(0,s.duration;length=N+1))
        frames=unique(round.(Int,range(1,N+1;length=min(s.frame_count,N+1))))
        nf=length(frames); n=size(a.M,1)
        p=zeros(size(a.fxy,2),nf); ux=zeros(size(a.sxy,2),nf); uy=similar(ux)
        traces=zeros(length(s.sensors),N+1)
        energy=zeros(N+1); solid_energy=zeros(N+1); work=zeros(N+1); absorbed=zeros(N+1)
        C=a.Ca+a.Ci
        A=lu(a.M+(dt/2)*C+(dt^2/4)*a.K)
        d=zeros(n); v=zeros(n); acc=zeros(n); fprev=zeros(n)
        solid=(a.nf+1):n
        Ms=a.M[solid,solid]; Ks=a.K[solid,solid]
        fi=1; peak_displacement=0.0
        solve_start=time_ns()
        for step in 1:(N+1)
            if step > 1
                dp=d+dt*v+(dt^2/4)*acc
                vp=v+(dt/2)*acc
                f=a.b*pulse(s,times[step])
                anew=A \ (f-C*vp-a.K*dp)
                dnew=dp+(dt^2/4)*anew
                vnew=vp+(dt/2)*anew
                vm=(v+vnew)/2
                work[step]=work[step-1]+dt*dot(vm,(fprev+f)/2)
                absorbed[step]=absorbed[step-1]+dt*dot(vm,a.Ca*vm)
                d=dnew; v=vnew; acc=anew; fprev=f
            end
            traces[:,step]=-s.water_density*a.F*(a.P*v)
            energy[step]=(dot(v,a.M*v)+dot(d,a.K*d))/2
            if !isempty(solid)
                solid_energy[step]=(dot(v[solid],Ms*v[solid])+dot(d[solid],Ks*d[solid]))/2
                peak_displacement=max(peak_displacement,
                    a.L*maximum(hypot.(d[a.sdofs[1,:]],d[a.sdofs[2,:]])))
            end
            if fi <= nf && step == frames[fi]
                p[:,fi]=-s.water_density*a.F*v[vec(a.fdofs)]
                if !isempty(a.sdofs)
                    ux[:,fi]=a.L*d[a.sdofs[1,:]]; uy[:,fi]=a.L*d[a.sdofs[2,:]]
                end
                fi+=1
            end
        end
        denom=max(maximum(abs.(work)),eps(Float64))
        diagnostics=Dict(
            "assembly_seconds"=>assembled,
            "time_loop_seconds"=>(time_ns()-solve_start)/1e9,
            "total_seconds"=>(time_ns()-started)/1e9,
            "dofs"=>Float64(n), "time_step_s"=>dt,
            "energy_balance_relative"=>maximum(abs.(energy+absorbed-work))/denom,
            "interface_skew_relative"=>norm(a.Ci+transpose(a.Ci))/max(norm(a.Ci),eps()),
            "interface_norm"=>norm(a.Ci),
            "peak_displacement_m"=>peak_displacement,
            "peak_sensor_pressure_Pa"=>maximum(abs.(traces)),
            "water_wavelength_elements"=>s.water_speed/(s.frequency*s.mesh_size),
            "cylinder_count"=>Float64(s.scenario == :empty || !isnothing(s.shape) ? 0 : length(s.centers)),
            "inclusion_count"=>Float64(s.scenario == :empty ? 0 : length(s.centers)),
        )
        CylinderResult(s,times,traces,times[frames],a.fxy,a.ftri,p,a.sxy,a.stri,
                       ux,uy,energy,solid_energy,work,absorbed,diagnostics)
    finally
        rm(tmp;recursive=true,force=true)
    end
end

include("solid_lens_fem.jl")

end
