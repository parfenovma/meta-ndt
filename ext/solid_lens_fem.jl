# Included inside MetaNDTGridapExt: reuse its meshing lock and FE-array adapters.
function solid_mesh(s,path)
    lock(MESH_LOCK) do
        gmsh.isInitialized()==0 || throw(ArgumentError("finish the caller-owned Gmsh session first"))
        gmsh.initialize()
        try
            gmsh.option.setNumber("General.Terminal",0)
            gmsh.model.add("connected_solid_lens"); o=gmsh.model.occ
            L=first(s.bricks).length
            if s.scenario == :solid
                domains=[(2,o.addRectangle(-s.width/2,-s.height/2,0,s.width,s.height))]
            else
                left=o.addRectangle(-s.width/2,-s.height/2,0,(s.width-L)/2,s.height)
                right=o.addRectangle(L/2,-s.height/2,0,(s.width-L)/2,s.height)
                ligaments=[(2,brick_surface(o,p,s.near_mesh_size)) for p in solid_lens_outlines(s)]
                domains,_=o.fuse([(2,left),(2,right)],ligaments)
            end
            o.synchronize()
            gmsh.model.addPhysicalGroup(2,last.(domains),-1,"solid")
            outer=Int[]; free=Int[]
            for (_,tag) in gmsh.model.getBoundary(domains,true,false,false)
                bb=gmsh.model.getBoundingBox(1,tag); tol=1e-6
                vertical=any(x->abs(bb[1]-x)<tol && abs(bb[4]-x)<tol,(-s.width/2,s.width/2))
                horizontal=any(y->abs(bb[2]-y)<tol && abs(bb[5]-y)<tol,(-s.height/2,s.height/2))
                push!(vertical || (s.absorb_sides && horizontal) ? outer : free,tag)
            end
            gmsh.model.addPhysicalGroup(1,outer,-1,"outer")
            !isempty(free) && gmsh.model.addPhysicalGroup(1,free,-1,"free")
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",0)
            gmsh.option.setNumber("Mesh.MeshSizeMin",s.near_mesh_size)
            gmsh.option.setNumber("Mesh.MeshSizeMax",s.mesh_size)
            field=gmsh.model.mesh.field.add("Box")
            for (key,value) in (("VIn",s.near_mesh_size),("VOut",s.mesh_size),
                ("XMin",-L/2-.002),("XMax",L/2+.002),
                ("YMin",-s.height/2),("YMax",s.height/2),("Thickness",.002))
                gmsh.model.mesh.field.setNumber(field,key,value)
            end
            gmsh.model.mesh.field.setAsBackgroundMesh(field)
            gmsh.option.setNumber("Mesh.ElementOrder",1)
            gmsh.option.setNumber("Mesh.RandomSeed",1)
            gmsh.model.mesh.generate(2); gmsh.write(path)
            GmshDiscreteModel(gmsh)
        finally
            gmsh.finalize()
        end
    end
end

function solid_assemble(s,path)
    model=solid_mesh(s,path)
    Ω=Triangulation(model;tags="solid"); dv=Measure(Ω,4)
    Γ=BoundaryTriangulation(model;tags="outer"); db=Measure(Γ,4); normal=get_normal_vector(Γ)
    reffe=ReferenceFE(TRI,lagrangian,VectorValue{2,Float64},1)
    V=TestFESpace(Ω,reffe;conformity=:H1); U=TrialFESpace(V)
    μ=s.young_modulus/(2*(1+s.poisson_ratio))
    λ=s.young_modulus*s.poisson_ratio/((1+s.poisson_ratio)*(1-2s.poisson_ratio))
    cp,cs=solid_wave_speeds(s)
    mass(u,v)=∫(s.density*(u⋅v))*dv
    stiff(u,v)=∫(λ*tr(ε(u))*tr(ε(v))+2μ*(ε(u)⊙ε(v)))*dv
    damping(u,v)=∫(s.density*(cs*(u⋅v)+(cp-cs)*(u⋅normal)*(v⋅normal)))*db
    # A finite-width plane-like body force, inside the left solid bulk.
    g(x)=exp(-(x[1]-s.source_x)^2/(2s.source_width^2)) *
        (s.absorb_sides ? cospi(x[2]/s.height)^2 : 1.0)
    profile=CellField(g,Ω); normalization=sum(∫(profile)*dv)
    direction=s.polarization == :shear ? VectorValue(0.0,1.0) : VectorValue(1.0,0.0)
    load(v)=∫((profile/normalization)*(direction⋅v))*dv
    M=assemble_matrix(mass,U,V); K=assemble_matrix(stiff,U,V)
    C=assemble_matrix(damping,U,V); force=assemble_vector(load,V)
    xy,tri,dofs=nodal_data(V,Ω,reffe)
    ys=collect(range(-.45s.height,.45s.height;length=301))
    positions=[(first(s.bricks).length/2+s.focus_distance,y) for y in ys]
    Px=sensor_sampler(xy,tri,view(dofs,1:1,:),positions,size(M,1))
    Py=sensor_sampler(xy,tri,view(dofs,2:2,:),positions,size(M,1))
    (;model,M,K,C,force,xy,tri,dofs,ys,Px,Py,area=sum(∫(1.0)*dv))
end

function MetaNDT.simulate(s::SolidLensExperiment;mesh_path=nothing)
    s=remake_experiment(s); started=time_ns(); tmp=mktempdir()
    try
        a=solid_assemble(s,isnothing(mesh_path) ? joinpath(tmp,"solid.msh") : mesh_path)
        ω=2π*s.frequency; f=s.source_strength*a.force
        u=(a.K-ω^2*a.M-im*ω*a.C)\f # exp(-iωt)
        power=.5real(im*ω*dot(u,f)); outgoing=.5ω^2*real(dot(u,a.C*u))
        profile=.25s.density*ω^2*(abs2.(a.Px*u)+abs2.(a.Py*u))
        right_face=first(s.bricks).length/2; available=s.width/2-right_face
        xs=collect(range(right_face+min(.005,s.focus_distance/2),
            s.width/2-min(.015,available/4);length=41))
        points=[(x,y) for y in a.ys for x in xs]
        SX=sensor_sampler(a.xy,a.tri,view(a.dofs,1:1,:),points,length(u))
        SY=sensor_sampler(a.xy,a.tri,view(a.dofs,2:2,:),points,length(u))
        field=reshape(.25s.density*ω^2*(abs2.(SX*u)+abs2.(SY*u)),length(xs),length(a.ys))
        widths=[MetaNDT.profile_width(a.ys,view(field,j,:)) for j in eachindex(xs)]
        peaks=vec(maximum(field;dims=2))
        residual=norm((a.K-ω^2*a.M-im*ω*a.C)*u-f)/max(norm(f),eps())
        diag=Dict("input_power_W_m"=>power,"outgoing_power_W_m"=>outgoing,
            "power_balance_relative"=>abs(power-outgoing)/max(abs(power),eps()),
            "linear_residual_relative"=>residual,"solid_area_m2"=>a.area,
            "dofs"=>Float64(length(u)),"total_seconds"=>(time_ns()-started)/1e9,
            "profile_peak_y_m"=>a.ys[argmax(profile)],
            "profile_fwhm_m"=>MetaNDT.profile_width(a.ys,profile),
            "profile_peak_per_input_power"=>maximum(profile)/max(power,eps()))
        SolidLensResult(s,a.xy,a.tri,u[a.dofs[1,:]],u[a.dofs[2,:]],a.ys,profile,xs,widths,peaks,diag)
    finally
        rm(tmp;recursive=true,force=true)
    end
end

function MetaNDT.design_stg_lens(s::SolidLensExperiment;
    amplitudes=collect(0.0:.0006:.0048),periods=[1,2,3])
    s=remake_experiment(s;scenario=:stg)
    b=first(s.bricks); candidates=SinusoidalBrick[]; transfer=ComplexF64[]
    # Calibration is a real one-cell elastic waveguide calculation, not a
    # hand-inserted delay. The complete wall still needs a separate full solve.
    base=remake_experiment(s;width=.060,height=b.height,bricks=[b],source_x=-.023,
        missing_cells=Int[],source_width=.001,focus_distance=.012,absorb_sides=false,scenario=:solid)
    baseline=simulate(base)
    component=s.polarization == :shear ? 2 : 1
    function sample(r)
        dofs=reshape(collect(1:size(r.nodes,2)),1,:)
        S=sensor_sampler(r.nodes,r.triangles,dofs,[(b.length/2+.012,0.0)],size(r.nodes,2))
        only(S*(component==2 ? r.displacement_y : r.displacement_x))
    end
    reference=sample(baseline)
    abs(reference)>0 || throw(ArgumentError("zero calibration source"))
    for N in periods, A in amplitudes
        cell=SinusoidalBrick(length=b.length,height=b.height,amplitude=A,periods=N,
            end_margin=b.end_margin,samples=24N)
        cal=simulate(remake_experiment(base;scenario=:stg,bricks=[cell]))
        push!(candidates,cell); push!(transfer,sample(cal)/reference)
    end
    cp,cs=solid_wave_speeds(s); k=2π*s.frequency/(component==2 ? cs : cp)
    centers=[(j-(length(s.bricks)+1)/2)*b.height for j in eachindex(s.bricks)]
    target=[-k*(hypot(s.focus_distance,y)-s.focus_distance) for y in centers]
    usable=findall(t->isfinite(abs(t)) && abs(t)>.15,transfer)
    isempty(usable) && throw(ArgumentError("no transmitting STG candidates"))
    best=Inf; selected=Int[]; errors=Float64[]
    # An irrelevant common phase is free; choose it to avoid blocked cells.
    for offset in range(-π,π;length=73)
        ids=[usable[argmin([abs(angle(transfer[i]*cis(-φ-offset)))^2+
            .04/max(abs(transfer[i]),.05)^2 for i in usable])] for φ in target]
        err=[angle(transfer[i]*cis(-φ-offset)) for (i,φ) in zip(ids,target)]
        score=sum(abs2,err)+sum(.04/max(abs(transfer[i]),.05)^2 for i in ids)
        if score<best; best=score; selected=ids; errors=err; end
    end
    designed=remake_experiment(s;bricks=candidates[selected])
    LensDesign(designed,candidates,transfer,selected,errors)
end

function MetaNDT.simulate_pulse(s::SolidLensExperiment;retain_history::Bool=false)
    s=remake_experiment(s); tmp=mktempdir(); started=time_ns()
    try
        a=solid_assemble(s,joinpath(tmp,"pulse.msh"))
        N=ceil(Int,s.duration/s.time_step); dt=s.duration/N
        times=collect(range(0,s.duration;length=N+1))
        frames=unique(round.(Int,range(1,N+1;length=min(s.frame_count,N+1))))
        ux=zeros(size(a.xy,2),length(frames)); uy=similar(ux)
        hx=zeros(Float32,size(a.xy,2),retain_history ? N+1 : 0); hy=similar(hx)
        energy=zeros(N+1); work=zeros(N+1); outgoing=zeros(N+1); section=zeros(length(a.ys))
        A=lu(a.M+dt/2*a.C+dt^2/4*a.K)
        d=zeros(size(a.M,1)); v=copy(d); acc=copy(d); fprev=copy(d); fi=1
        for step in 1:N+1
            if step>1
                dp=d+dt*v+dt^2/4*acc; vp=v+dt/2*acc
                drive=MetaNDT.solid_pulse_drive(s,times[step])
                f=s.source_strength*drive*a.force
                anew=A\(f-a.C*vp-a.K*dp)
                dn=dp+dt^2/4*anew; vn=vp+dt/2*anew; vm=(v+vn)/2
                work[step]=work[step-1]+dt*dot(vm,(fprev+f)/2)
                outgoing[step]=outgoing[step-1]+dt*dot(vm,a.C*vm)
                section+=dt*.25s.density*((a.Px*v).^2+(a.Py*v).^2+(a.Px*vn).^2+(a.Py*vn).^2)
                d=dn; v=vn; acc=anew; fprev=f
            end
            energy[step]=.5(dot(v,a.M*v)+dot(d,a.K*d))
            if retain_history
                hx[:,step]=d[a.dofs[1,:]]; hy[:,step]=d[a.dofs[2,:]]
            end
            if fi<=length(frames) && step==frames[fi]
                ux[:,fi]=d[a.dofs[1,:]]; uy[:,fi]=d[a.dofs[2,:]]; fi+=1
            end
        end
        diag=Dict("energy_balance_relative"=>maximum(abs.(energy+outgoing-work))/max(maximum(abs,work),eps()),
            "time_step_s"=>dt,"input_work_J_m"=>last(work),"total_seconds"=>(time_ns()-started)/1e9,
            "dofs"=>Float64(length(d)),"retained_history_bytes"=>Float64(sizeof(hx)+sizeof(hy)))
        SolidPulseResult(s,a.xy,a.tri,times[frames],ux,uy,times,energy,work,outgoing,a.ys,section,diag,hx,hy)
    finally
        rm(tmp;recursive=true,force=true)
    end
end
