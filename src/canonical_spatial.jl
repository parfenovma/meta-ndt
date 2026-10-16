# Positive-tensor acoustic FEM for the ground cloak.

function _cloak_assembly(x,y,variant)
    nx=length(x); ny=length(y); mass=zeros(nx*ny); I=Int[]; J=Int[]; V=Float64[]
    ρ0=1.2; B0=ρ0*343^2; a=.056; b=.168; c=.112
    for j in 1:ny-1,i in 1:nx-1
        ll=i+(j-1)*nx; lr=ll+1; ul=ll+nx; ur=ul+1
        for ids in ((ll,lr,ur),(ll,ur,ul))
            xx=[x[(q-1)%nx+1] for q in ids]; yy=[y[(q-1)÷nx+1] for q in ids]
            xc=sum(xx)/3; yc=sum(yy)/3
            variant!=:flat && yc<a*max(0,1-abs(xc)/b) && continue
            area=abs((xx[2]-xx[1])*(yy[3]-yy[1])-(xx[3]-xx[1])*(yy[2]-yy[1]))/2
            G=[yy[2]-yy[3] yy[3]-yy[1] yy[1]-yy[2];
               xx[3]-xx[2] xx[1]-xx[3] xx[2]-xx[1]]/(2area)
            invρ=Matrix{Float64}(LinearAlgebra.I,2,2)/ρ0; B=B0
            if variant==:cloak && yc<c*max(0,1-abs(xc)/b)
                # Popa Eq.(3): principal axes rotate oppositely in each half.
                θ=sign(xc)*deg2rad(66.9)
                R=[cos(θ) -sin(θ);sin(θ) cos(θ)]
                invρ=R*Diagonal([1/(5.7ρ0),1/(1.1ρ0)])*R'; B=1.25B0
            end
            localK=area*(G'*invρ*G)
            for u in 1:3
                mass[ids[u]]+=area/(3B)
                for v in 1:3
                    push!(I,ids[u]); push!(J,ids[v]); push!(V,localK[u,v])
                end
            end
        end
    end
    sparse(I,J,V,nx*ny,nx*ny),mass
end

function _cloak_record(variant;nx=161,ny=101)
    variant in (:flat,:bump,:cloak) || throw(ArgumentError("invalid ground specimen"))
    x=collect(range(-.6,.6;length=nx)); y=collect(range(0,.75;length=ny))
    K,M=_cloak_assembly(x,y,variant); valid=M.>0; invM=[m>0 ? inv(m) : 0. for m in M]
    h=min(x[2]-x[1],y[2]-y[1]); dt=.22h/343; nt=ceil(Int,.003/dt)
    # Initial incoming finite Gaussian packet, identical for all three specimens.
    # Its center reaches the bump at about 1.2 ms; it is not a harmonic phase film.
    p=zeros(nx*ny); v=zeros(nx*ny); damping=zeros(nx*ny); c=343.; f=3000.
    direction=[inv(sqrt(2)),-inv(sqrt(2))]; center=[-.31,.35]
    for j in 1:ny,i in 1:nx
        q=i+(j-1)*nx; dx=x[i]-center[1]; dy=y[j]-center[2]
        s=direction[1]*dx+direction[2]*dy; z=-direction[2]*dx+direction[1]*dy
        envelope=exp(-(s/.07)^2-(z/.105)^2)
        p[q]=envelope*cos(2π*f*s/c)
        v[q]=c*envelope*(2s/.07^2*cos(2π*f*s/c)+2π*f/c*sin(2π*f*s/c))
        distance=min(x[i]-first(x),last(x)-x[i],last(y)-y[j])
        damping[q]=5000*max(0,1-distance/.11)^2
        valid[q] || (p[q]=0;v[q]=0)
    end
    ix=collect(1:2:nx); iy=collect(1:2:ny); selected=[i+(j-1)*nx for j in iy for i in ix]
    stride=4; steps=collect(0:stride:nt); fields=zeros(Float32,length(selected),length(steps))
    fields[:,1]=p[selected]; force=zeros(length(p)); mul!(force,K,p)
    vhalf=v-(dt/2).*(-force.*invM-damping.*v)
    maxenergy=0.; initial=(sum(M.*v.^2)+dot(p,K*p))/2
    for step in 1:nt
        mul!(force,K,p)
        @. vhalf=((1-dt*damping/2)*vhalf-dt*force*invM)/(1+dt*damping/2)
        @. p=p+dt*vhalf
        if step%stride==0
            fields[:,step÷stride+1]=p[selected]
            maxenergy=max(maxenergy,sum(M.*vhalf.^2)/2+dot(p,K*p)/2)
        end
    end
    diag=Dict("frequency_Hz"=>f,"time_step_s"=>dt,"measurement_step_s"=>stride*dt,
        "initial_wave_energy"=>initial,"max_wave_energy_ratio"=>maxenergy/initial,
        "mesh_nodes"=>Float64(nx*ny))
    CanonicalRecord(x[ix],y[iy],steps.*dt,fields,BitVector(valid[selected]),zeros(Float32,0,0),diag)
end

function canonical_spatial_record(kind::Symbol,variant::Symbol;kwargs...)
    kind==:cloak && return _cloak_record(variant;kwargs...)
    throw(ArgumentError("unknown spatial experiment"))
end
