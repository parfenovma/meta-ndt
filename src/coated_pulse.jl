export CoatedPulseExperiment, CoatedPulseResult, coated_pulse_comparison,
       coated_trace, coated_incident_pressure, coated_pulse_drive, coated_pulse_profile

"""Direct time-domain concentric acoustic–elastic FEM, SI, plane strain.

Water is represented as scattered velocity potential plus an exact finite
incident plane pulse. Air is a compressible FEM domain, not an empty cavity.
Kelvin–Voigt shear viscosity matches the harmonic loss tangent at the FIXED
`loss_reference_frequency`, not at each excitation frequency. No FFT synthesis.
"""
Base.@kwdef struct CoatedPulseExperiment
    model::CoatedCylinderExperiment = CoatedCylinderExperiment()
    pulse_cycles::Float64 = 3.0
    source_plane::Float64 = -.045
    duration::Float64 = pulse_cycles/model.frequency+160e-6
    time_step::Float64 = min(.5e-6,1/(80model.frequency))
    frame_count::Int = 141
    water_mesh_size::Float64 = min(.0015,model.water_speed/(32model.frequency))
    water_interface_divisions::Int = 16
    air_elements::Int = 48
    retained_radius::Float64 = .040
    outer_radius::Float64 = retained_radius+model.water_speed*duration
    loss_reference_frequency::Float64 = 10e3
end

function validate_experiment(s::CoatedPulseExperiment)
    validate_experiment(s.model)
    for k in (:pulse_cycles,:duration,:time_step,:water_mesh_size,:retained_radius,:outer_radius,:loss_reference_frequency)
        v=getfield(s,k); isfinite(v)&&v>0 || throw(ArgumentError("$k must be positive and finite"))
    end
    s.pulse_cycles>=2 || throw(ArgumentError("use at least two complete pulse cycles"))
    isinteger(s.pulse_cycles) || throw(ArgumentError("use an integer number of pulse cycles"))
    s.duration>s.pulse_cycles/s.model.frequency || throw(ArgumentError("time window must exceed source duration"))
    s.time_step*s.model.frequency<=.02 || throw(ArgumentError("at least fifty steps per carrier period required"))
    s.water_mesh_size*s.model.frequency<=s.model.water_speed/24 || throw(ArgumentError("need 24 water elements per carrier wavelength"))
    s.air_elements>=24 && s.frame_count>=2 && s.water_interface_divisions>=12 || throw(ArgumentError("insufficient air/interface/frame sampling"))
    s.retained_radius>s.model.outer_radius+s.model.coating_thickness || throw(ArgumentError("retained water domain too small"))
    isfinite(s.source_plane) && s.source_plane < -s.retained_radius ||
        throw(ArgumentError("source plane must precede retained water so every shown point starts at rest"))
    2(s.outer_radius-s.retained_radius)/s.model.water_speed>s.duration ||
        throw(ArgumentError("outer boundary return would enter the retained observation window"))
    s
end

function remake_experiment(s::CoatedPulseExperiment;kwargs...)
    values=NamedTuple{fieldnames(CoatedPulseExperiment)}(Tuple(getfield(s,k) for k in fieldnames(CoatedPulseExperiment)))
    validate_experiment(CoatedPulseExperiment(;merge(values,(;kwargs...))...))
end

function experiment_toml(s::CoatedPulseExperiment)
    validate_experiment(s)
    d=Dict{String,Any}(String(k)=>getfield(s,k) for k in fieldnames(CoatedPulseExperiment) if k!=:model)
    d["model"]=TOML.parse(experiment_toml(s.model))
    d["physics"]="coated_pulse";d["schema_version"]=1
    io=IOBuffer();TOML.print(io,d;sorted=true);String(take!(io))
end
save_experiment(path::AbstractString,s::CoatedPulseExperiment)=write(path,experiment_toml(s))

struct CoatedPulseResult
    experiment::CoatedPulseExperiment
    times::Vector{Float64}
    frame_steps::Vector{Int}
    water_nodes::Vector{Float64}
    # Radial nodes × EVERY integration step × cosine angular modes. Pa.
    scattered_pressure::Array{Float32,3}
    energy::Vector{Float64}
    solid_energy::Vector{Float64}
    air_energy::Vector{Float64}
    input_work::Vector{Float64}
    dissipated_energy::Vector{Float64}
    diagnostics::Dict{String,Float64}
end

function coated_pulse_drive(s::CoatedPulseExperiment,t::Real)
    T=s.pulse_cycles/s.model.frequency
    0<t<T ? sin(2π*s.model.frequency*t)*sinpi(t/T)^2 : 0.0
end

function coated_incident_pressure(s::CoatedPulseExperiment,point,t::Real)
    m=s.model; x=point[1]*cos(m.incident_angle)+point[2]*sin(m.incident_angle)
    m.incident_pressure*coated_pulse_drive(s,t-(x-s.source_plane)/m.water_speed)
end

# Scalar P1 radial velocity-potential forms: regular at r=0. For n>0 impose
# phi(0)=0; mode zero has its natural symmetry condition and is NOT clamped.
function _radial_acoustic_forms(nodes,n,ρ,c,F;regular=false)
    I=Int[];J=Int[];mv=Float64[];kv=Float64[]
    for e in 1:length(nodes)-1
        lo,hi=nodes[e:e+1]; h=hi-lo; me=zeros(2,2);ke=zeros(2,2)
        for (ξ,w) in zip(_RADIAL_GAUSS_X,_RADIAL_GAUSS_W)
            r=(lo+hi)/2+h*ξ/2; N=[(1-ξ)/2,(1+ξ)/2]; dN=[-1/h,1/h]
            me .+= (w*h*r/2*ρ*F^2/c^2)*(N*transpose(N))
            ke .+= (w*h/2*ρ*F^2)*(r*dN*transpose(dN)+(n^2/r)*(N*transpose(N)))
        end
        for j in 1:2,i in 1:2
            push!(I,e+i-1);push!(J,e+j-1);push!(mv,me[i,j]);push!(kv,ke[i,j])
        end
    end
    count=length(nodes); active=regular&&n>0 ? (2:count) : (1:count)
    M=sparse(I,J,mv,count,count);K=sparse(I,J,kv,count,count)
    (;M=M[active,active],K=K[active,active])
end

function _coated_pulse_system(s,n,water_nodes)
    m=s.model; solid_nodes,els,mats=_radial_mesh(m)
    a=_radial_elastic_forms(m,n,solid_nodes,els,mats)
    # Eliminate the static volumetric stress, preserving precisely the same
    # mixed P2/P1 form as the harmonic backend (not a locking-prone P1 solid).
    L=last(solid_nodes);F=m.water_speed*L
    Ks=sparse(Symmetric(Matrix(a.Kdev+a.B*(a.Q\Matrix(transpose(a.B))))))*L^2
    Ms=a.M*L^2;Cs=a.Kcoating*(L^2*m.coating_shear_loss/(2π*s.loss_reference_frequency))
    air_nodes=collect(range(0,m.inner_radius;length=s.air_elements+1))
    w=_radial_acoustic_forms(water_nodes,n,m.water_density,m.water_speed,F)
    air=_radial_acoustic_forms(air_nodes,n,m.air_density,m.air_speed,F;regular=true)
    nw,ns,na=size(w.M,1),size(Ms,1),size(air.M,1);N=nw+ns+na
    M=blockdiag(w.M,Ms,air.M);K=blockdiag(w.K,Ks,air.K)
    Cd=blockdiag(spzeros(nw,nw),Cs,spzeros(na,na))
    Cd[nw,nw]=m.water_density*F^2*s.outer_radius/m.water_speed
    Ci=spzeros(N,N)
    # Fluid normal velocity ↔ solid pressure traction: skew interface damping.
    outer=nw+a.nv;inner=nw+1;wa=1;aa=N
    γw=m.water_density*F*L*L;γa=m.air_density*F*L*m.inner_radius
    Ci[wa,outer]=γw;Ci[outer,wa]=-γw
    Ci[aa,inner]=-γa;Ci[inner,aa]=γa
    (;M,K,Cd,Ci,C=Cd+Ci,nw,ns,na,F,L,outer,inner,solid_nodes)
end

function _coated_pulse_forcing(s,R,times)
    # Fourier projection of the exact delayed plane pulse on the OUTER solid
    # interface, including both incident pressure and incident normal velocity.
    nθ=max(128,8(s.model.modes+1));θ=(0:nθ-1).*(2π/nθ);cs=cos.(θ)
    W=[(n==0 ? 1.0 : 2.0)*cos(n*t)/nθ for n in 0:s.model.modes,t in θ]
    p=zeros(s.model.modes+1,length(times));v=similar(p)
    for (j,t) in enumerate(times)
        values=[s.model.incident_pressure*coated_pulse_drive(s,t-(R*x-s.source_plane)/s.model.water_speed) for x in cs]
        p[:,j]=W*values;v[:,j]=W*(values.*cs)/s.model.water_speed
    end
    (;p,v)
end

"""Direct Newmark integration of coupled water/steel/coating/air radial FEM.

Full-step pressure histories are retained only up to `retained_radius`; the
outer domain is large enough that no reflected characteristic can return to
this region during the requested window. No harmonic solutions are used.
"""
function simulate_pulse(s::CoatedPulseExperiment)
    s=remake_experiment(s);m=s.model;started=time_ns()
    ratio=s.duration/s.time_step;N=ceil(Int,ratio-8eps(ratio));dt=s.duration/N
    times=collect(range(0,s.duration;length=N+1))
    frames=unique(round.(Int,range(1,N+1;length=min(s.frame_count,N+1))))
    energy=zeros(N+1);es=zeros(N+1);ea=zeros(N+1);work=zeros(N+1);loss=zeros(N+1)
    if m.scenario==:empty
        return CoatedPulseResult(s,times,frames,Float64[],zeros(Float32,0,N+1,m.modes+1),
            energy,es,ea,work,loss,Dict("energy_balance_relative"=>0.,"interface_skew_relative"=>0.,
            "peak_displacement_m"=>0.,"peak_air_pressure_Pa"=>0.,"time_step_s"=>dt,"total_seconds"=>(time_ns()-started)/1e9))
    end
    R=m.outer_radius+(m.scenario==:coated ? m.coating_thickness : 0.)
    water=Float64[R]
    while last(water)<s.outer_radius
        push!(water,min(s.outer_radius,last(water)+min(s.water_mesh_size,last(water)/s.water_interface_divisions)))
    end
    kept=min(length(water),searchsortedfirst(water,s.retained_radius))
    history=zeros(Float32,kept,N+1,m.modes+1);source=_coated_pulse_forcing(s,R,times)
    peak=0.;peak_air=0.;skew=0.
    for n in 0:m.modes
        a=_coated_pulse_system(s,n,water);count=size(a.M,1)
        solver=lu(a.M+dt/2*a.C+dt^2/4*a.K)
        d=zeros(count);v=copy(d);acc=copy(d);fprev=copy(d)
        solid=(a.nw+1):(a.nw+a.ns);air=(last(solid)+1):count
        Ms=a.M[solid,solid];Ks=a.K[solid,solid];Ma=a.M[air,air];Ka=a.K[air,air]
        weight=n==0 ? 2π : π
        skew=max(skew,norm(a.Ci+transpose(a.Ci))/max(norm(a.Ci),eps()))
        for step in 1:N+1
            if step>1
                dp=d+dt*v+dt^2/4*acc;vp=v+dt/2*acc
                f=zeros(count);f[1]=a.F*R*source.v[n+1,step];f[a.outer]=-a.L*R*source.p[n+1,step]
                anew=solver\(f-a.C*vp-a.K*dp)
                dn=dp+dt^2/4*anew;vn=vp+dt/2*anew;vm=(v+vn)/2
                # Sum modal energy increments; cumulative sums follow below.
                work[step]+=weight*dt*dot(vm,(fprev+f)/2)
                loss[step]+=weight*dt*dot(vm,a.Cd*vm)
                d=dn;v=vn;acc=anew;fprev=f
            end
            energy[step]+=weight*.5(dot(v,a.M*v)+dot(d,a.K*d))
            es[step]+=weight*.5(dot(v[solid],Ms*v[solid])+dot(d[solid],Ks*d[solid]))
            ea[step]+=weight*.5(dot(v[air],Ma*v[air])+dot(d[air],Ka*d[air]))
            history[:,step,n+1]=-m.water_density*a.F*v[1:kept]
            peak=max(peak,a.L*maximum(abs,view(d,solid)))
            peak_air=max(peak_air,m.air_density*a.F*maximum(abs,view(v,air)))
        end
    end
    cumsum!(work,work);cumsum!(loss,loss)
    diag=Dict("energy_balance_relative"=>maximum(abs.(energy+loss-work))/max(maximum(abs,work),floatmin(Float64)),
        "interface_skew_relative"=>skew,"peak_displacement_m"=>peak,"peak_air_pressure_Pa"=>peak_air,
        "time_step_s"=>dt,"total_seconds"=>(time_ns()-started)/1e9,"retained_history_bytes"=>Float64(sizeof(history)),
        "outer_return_lower_bound_s"=>2(s.outer_radius-s.retained_radius)/m.water_speed,
        "kelvin_voigt_time_s"=>m.coating_shear_loss/(2π*s.loss_reference_frequency))
    CoatedPulseResult(s,times,frames,water[1:kept],history,energy,es,ea,work,loss,diag)
end

"""Ideal water microphone, Pa, at every integration step by default.

Exact known incident pulse plus P1 interpolation of calculated scattered
pressure modes; no display-pixel lookup, remeshing, phasor or repeat solve.
`indices` can select display frames without changing the underlying history.
"""
function coated_trace(r::CoatedPulseResult,point;indices=eachindex(r.times))
    length(point)==2&&all(isfinite,point) || throw(ArgumentError("invalid water point"))
    s=r.experiment;m=s.model;radius=hypot(point...)
    R=m.outer_radius+(m.scenario==:coated ? m.coating_thickness : 0.)
    radius<=s.retained_radius || throw(ArgumentError("point outside retained pressure history"))
    m.scenario==:empty||radius>=R || throw(ArgumentError("microphone must be in exterior water"))
    ids=collect(indices);all(i->i in eachindex(r.times),ids) || throw(ArgumentError("invalid time indices"))
    scattered=zeros(length(ids))
    if m.scenario!=:empty
        nodes=r.water_nodes;e=clamp(searchsortedlast(nodes,radius),1,length(nodes)-1)
        w=(radius-nodes[e])/(nodes[e+1]-nodes[e]);θ=atan(point[2],point[1])-m.incident_angle
        for n in 0:m.modes,j in eachindex(ids)
            scattered[j]+=cos(n*θ)*((1-w)*r.scattered_pressure[e,ids[j],n+1]+w*r.scattered_pressure[e+1,ids[j],n+1])
        end
    end
    times=r.times[ids];incident=[coated_incident_pressure(s,point,t) for t in times]
    (;times,pressure=incident+scattered,scattered)
end

coated_pulse_comparison(s::CoatedPulseExperiment=CoatedPulseExperiment())=
    [simulate_pulse(remake_experiment(s;model=remake_experiment(s.model;scenario=c))) for c in (:empty,:bare,:coated)]

"""Finite-window angular profile ∫p_scattered²dt / ∫p_incident(origin)²dt.

This is a microphone-pressure observable at a specified contour, NOT radiated
energy or stationary far-field scattering width. No independent peak scaling.
"""
function coated_pulse_profile(r::CoatedPulseResult,angles;radius=r.experiment.model.observation_radius)
    isfinite(radius)&&0<radius<=r.experiment.retained_radius || throw(ArgumentError("invalid observation radius"))
    all(isfinite,angles) || throw(ArgumentError("invalid observation angles"))
    window=(0.,last(r.times))
    normin=pressure_integral(r.times,[coated_incident_pressure(r.experiment,(0.,0.),t) for t in r.times],window)
    normin>0 || throw(ArgumentError("incident pulse has not reached reference point within this time window"))
    [pressure_integral(r.times,coated_trace(r,(radius*cos(θ),radius*sin(θ))).scattered,window)/normin for θ in angles]
end
