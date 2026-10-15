using SparseArrays
using SpecialFunctions: besselj, bessely, hankelh1

export CoatedCylinderExperiment, CoatedCylinderResult, CoatedCylinderSweep,
       coated_comparison, coated_frequency_sweep, coated_sweep_csv,
       coated_pressure, coated_displacement, coated_cavity_resonances

"""Concentric, bonded steel/silicone tube in water; SI, plane strain, exp(-iωt).
This fixed-geometry backend uses angular Fourier modes and mixed radial FEM.
The sealed cavity contains acoustic air, not water and not a fixed boundary.
"""
Base.@kwdef struct CoatedCylinderExperiment
    scenario::Symbol = :coated
    outer_radius::Float64 = 0.004
    inner_radius::Float64 = 0.004sqrt((7850-1000)/(7850-1.2))
    coating_thickness::Float64 = 0.001
    water_density::Float64 = 1000.0
    water_speed::Float64 = 1480.0
    air_density::Float64 = 1.2
    air_speed::Float64 = 343.0
    steel_density::Float64 = 7850.0
    steel_young::Float64 = 210e9
    steel_poisson::Float64 = 0.3
    coating_density::Float64 = 1020.0
    coating_pressure_speed::Float64 = 980.0
    coating_shear_speed::Float64 = 30.0
    coating_shear_loss::Float64 = 0.15
    frequency::Float64 = 10e3
    incident_pressure::Float64 = 1.0
    incident_angle::Float64 = 0.0
    steel_elements::Int = 12
    coating_elements::Int = 48
    modes::Int = 12
    observation_radius::Float64 = 0.030
    angular_samples::Int = 181
end

function validate_experiment(s::CoatedCylinderExperiment)
    s.scenario in (:empty,:bare,:coated) || throw(ArgumentError("unknown coated-cylinder case"))
    for k in (:outer_radius,:inner_radius,:coating_thickness,:water_density,:water_speed,
              :air_density,:air_speed,:steel_density,:steel_young,:coating_density,
              :coating_pressure_speed,:coating_shear_speed,:frequency,
              :incident_pressure,:observation_radius)
        v=getfield(s,k)
        isfinite(v) && v > 0 || throw(ArgumentError("$k must be finite and positive"))
    end
    s.inner_radius < s.outer_radius || throw(ArgumentError("tube must have a positive wall thickness"))
    0 <= s.steel_poisson < 0.5 || throw(ArgumentError("steel Poisson ratio must be in [0,0.5)"))
    isfinite(s.coating_shear_loss) && s.coating_shear_loss >= 0 || throw(ArgumentError("negative/nonfinite shear loss"))
    s.coating_pressure_speed^2 > 4s.coating_shear_speed^2/3 || throw(ArgumentError("nonpositive coating bulk modulus"))
    isfinite(s.incident_angle) || throw(ArgumentError("nonfinite incident angle"))
    s.steel_elements >= 2 && s.coating_elements >= 2 || throw(ArgumentError("at least two radial elements per layer"))
    s.modes >= 3 || throw(ArgumentError("at least modes 0 through 3 required"))
    s.angular_samples >= 4s.modes+1 || throw(ArgumentError("angular sampling cannot resolve retained modes"))
    s.observation_radius > s.outer_radius+s.coating_thickness || throw(ArgumentError("observation contour must surround entire coating"))
    s.coating_shear_speed*s.coating_elements/(s.frequency*s.coating_thickness) >= 10 ||
        throw(ArgumentError("need at least ten radial elements per coating shear wavelength"))
    2π*s.frequency*(s.outer_radius+s.coating_thickness)/s.water_speed < s.modes-2 ||
        throw(ArgumentError("too few angular modes for this frequency/size"))
    s
end

function remake_experiment(s::CoatedCylinderExperiment;kwargs...)
    parameters=NamedTuple{fieldnames(CoatedCylinderExperiment)}(Tuple(getfield(s,k) for k in fieldnames(CoatedCylinderExperiment)))
    validate_experiment(CoatedCylinderExperiment(;merge(parameters,(;kwargs...))...))
end

function experiment_toml(s::CoatedCylinderExperiment)
    validate_experiment(s)
    d=Dict(String(k)=>getfield(s,k) for k in fieldnames(CoatedCylinderExperiment))
    d["physics"]="coated_cylinder"; d["schema_version"]=1; d["scenario"]=String(s.scenario)
    io=IOBuffer(); TOML.print(io,d;sorted=true); String(take!(io))
end
save_experiment(path::AbstractString,s::CoatedCylinderExperiment)=write(path,experiment_toml(s))

struct CoatedCylinderResult
    experiment::CoatedCylinderExperiment
    coefficients::Vector{ComplexF64}
    radial_nodes::Vector{Float64}
    displacement_radial::Matrix{ComplexF64}
    displacement_tangential::Matrix{ComplexF64}
    angles::Vector{Float64}
    scattered_pressure::Vector{ComplexF64}
    far_field_amplitude::Vector{ComplexF64}
    scattering_width::Float64
    diagnostics::Dict{String,Float64}
end

# Derivatives with respect to the Bessel argument; recurrence avoids differences.
_bessel_derivative(f,n,z)=(f(n-1,z)-f(n+1,z))/2
_modal_weight(n)=n == 0 ? 1 : 2

function _coated_materials(s)
    μs=s.steel_young/(2(1+s.steel_poisson))
    λs=s.steel_young*s.steel_poisson/((1+s.steel_poisson)*(1-2s.steel_poisson))
    μc=s.coating_density*s.coating_shear_speed^2*(1-im*s.coating_shear_loss)
    Kc=s.coating_density*(s.coating_pressure_speed^2-4s.coating_shear_speed^2/3)
    (steel=(s.steel_density,complex(λs+2μs/3),complex(μs)),
     coating=(s.coating_density,complex(Kc),μc))
end

function _radial_mesh(s)
    steel=collect(range(s.inner_radius,s.outer_radius;length=2s.steel_elements+1))
    if s.scenario == :bare
        return steel,[(i,i+1,i+2) for i in 1:2:length(steel)-2],fill(:steel,s.steel_elements)
    end
    coating=collect(range(s.outer_radius,s.outer_radius+s.coating_thickness;length=2s.coating_elements+1))
    nodes=vcat(steel,coating[2:end]); count=s.steel_elements+s.coating_elements
    nodes,[(i,i+1,i+2) for i in 1:2:2count-1],vcat(fill(:steel,s.steel_elements),fill(:coating,s.coating_elements))
end

# Five-point quadrature on [-1,1]. P2 displacement / P1 volumetric stress,
# with separate stress DOFs on the two sides of the bonded material interface.
const _RADIAL_GAUSS_X=(-0.906179845938664,-0.538469310105683,0.0,0.538469310105683,0.906179845938664)
const _RADIAL_GAUSS_W=(0.236926885056189,0.478628670499366,0.568888888888889,0.478628670499366,0.236926885056189)

include("radial_elasticity.jl")

function _radial_system(s,n,nodes,elements,materials)
    a=_radial_elastic_forms(s,n,nodes,elements,materials); nv=a.nv
    nu,np=size(a.B); ω=2π*s.frequency
    K=a.Kdev-im*s.coating_shear_loss*a.Kcoating-ω^2*a.M
    A=[K a.B spzeros(nu,1); transpose(a.B) -a.Q spzeros(np,1);
        spzeros(1,nu) spzeros(1,np) spzeros(1,1)]
    rhs=zeros(ComplexF64,size(A,1))
    b,R=first(nodes),last(nodes); k=ω/s.water_speed; ka=ω/s.air_speed
    # Keep cavity amplitude explicitly: dividing by Jn′ would be singular at
    # rigid-cavity Neumann frequencies, which lie inside the displayed band.
    Jair=besselj(n,ka*b); dJair=ka*_bessel_derivative(besselj,n,ka*b)
    H=hankelh1(n,k*R); dH=k*_bessel_derivative(hankelh1,n,k*R)
    incident=_modal_weight(n)*im^n*s.incident_pressure
    p0=incident*besselj(n,k*R); dp0=incident*k*_bessel_derivative(besselj,n,k*R)
    waterZ=s.water_density*ω^2*H/dH; g=p0-H*dp0/dH
    elasticA=copy(A)
    A[1,end]=-b*Jair
    A[end,1]=-b; A[end,end]=b*dJair/(s.air_density*ω^2)
    A[nv,nv]+=R*waterZ; rhs[nv]=-R*g
    (;A,rhs,dp0,dH,nv,elasticA)
end

"""Full tensor, coupled plane-strain solve, with an exact outgoing water condition.
No outer truncation boundary/PML is involved. Includes compressible cavity air.
"""
function simulate(s::CoatedCylinderExperiment)
    validate_experiment(s); started=time_ns()
    angles=collect(range(0,2π;length=s.angular_samples))
    coefficients=zeros(ComplexF64,s.modes+1); residual=0.0; absorbed=0.0
    if s.scenario == :empty
        nodes=Float64[]; ur=zeros(ComplexF64,0,s.modes+1); ut=copy(ur); dofs=0
    else
        nodes,elements,materials=_radial_mesh(s)
        ur=zeros(ComplexF64,length(nodes),s.modes+1); ut=copy(ur); dofs=0
        for n in 0:s.modes
            a=_radial_system(s,n,nodes,elements,materials)
            d=a.A\a.rhs
            # Rowwise scaled residual: raw mixed-system rows have different units.
            res=abs.(a.A*d-a.rhs); scale=abs.(a.A)*abs.(d)+abs.(a.rhs)
            residual=max(residual,maximum(res./max.(scale,eps(Float64))))
            absorbed-=(2π*s.frequency/2)*(n == 0 ? 2π : π)*imag(dot(d,a.elasticA*d))
            ur[:,n+1]=d[1:a.nv]
            n > 0 && (ut[:,n+1]=d[a.nv+1:2a.nv])
            coefficients[n+1]=(s.water_density*(2π*s.frequency)^2*ur[end,n+1]-a.dp0)/a.dH
            dofs=max(dofs,length(d))
        end
    end
    k=2π*s.frequency/s.water_speed
    scattered=[sum(coefficients[n+1]*hankelh1(n,k*s.observation_radius)*cos(n*(θ-s.incident_angle)) for n in 0:s.modes) for θ in angles]
    far=[sqrt(2/(π*k))*sum(coefficients[n+1]*exp(-im*n*π/2)*cos(n*(θ-s.incident_angle)) for n in 0:s.modes) for θ in angles]
    width=(2/k)*(2abs2(coefficients[1])+sum(abs2,coefficients[2:end]))/s.incident_pressure^2
    extinction=-(4/k)*real(sum(coefficients[n+1]*exp(-im*n*π/2) for n in 0:s.modes))/s.incident_pressure
    absorbed_width=absorbed/(s.incident_pressure^2/(2s.water_density*s.water_speed))
    tail=(2/k)*sum(abs2,coefficients[end-1:end])/s.incident_pressure^2
    diag=Dict("total_seconds"=>(time_ns()-started)/1e9,"dofs_per_mode"=>Float64(dofs),
        "scaled_residual"=>residual,"modal_tail_fraction"=>tail/max(width,eps()),
        "extinction_width_m"=>extinction,"absorption_width_m"=>absorbed_width,
        "power_balance_relative"=>abs(extinction-width-absorbed_width)/max(extinction,width,eps()),
        "ka"=>k*s.outer_radius,"coating_shear_wavelength_elements"=>s.coating_shear_speed*s.coating_elements/(s.frequency*s.coating_thickness),
        "steel_average_density"=>s.steel_density*(1-(s.inner_radius/s.outer_radius)^2)+s.air_density*(s.inner_radius/s.outer_radius)^2)
    CoatedCylinderResult(s,coefficients,nodes,ur,ut,angles,scattered,far,width,diag)
end

"""Complex water pressure, Pa. `scattered=true` removes the known incident wave.
Points in the air/steel/coating are rejected; they are not water-pressure pixels.
"""
function coated_pressure(r::CoatedCylinderResult,point;scattered::Bool=false)
    length(point) == 2 && all(isfinite,point) || throw(ArgumentError("invalid water point"))
    s=r.experiment; radius=hypot(point...)
    R=s.outer_radius+(s.scenario == :coated ? s.coating_thickness : 0)
    s.scenario == :empty || radius >= R || throw(ArgumentError("pressure point is not in exterior water"))
    k=2π*s.frequency/s.water_speed; θ=atan(point[2],point[1])-s.incident_angle
    p=s.scenario == :empty ? 0.0im : sum(r.coefficients[n+1]*hankelh1(n,k*radius)*cos(n*θ) for n in 0:s.modes)
    scattered ? p : p+s.incident_pressure*exp(im*k*(point[1]*cos(s.incident_angle)+point[2]*sin(s.incident_angle)))
end

"""Interpolate the calculated P2 displacement at a point in either solid layer."""
function coated_displacement(r::CoatedCylinderResult,point)
    length(point) == 2 && all(isfinite,point) || throw(ArgumentError("invalid solid point"))
    radius=hypot(point...); nodes=r.radial_nodes
    !isempty(nodes) && first(nodes) <= radius <= last(nodes) || throw(ArgumentError("point is not in the solid"))
    e=clamp(cld(searchsortedlast(nodes,radius),2),1,(length(nodes)-1)÷2)
    ids=(2*e-1):(2*e+1); lo,hi=nodes[first(ids)],nodes[last(ids)]; ξ=2(radius-lo)/(hi-lo)-1
    N=[ξ*(ξ-1)/2,1-ξ^2,ξ*(ξ+1)/2]; θ=atan(point[2],point[1])-r.experiment.incident_angle
    ur=sum(sum(N.*r.displacement_radial[ids,n+1])*cos(n*θ) for n in 0:r.experiment.modes)
    ut=sum(sum(N.*r.displacement_tangential[ids,n+1])*sin(n*θ) for n in 0:r.experiment.modes)
    angle=atan(point[2],point[1]); (ur*cos(angle)-ut*sin(angle),ur*sin(angle)+ut*cos(angle))
end

coated_comparison(s::CoatedCylinderExperiment=CoatedCylinderExperiment())=
    [simulate(remake_experiment(s;scenario=c)) for c in (:empty,:bare,:coated)]

struct CoatedCylinderSweep
    experiments::Vector{CoatedCylinderExperiment}
    frequencies::Vector{Float64}
    bare_width::Vector{Float64}
    coated_width::Vector{Float64}
    ratio::Vector{Float64}
    cavity_frequencies::Vector{Float64}
end

"""Rigid-wall air-cavity Neumann frequencies within a band, Hz.
These are sampling seeds, NOT exact resonances of the coupled elastic object.
"""
function coated_cavity_resonances(s::CoatedCylinderExperiment,lower::Real,upper::Real)
    validate_experiment(s)
    isfinite(lower) && isfinite(upper) && 0 < lower < upper || throw(ArgumentError("invalid cavity frequency band"))
    scale=s.air_speed/(2π*s.inner_radius); roots=Float64[]
    for n in 0:s.modes
        previous=lower/scale; value=_bessel_derivative(besselj,n,previous)
        for x in range(previous,upper/scale;length=max(2,ceil(Int,(upper-lower)/scale/.1)+1))[2:end]
            current=_bessel_derivative(besselj,n,x)
            if signbit(current) != signbit(value)
                lo,hi=previous,x
                for _ in 1:48
                    mid=(lo+hi)/2
                    if signbit(_bessel_derivative(besselj,n,mid)) == signbit(value)
                        lo=mid
                    else
                        hi=mid
                    end
                end
                push!(roots,(lo+hi)/2*scale)
            end
            previous=x; value=current
        end
    end
    sort!(unique!(roots))
end

function coated_frequency_sweep(s::CoatedCylinderExperiment,frequencies;resolve_cavity::Bool=false)
    f=Float64.(collect(frequencies))
    !isempty(f) && all(isfinite,f) && all(f.>0) && all(diff(f).>0) || throw(ArgumentError("positive increasing frequencies required"))
    seeds=length(f) > 1 ? coated_cavity_resonances(s,first(f),last(f)) : Float64[]
    if resolve_cavity
        # Fixed targeted sampling, not a claim that all resonances are resolved.
        # ~2.5-Hz-wide default cavity peaks need much finer than a 5-kHz grid.
        enriched=[seed+δ for seed in seeds for δ in vcat(collect(-30:.5:30),collect(-1000:100:1000))
            if first(f) <= seed+δ <= last(f)]
        f=sort!(unique!(vcat(f,enriched)))
    end
    specs=[remake_experiment(s;frequency=x) for x in f] # validate all before solving
    bare=Float64[]; coated=Float64[]
    for spec in specs
        push!(bare,simulate(remake_experiment(spec;scenario=:bare)).scattering_width)
        push!(coated,simulate(remake_experiment(spec;scenario=:coated)).scattering_width)
    end
    all(bare.>0) || throw(ArgumentError("bare scattering reference is zero"))
    CoatedCylinderSweep(specs,f,bare,coated,coated./bare,seeds)
end

function coated_sweep_csv(r::CoatedCylinderSweep)
    io=IOBuffer(); println(io,"frequency_Hz,bare_scattering_width_m,coated_scattering_width_m,coated_over_bare,ratio_dB")
    for i in eachindex(r.frequencies)
        println(io,join((r.frequencies[i],r.bare_width[i],r.coated_width[i],r.ratio[i],10log10(max(r.ratio[i],eps()))),','))
    end
    String(take!(io))
end
