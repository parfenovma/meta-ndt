export SolidLensExperiment, SolidLensResult, LensDesign, solid_lens_outlines,
       solid_lens_comparison, design_stg_lens, simulate_pulse, solid_profile_csv,
       lens_design_csv, solid_wave_speeds, SolidPulseResult, solid_examples, solid_preset,
       solid_pulse_drive

"""Connected solid plate with a wall of STG ligaments and traction-free voids. SI, plane strain."""
Base.@kwdef struct SolidLensExperiment
    scenario::Symbol = :stg
    bricks::Vector{SinusoidalBrick} = [SinusoidalBrick() for _ in 1:13]
    missing_cells::Vector{Int} = Int[]
    width::Float64 = .180
    height::Float64 = .112
    density::Float64 = 2700.0
    young_modulus::Float64 = 70e9
    poisson_ratio::Float64 = .33
    frequency::Float64 = 140e3
    polarization::Symbol = :shear
    source_x::Float64 = -.060
    source_width::Float64 = .002
    source_strength::Float64 = 1.0
    focus_distance::Float64 = .040
    mesh_size::Float64 = .0012
    near_mesh_size::Float64 = .0005
    absorb_sides::Bool = true
    pulse_cycles::Float64 = 6.0
    duration::Float64 = 140e-6
    time_step::Float64 = .1e-6
    frame_count::Int = 81
end

function solid_wave_speeds(s::SolidLensExperiment)
    μ=s.young_modulus/(2*(1+s.poisson_ratio))
    λ=s.young_modulus*s.poisson_ratio/((1+s.poisson_ratio)*(1-2s.poisson_ratio))
    (sqrt((λ+2μ)/s.density),sqrt(μ/s.density))
end

function validate_experiment(s::SolidLensExperiment)
    s.scenario in (:solid,:straight,:stg) || throw(ArgumentError("unknown solid lens scenario"))
    s.polarization in (:shear,:pressure) || throw(ArgumentError("choose shear or pressure incidence"))
    all(x->isfinite(x) && x>0,(s.width,s.height,s.density,s.young_modulus,s.frequency,
        s.source_width,s.focus_distance,s.mesh_size,s.near_mesh_size,s.pulse_cycles,s.duration,s.time_step)) ||
        throw(ArgumentError("solid lens material/length/time values must be finite and positive"))
    isfinite(s.source_strength) && isfinite(s.source_x) || throw(ArgumentError("nonfinite solid source"))
    -1 < s.poisson_ratio < .5 || throw(ArgumentError("invalid Poisson ratio"))
    !isempty(s.bricks) && isodd(length(s.bricks)) || throw(ArgumentError("use an odd, nonempty wall"))
    length(unique(s.missing_cells)) == length(s.missing_cells) &&
        all(i->i in eachindex(s.bricks),s.missing_cells) && length(s.missing_cells)<length(s.bricks) ||
        throw(ArgumentError("invalid removed wall cells"))
    foreach(validate_shape,s.bricks)
    b=first(s.bricks)
    all(c->c.length==b.length && c.height==b.height,s.bricks) ||
        throw(ArgumentError("wall cells must have common length and pitch/height"))
    length(s.bricks)*b.height <= s.height || throw(ArgumentError("wall aperture exceeds plate height"))
    -s.width/2+3s.source_width < s.source_x < -b.length/2-3s.source_width ||
        throw(ArgumentError("source must be wholly inside the left bulk"))
    b.length/2+s.focus_distance+3s.mesh_size < s.width/2 || throw(ArgumentError("focus too near outlet"))
    s.near_mesh_size <= s.mesh_size && s.near_mesh_size <= minimum(c.height-c.amplitude for c in s.bricks)/3 ||
        throw(ArgumentError("need at least three elements through narrowest ligament"))
    s.mesh_size*s.frequency <= last(solid_wave_speeds(s))/12 ||
        throw(ArgumentError("need twelve nominal elements per shear carrier wavelength"))
    s.frame_count >= 2 && s.time_step*s.frequency <= .04 &&
        s.duration > s.pulse_cycles/s.frequency || throw(ArgumentError("invalid pulse time preset"))
    s
end

function remake_experiment(s::SolidLensExperiment;kwargs...)
    values=NamedTuple{fieldnames(SolidLensExperiment)}(Tuple(getfield(s,k) for k in fieldnames(SolidLensExperiment)))
    p=merge(values,(;kwargs...)); p=merge(p,(bricks=copy(p.bricks),missing_cells=copy(p.missing_cells)))
    validate_experiment(SolidLensExperiment(;p...))
end

function solid_lens_outlines(s::SolidLensExperiment)
    validate_experiment(s)
    N=length(s.bricks); pitch=first(s.bricks).height
    [[(p[1],p[2]+(j-(N+1)/2)*pitch) for p in brick_polygon(
        s.scenario == :straight ? RectangularBrick(length=b.length,height=specimen_area(b)/b.length) : b)]
        for (j,b) in enumerate(s.bricks) if !(j in s.missing_cells)]
end

struct SolidLensResult
    experiment::SolidLensExperiment
    nodes::Matrix{Float64}
    triangles::Matrix{Int}
    displacement_x::Vector{ComplexF64}
    displacement_y::Vector{ComplexF64}
    profile_y::Vector{Float64}
    profile_kinetic::Vector{Float64}
    scan_x::Vector{Float64}
    scan_fwhm::Vector{Float64}
    scan_peak::Vector{Float64}
    diagnostics::Dict{String,Float64}
end

struct LensDesign
    experiment::SolidLensExperiment
    candidates::Vector{SinusoidalBrick}
    transmission::Vector{ComplexF64}
    selected::Vector{Int}
    phase_error::Vector{Float64}
end

struct SolidPulseResult
    experiment::SolidLensExperiment
    nodes::Matrix{Float64}
    triangles::Matrix{Int}
    frame_times::Vector{Float64}
    displacement_x::Matrix{Float64}
    displacement_y::Matrix{Float64}
    times::Vector{Float64}
    energy::Vector{Float64}
    input_work::Vector{Float64}
    outgoing_energy::Vector{Float64}
    profile_y::Vector{Float64}
    profile_integrated_kinetic::Vector{Float64}
    diagnostics::Dict{String,Float64}
    # Optional nodal displacement at EVERY integration step. Solver arithmetic
    # stays Float64; Float32 retention avoids making interactive probes depend
    # on the much sparser display-frame sampling.
    history_x::Matrix{Float32}
    history_y::Matrix{Float32}
end

function design_stg_lens end
"""Transient P1 elastodynamics, zero initial state, finite Hann-windowed force.

`retain_history=true` retains Float32 nodal displacement at every time step for
subsequent arbitrary-point `solid_trace` queries without another solve. Memory
cost is 8 × number of nodes × number of time steps bytes, besides display frames.
"""
function simulate_pulse end

"""Dimensionless source waveform used by the transient solver (SI time)."""
function solid_pulse_drive(s::SolidLensExperiment,t::Real)
    T=s.pulse_cycles/s.frequency
    0<t<T ? sin(2π*s.frequency*t)*sinpi(t/T)^2 : 0.0
end

solid_lens_comparison(s::SolidLensExperiment) =
    [simulate(remake_experiment(s;scenario=c)) for c in (:solid,:straight,:stg)]

straight_reference(s::SolidLensExperiment)=remake_experiment(s;scenario=:straight)

function experiment_toml(s::SolidLensExperiment)
    validate_experiment(s)
    d=Dict(String(k)=>getfield(s,k) for k in fieldnames(SolidLensExperiment))
    d["schema_version"]=1; d["physics"]="solid_lens"
    d["scenario"]=String(s.scenario); d["polarization"]=String(s.polarization)
    d["bricks"]=shape_toml.(s.bricks)
    io=IOBuffer(); TOML.print(io,d;sorted=true); String(take!(io))
end

save_experiment(path::AbstractString,s::SolidLensExperiment)=write(path,experiment_toml(s))

function solid_profile_csv(r::SolidLensResult)
    io=IOBuffer(); println(io,"y_m,kinetic_energy_density_J_m3,input_power_W_m")
    for j in eachindex(r.profile_y)
        println(io,join((r.profile_y[j],r.profile_kinetic[j],r.diagnostics["input_power_W_m"]),','))
    end
    String(take!(io))
end

function lens_design_csv(d::LensDesign)
    io=IOBuffer(); println(io,"cell_index,amplitude_m,periods,phase_error_rad,calibrated_transmission_magnitude,calibrated_transmission_phase_rad,design_frequency_Hz,design_polarization")
    for (j,i) in enumerate(d.selected)
        b=d.candidates[i]
        println(io,join((j,b.amplitude,b.periods,d.phase_error[j],abs(d.transmission[i]),
            angle(d.transmission[i]),d.experiment.frequency,d.experiment.polarization),','))
    end
    String(take!(io))
end

function profile_width(y,value)
    m=argmax(value); half=value[m]/2
    value[m]>0 || return NaN
    l=m; r=m
    while l>1 && value[l-1]>=half; l-=1; end
    while r<length(y) && value[r+1]>=half; r+=1; end
    (l==1 || r==length(y)) && return NaN # beam not resolved inside section
    yl=y[l-1]+(half-value[l-1])*(y[l]-y[l-1])/(value[l]-value[l-1])
    yr=y[r]+(half-value[r])*(y[r+1]-y[r])/(value[r+1]-value[r])
    yr-yl
end

solid_examples()=[
    (id=:bulk,title="1 · Сплошное тело",question="Как выглядит волна без стенки?"),
    (id=:straight,title="2 · Прямые перемычки",question="Что дают прорези без извилистости?"),
    (id=:uniform,title="3 · Одинаковые STG",question="Достаточно ли одинаковых ячеек для фокуса?"),
    (id=:focused,title="4 · STG-линза",question="Собирает ли рассчитанная фазовая стенка волну?"),
    (id=:scrambled,title="5 · Перемешанные ячейки",question="Что меняет порядок тех же самых ячеек?"),
    (id=:missing,title="6 · Дефект линзы",question="Что происходит, если одна перемычка отсутствует?")]

function solid_preset(d::LensDesign,id::Symbol;frequency=d.experiment.frequency,polarization=d.experiment.polarization)
    s=remake_experiment(d.experiment;frequency,polarization)
    id == :bulk && return remake_experiment(s;scenario=:solid)
    id == :straight && return remake_experiment(s;scenario=:straight)
    id == :focused && return s
    if id == :uniform
        b=first(s.bricks)
        uniform=SinusoidalBrick(length=b.length,height=b.height,amplitude=5b.height/14,
            end_margin=b.end_margin,periods=2,samples=80)
        return remake_experiment(s;bricks=fill(uniform,length(s.bricks)))
    end
    id == :scrambled && return remake_experiment(s;bricks=circshift(s.bricks,3))
    id == :missing && return remake_experiment(s;missing_cells=[(length(s.bricks)+1)÷2])
    throw(ArgumentError("unknown solid experiment preset"))
end
