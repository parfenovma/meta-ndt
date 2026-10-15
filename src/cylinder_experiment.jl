using TOML

export CylinderExperiment, CylinderResult, validate_experiment, remake_experiment,
       simulate, simulate_comparison, experiment_toml, save_experiment,
       load_experiment, sensor_csv, sweep_experiments, cylinder_lattice,
       FrequencySweepResult, frequency_sweep, frequency_sweep_csv

"""Linear 2D water/solid experiment. All inputs use SI units; solid is plane strain."""
Base.@kwdef struct CylinderExperiment
    scenario::Symbol = :elastic
    radius::Float64 = 0.010
    centers::Vector{NTuple{2,Float64}} = [(0.0, 0.0)]
    shape::Union{Nothing,RectangularBrick,SinusoidalBrick} = nothing
    width::Float64 = 0.160
    height::Float64 = 0.100
    water_density::Float64 = 1000.0
    water_speed::Float64 = 1480.0
    solid_density::Float64 = 2700.0
    young_modulus::Float64 = 70e9
    poisson_ratio::Float64 = 0.33
    source::NTuple{2,Float64} = (-0.045, 0.0)
    source_width::Float64 = 0.003
    source_strength::Float64 = 1e-5
    frequency::Float64 = 60e3
    pulse_cycles::Float64 = 2.0
    sensors::Vector{NTuple{2,Float64}} = [(0.040, 0.0), (-0.030, 0.020), (0.025, -0.025)]
    mesh_size::Float64 = 0.002
    time_step::Float64 = 0.25e-6
    duration::Float64 = 150e-6
    frame_count::Int = 81
    absorbing_boundary::Bool = true
end

function validate_experiment(s::CylinderExperiment)
    s.scenario in (:empty, :rigid, :elastic) || throw(ArgumentError("unknown scenario"))
    positive = (:radius, :width, :height, :water_density, :water_speed, :solid_density,
                :young_modulus, :source_width, :frequency, :pulse_cycles,
                :mesh_size, :time_step, :duration)
    all(k -> isfinite(getfield(s,k)) && getfield(s,k) > 0, positive) ||
        throw(ArgumentError("material, length and time parameters must be finite and positive"))
    isfinite(s.source_strength) || throw(ArgumentError("source strength must be finite"))
    -1 < s.poisson_ratio < 0.5 || throw(ArgumentError("Poisson ratio must be in (-1, 0.5)"))
    !isnothing(s.shape) && validate_shape(s.shape)
    ex,ey=specimen_half_extents(s)
    max(ex,ey) < min(s.width,s.height)/4 || throw(ArgumentError("inclusion too large for domain"))
    s.frame_count >= 2 || throw(ArgumentError("at least two frames required"))
    s.mesh_size <= specimen_min_thickness(s)/4 || throw(ArgumentError("mesh needs at least four elements across inclusion thickness"))
    s.time_step*s.frequency <= 0.04 || throw(ArgumentError("need at least 25 time steps per carrier period"))
    s.duration > s.pulse_cycles/s.frequency || throw(ArgumentError("duration must exceed pulse duration"))
    isempty(s.sensors) && throw(ArgumentError("at least one sensor required"))
    isempty(s.centers) && throw(ArgumentError("at least one cylinder center required"))
    for p in s.centers
        all(isfinite,p) || throw(ArgumentError("nonfinite cylinder center"))
        abs(p[1])+ex+2s.mesh_size < s.width/2 &&
            abs(p[2])+ey+2s.mesh_size < s.height/2 ||
            throw(ArgumentError("inclusion must fit inside the domain with a mesh clearance"))
    end
    for i in eachindex(s.centers), j in 1:i-1
        a,b=s.centers[i],s.centers[j]
        clear=isnothing(s.shape) ? hypot(a[1]-b[1],a[2]-b[2]) >= 2s.radius+2s.mesh_size :
            abs(a[1]-b[1]) >= 2ex+2s.mesh_size || abs(a[2]-b[2]) >= 2ey+2s.mesh_size
        clear || throw(ArgumentError("inclusions overlap or their bounding-box water gap is under two mesh sizes"))
    end
    for (i,p) in enumerate([s.source, s.sensors...])
        all(isfinite,p) || throw(ArgumentError("nonfinite source/sensor coordinates"))
        abs(p[1]) < s.width/2 && abs(p[2]) < s.height/2 ||
            throw(ArgumentError("source and sensors must be inside the water domain"))
        if s.scenario != :empty
            all(c -> distance_to_inclusion(s,p,c) >
                (i == 1 ? 3s.source_width : s.mesh_size/2), s.centers) ||
                throw(ArgumentError("source/sensor overlaps an inclusion"))
        end
    end
    abs(s.source[1])+3s.source_width < s.width/2 &&
        abs(s.source[2])+3s.source_width < s.height/2 ||
        throw(ArgumentError("source must be at least three source widths from outer boundary"))
    return s
end

function remake_experiment(s::CylinderExperiment; kwargs...)
    values = NamedTuple{fieldnames(CylinderExperiment)}(Tuple(getfield(s,k) for k in fieldnames(CylinderExperiment)))
    parameters = merge(values, (;kwargs...))
    # Each experiment owns its geometry and sensor lists.
    parameters = merge(parameters, (sensors=copy(parameters.sensors),centers=copy(parameters.centers)))
    validate_experiment(CylinderExperiment(; parameters...))
end

"""Plain arrays: FE mesh, actual transient snapshots, pressure traces and energy diagnostics."""
struct CylinderResult
    experiment::CylinderExperiment
    times::Vector{Float64}
    sensor_pressure::Matrix{Float64}
    frame_times::Vector{Float64}
    fluid_nodes::Matrix{Float64}
    fluid_triangles::Matrix{Int}
    pressure::Matrix{Float64}
    solid_nodes::Matrix{Float64}
    solid_triangles::Matrix{Int}
    displacement_x::Matrix{Float64}
    displacement_y::Matrix{Float64}
    energy::Vector{Float64}
    solid_energy::Vector{Float64}
    input_work::Vector{Float64}
    absorbed_energy::Vector{Float64}
    diagnostics::Dict{String,Float64}
end

"""Load Gridap, GridapGmsh and Gmsh to enable the optional real FEM implementation."""
function simulate end

function simulate_comparison(s::CylinderExperiment; kwargs...)
    [simulate(remake_experiment(s; scenario=c); kwargs...) for c in (:empty,:rigid,:elastic)]
end

function experiment_toml(s::CylinderExperiment)
    validate_experiment(s)
    d = Dict(String(k) => getfield(s,k) for k in fieldnames(CylinderExperiment))
    d["schema_version"] = 3
    d["scenario"] = String(s.scenario)
    d["source"] = collect(s.source)
    d["sensors"] = collect.(s.sensors)
    d["centers"] = collect.(s.centers)
    d["shape"] = shape_toml(s.shape)
    io = IOBuffer()
    TOML.print(io,d; sorted=true)
    String(take!(io))
end

save_experiment(path::AbstractString,s::CylinderExperiment) = write(path,experiment_toml(s))

function load_experiment(path::AbstractString)
    d = TOML.parsefile(path)
    if get(d,"physics","")=="coated_pulse"
        pop!(d,"physics");pop!(d,"schema_version",0)==1 || throw(ArgumentError("unsupported coated pulse schema"))
        model=pop!(d,"model")
        pop!(model,"physics","")=="coated_cylinder" && pop!(model,"schema_version",0)==1 ||
            throw(ArgumentError("unsupported nested coated-cylinder schema"))
        model["scenario"]=Symbol(model["scenario"])
        m=CoatedCylinderExperiment(;(Symbol(k)=>v for (k,v) in model)...)
        return validate_experiment(CoatedPulseExperiment(;model=m,(Symbol(k)=>v for (k,v) in d)...))
    end
    if get(d,"physics","") == "coated_cylinder"
        pop!(d,"physics"); pop!(d,"schema_version",0) == 1 || throw(ArgumentError("unsupported coated cylinder schema"))
        d["scenario"]=Symbol(d["scenario"])
        return validate_experiment(CoatedCylinderExperiment(; (Symbol(k)=>v for (k,v) in d)...))
    end
    if get(d,"physics","") == "solid_lens"
        pop!(d,"physics"); pop!(d,"schema_version",0) == 1 || throw(ArgumentError("unsupported solid lens schema"))
        d["scenario"]=Symbol(d["scenario"]); d["polarization"]=Symbol(d["polarization"])
        d["bricks"]=SinusoidalBrick[shape_from_toml(b) for b in d["bricks"]]
        return validate_experiment(SolidLensExperiment(; (Symbol(k)=>v for (k,v) in d)...))
    end
    version=pop!(d,"schema_version",0)
    version in (1,2,3) || throw(ArgumentError("unsupported experiment schema"))
    if version == 1
        d["centers"] = [[0.0,0.0]]
    end
    d["shape"]=version < 3 ? nothing : shape_from_toml(d["shape"])
    d["scenario"] = Symbol(d["scenario"])
    d["source"] = Tuple(Float64.(d["source"]))
    d["sensors"] = [Tuple(Float64.(p)) for p in d["sensors"]]
    d["centers"] = [Tuple(Float64.(p)) for p in d["centers"]]
    validate_experiment(CylinderExperiment(; (Symbol(k)=>v for (k,v) in d)...))
end

shape_toml(::Nothing) = Dict("kind"=>"circle")
shape_toml(b::RectangularBrick) = Dict("kind"=>"rectangle","length"=>b.length,"height"=>b.height)
shape_toml(b::SinusoidalBrick) = merge(Dict(String(k)=>getfield(b,k) for k in fieldnames(SinusoidalBrick)),Dict("kind"=>"stg"))

function shape_from_toml(value)
    d=copy(value); kind=pop!(d,"kind","")
    kind == "circle" && (isempty(d) || throw(ArgumentError("unexpected circle shape fields")); return nothing)
    parameters=(; (Symbol(k)=>v for (k,v) in d)...)
    kind == "rectangle" && return validate_shape(RectangularBrick(;parameters...))
    kind == "stg" && return validate_shape(SinusoidalBrick(;parameters...))
    throw(ArgumentError("unsupported specimen shape"))
end

"""Centers of a rectangular lattice; columns run along x, rows along y, pitch in m."""
function cylinder_lattice(columns::Integer,rows::Integer,pitch::Real; center=(0.0,0.0))
    columns > 0 && rows > 0 || throw(ArgumentError("positive lattice dimensions required"))
    isfinite(pitch) && pitch > 0 || throw(ArgumentError("positive finite pitch required"))
    length(center) == 2 && all(isfinite,center) || throw(ArgumentError("invalid lattice center"))
    [(Float64(center[1]+(i-(columns+1)/2)*pitch),
      Float64(center[2]+(j-(rows+1)/2)*pitch)) for j in 1:rows for i in 1:columns]
end

"""Paired pulse sweeps. Ratios are sensor pressure L2 norms, NOT energy transmission."""
struct FrequencySweepResult
    frequencies::Vector{Float64}
    comparisons::Vector{Vector{CylinderResult}}
    ratios::Matrix{Float64}
    reference_integrals::Vector{Float64}
    sensor::Int
    time_window::NTuple{2,Float64}
end

function pressure_integral(times,pressure,time_window)
    lo,hi=time_window
    integral=0.0
    for i in 1:length(times)-1
        a,b=max(lo,times[i]),min(hi,times[i+1])
        b <= a && continue
        dt=times[i+1]-times[i]
        pa=pressure[i]+(pressure[i+1]-pressure[i])*(a-times[i])/dt
        pb=pressure[i]+(pressure[i+1]-pressure[i])*(b-times[i])/dt
        integral+=(b-a)*(pa^2+pa*pb+pb^2)/3
    end
    integral
end

"""
Run empty/rigid/elastic cases at each frequency with a common mesh and time step.
The common mesh is refined to >=14 nominal water elements per shortest wavelength;
the time step is shortened if needed for >=25 samples per highest carrier period.
Only two field frames per run are retained; traces and peak diagnostics are complete.
Return sqrt(integral(p_case^2 dt)/integral(p_empty^2 dt)) at the chosen sensor/window.
This local, finite-pulse observable can exceed one and is not a power coefficient.
"""
function frequency_sweep(s::CylinderExperiment,frequencies; sensor::Integer=1,
                         time_window=(0.0,s.duration),include_straight_reference::Bool=false)
    validate_experiment(s)
    f=Float64.(collect(frequencies))
    !isempty(f) && all(x -> isfinite(x) && x > 0,f) && all(diff(f) .> 0) ||
        throw(ArgumentError("frequencies must be positive, finite and strictly increasing"))
    sensor in eachindex(s.sensors) || throw(ArgumentError("invalid sweep sensor"))
    length(time_window) == 2 && all(isfinite,time_window) &&
        0 <= time_window[1] < time_window[2] <= s.duration ||
        throw(ArgumentError("time window must lie within the experiment"))
    s.source_strength != 0 || throw(ArgumentError("a sweep needs a nonzero reference signal"))
    !include_straight_reference || s.shape isa SinusoidalBrick ||
        throw(ArgumentError("straight sweep control requires an STG specimen"))
    window=(Float64(time_window[1]),Float64(time_window[2]))
    common=remake_experiment(s; mesh_size=min(s.mesh_size,s.water_speed/(14maximum(f))),
        time_step=min(s.time_step,0.04/maximum(f)),frame_count=2)
    # Validate every point before spending time on any solver call.
    specs=[remake_experiment(common;frequency=frequency) for frequency in f]
    comparisons=Vector{CylinderResult}[]
    ratios=zeros(include_straight_reference ? 4 : 3,length(f)); references=zeros(length(f))
    for (j,spec) in enumerate(specs)
        cases=simulate_comparison(spec)
        include_straight_reference && push!(cases,simulate(straight_reference(cases[3].experiment)))
        empty=cases[1]
        ref=pressure_integral(empty.times,view(empty.sensor_pressure,sensor,:),window)
        full=pressure_integral(empty.times,view(empty.sensor_pressure,sensor,:),(0.0,s.duration))
        isfinite(ref) && ref > 0 && ref >= 1e-6*full ||
            throw(ArgumentError("reference signal in sweep window is too small; widen or move the window"))
        references[j]=ref
        for (i,r) in enumerate(cases)
            value=pressure_integral(r.times,view(r.sensor_pressure,sensor,:),window)
            isfinite(value) || throw(ArgumentError("nonfinite pressure integral"))
            ratios[i,j]=sqrt(value/ref)
        end
        push!(comparisons,cases)
    end
    FrequencySweepResult(f,comparisons,ratios,references,Int(sensor),window)
end

"""CSV of paired sensor norm ratios, with the measurement window and reference norm."""
function frequency_sweep_csv(r::FrequencySweepResult)
    io=IOBuffer()
    control_header=size(r.ratios,1) == 4 ? ",straight_ratio" : ""
    println(io,"frequency_Hz,sensor_index,window_start_s,window_end_s,reference_p2_integral_Pa2_s,empty_ratio,rigid_ratio,elastic_ratio$(control_header),mesh_size_m,time_step_s,max_energy_balance_relative")
    for j in eachindex(r.frequencies)
        cases=r.comparisons[j]; spec=cases[end].experiment
        println(io,join((r.frequencies[j],r.sensor,r.time_window...,r.reference_integrals[j],
            r.ratios[:,j]...,spec.mesh_size,cases[end].diagnostics["time_step_s"],
            maximum(c.diagnostics["energy_balance_relative"] for c in cases)),","))
    end
    String(take!(io))
end

function sensor_csv(r::CylinderResult)
    io = IOBuffer()
    println(io, "time_s,",join(["sensor_$(i)_pressure_Pa" for i in axes(r.sensor_pressure,1)],","))
    for j in eachindex(r.times)
        println(io, r.times[j],",",join(r.sensor_pressure[:,j],","))
    end
    String(take!(io))
end

"""A small sequential sweep using exactly the same solver as Pluto."""
function sweep_experiments(experiments; kwargs...)
    map(experiments) do s
        result=simulate(s; kwargs...)
        (experiment=result.experiment, result=result)
    end
end
