module MetaNDT

using LinearAlgebra

export AbstractDetector,
       AbstractNoiseModel,
       DetectionResult,
       GaussianReceiverNoise,
       MatchedFilterDetector,
       PairedMeasurement,
       DEFAULT_THRESHOLD_SIGMA,
       detect,
       detect_pair,
       detection_probability,
       matched_filter,
       normalized_noise,
       probability_curve,
       rms_amplitude,
       sample_noise,
       simulate_paired_measurement

abstract type AbstractDetector end
abstract type AbstractNoiseModel end

const DEFAULT_THRESHOLD_SIGMA = 1.6448536269514722

"""Deterministic white Gaussian receiver noise for controlled comparisons."""
struct GaussianReceiverNoise <: AbstractNoiseModel
    sigma::Float64
    base_seed::UInt32
    seed_step::UInt32

    function GaussianReceiverNoise(
        sigma::Real,
        base_seed::UInt32=UInt32(0x00c0ffee),
        seed_step::UInt32=UInt32(7919),
    )
        sigma >= 0 || throw(ArgumentError("noise sigma must be non-negative"))
        new(Float64(sigma), base_seed, seed_step)
    end
end

"""A known-template detector after an explicit scalar noise normalization."""
struct MatchedFilterDetector{T<:Number,V<:AbstractVector{T}} <: AbstractDetector
    template::V
    noise_sigma::Float64
    threshold_sigma::Float64

    function MatchedFilterDetector(
        template::V,
        noise_sigma::Real;
        threshold_sigma::Real=DEFAULT_THRESHOLD_SIGMA,
    ) where {T<:Number,V<:AbstractVector{T}}
        isempty(template) && throw(ArgumentError("template must not be empty"))
        norm(template) > 0 || throw(ArgumentError("template norm must be positive"))
        noise_sigma >= 0 || throw(ArgumentError("noise sigma must be non-negative"))
        new{T,V}(template, Float64(noise_sigma), Float64(threshold_sigma))
    end
end

"""One paired observation with the same additive noise in both paths."""
struct PairedMeasurement{P<:Number,C<:Number,N<:Number}
    primary::Vector{P}
    control::Vector{C}
    noise::Vector{N}
    noise_sigma::Float64
    signal_present::Bool
    measurement_index::Int
end

"""Fixed-threshold detector output and its model probability."""
Base.@kwdef struct DetectionResult
    score_sigma::Float64
    found::Bool
    probability::Float64
    probability_kind::Symbol
end

rms_amplitude(values::AbstractVector{<:Number}) = norm(values) / sqrt(length(values))

function _uniform_sequence(count::Integer, seed::UInt32)
    count >= 0 || throw(ArgumentError("count must be non-negative"))
    values = Vector{Float64}(undef, count)
    state = seed
    for index in eachindex(values)
        state += UInt32(0x6d2b79f5)
        value = state
        value = (value ⊻ (value >> 15)) * (value | UInt32(1))
        value = value ⊻ (value + (value ⊻ (value >> 7)) * (value | UInt32(61)))
        value = value ⊻ (value >> 14)
        values[index] = Float64(value) / 4294967296.0
    end
    return values
end

"""A deterministic N(0, 1) sequence independent of Julia's global RNG."""
function normalized_noise(
    count::Integer,
    measurement_index::Integer;
    base_seed::UInt32=UInt32(0x00c0ffee),
    seed_step::UInt32=UInt32(7919),
)
    count >= 0 || throw(ArgumentError("count must be non-negative"))
    measurement_index >= 0 || throw(ArgumentError("measurement_index must be non-negative"))
    pair_count = cld(count, 2)
    seed = base_seed + UInt32(measurement_index) * seed_step
    uniforms = _uniform_sequence(2 * pair_count, seed)
    noise = Float64[]
    sizehint!(noise, count)
    for pair in 1:pair_count
        first_uniform = max(uniforms[2pair - 1], eps(Float64))
        second_uniform = uniforms[2pair]
        radius = sqrt(-2 * log(first_uniform))
        push!(noise, radius * cos(2π * second_uniform))
        length(noise) < count && push!(noise, radius * sin(2π * second_uniform))
    end
    return noise
end

function sample_noise(model::GaussianReceiverNoise, count::Integer, measurement_index::Integer)
    model.sigma .* normalized_noise(
        count,
        measurement_index;
        base_seed=model.base_seed,
        seed_step=model.seed_step,
    )
end

"""
    simulate_paired_measurement(primary, control; noise_model, signal_present=true, measurement_index=0)

Build a controlled A/B measurement with one shared additive-noise realization.
The names `primary` and `control` deliberately carry no geometry assumption.
"""
function simulate_paired_measurement(
    primary_template::AbstractVector{<:Number},
    control_template::AbstractVector{<:Number};
    noise_model::AbstractNoiseModel,
    signal_present::Bool=true,
    measurement_index::Integer=0,
)
    length(primary_template) == length(control_template) ||
        throw(ArgumentError("paired templates must have equal length"))
    noise = sample_noise(noise_model, length(primary_template), measurement_index)
    signal_scale = signal_present ? 1.0 : 0.0
    primary = signal_scale .* primary_template .+ noise
    control = signal_scale .* control_template .+ noise
    return PairedMeasurement(
        collect(primary),
        collect(control),
        collect(noise),
        noise_model.sigma,
        signal_present,
        Int(measurement_index),
    )
end

function _erf_approx(value::Real)
    sign_value = value < 0 ? -1.0 : 1.0
    x = abs(Float64(value))
    a1, a2, a3 = 0.254829592, -0.284496736, 1.421413741
    a4, a5, p = -1.453152027, 1.061405429, 0.3275911
    t = 1 / (1 + p * x)
    polynomial = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
    return sign_value * (1 - polynomial * exp(-x^2))
end

_normal_cdf(value::Real) = 0.5 * (1 + _erf_approx(value / sqrt(2)))

function detection_probability(
    template::AbstractVector{<:Number},
    noise_sigma::Real;
    threshold_sigma::Real=DEFAULT_THRESHOLD_SIGMA,
)
    noise_sigma >= 0 || throw(ArgumentError("noise sigma must be non-negative"))
    noise_sigma == 0 && return 1.0
    matched_snr = norm(template) / noise_sigma
    return _normal_cdf(matched_snr - threshold_sigma)
end

"""
    detect(detector, measurement; hypothesis_present=true)

Project a real or complex, already conditioned trace onto a known template.
For complex traces the statistic is `real(dot(template, measurement))`. The
caller owns calibration and whitening; MetaNDT never guesses them.
"""
function detect(
    detector::MatchedFilterDetector,
    measurement::AbstractVector{<:Number};
    hypothesis_present::Bool=true,
)
    length(measurement) == length(detector.template) ||
        throw(ArgumentError("measurement/template length mismatch"))
    score = real(dot(detector.template, measurement)) / norm(detector.template)
    score_sigma = detector.noise_sigma == 0 ? (score > 0 ? Inf : 0.0) : score / detector.noise_sigma
    probability = if hypothesis_present
        detection_probability(
            detector.template,
            detector.noise_sigma;
            threshold_sigma=detector.threshold_sigma,
        )
    else
        detector.noise_sigma == 0 ? 0.0 : 1 - _normal_cdf(detector.threshold_sigma)
    end
    return DetectionResult(
        score_sigma=score_sigma,
        found=score_sigma >= detector.threshold_sigma,
        probability=probability,
        probability_kind=hypothesis_present ? :Pd : :Pfa,
    )
end

function matched_filter(
    measurement::AbstractVector{<:Number},
    template::AbstractVector{<:Number},
    noise_sigma::Real;
    hypothesis_present::Bool=true,
    threshold_sigma::Real=DEFAULT_THRESHOLD_SIGMA,
)
    detector = MatchedFilterDetector(template, noise_sigma; threshold_sigma=threshold_sigma)
    return detect(detector, measurement; hypothesis_present=hypothesis_present)
end

function detect_pair(
    primary_template::AbstractVector{<:Number},
    control_template::AbstractVector{<:Number},
    measurement::PairedMeasurement;
    threshold_sigma::Real=DEFAULT_THRESHOLD_SIGMA,
)
    return (
        primary=matched_filter(
            measurement.primary,
            primary_template,
            measurement.noise_sigma;
            hypothesis_present=measurement.signal_present,
            threshold_sigma=threshold_sigma,
        ),
        control=matched_filter(
            measurement.control,
            control_template,
            measurement.noise_sigma;
            hypothesis_present=measurement.signal_present,
            threshold_sigma=threshold_sigma,
        ),
    )
end

function probability_curve(primary_template, control_template, noise_sigmas)
    sigma = Float64.(collect(noise_sigmas))
    return (
        sigma=sigma,
        primary=[detection_probability(primary_template, value) for value in sigma],
        control=[detection_probability(control_template, value) for value in sigma],
    )
end

end
