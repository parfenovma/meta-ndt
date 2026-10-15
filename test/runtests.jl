using LinearAlgebra
using MetaNDT
using Test

# This backend uses only library dependencies, so the standalone library CI
# must exercise its physics without requiring the optional Gridap/Gmsh stack.
include("coated_fem.jl")
include("coated_pulse.jl")

@testset "noise model" begin
    model = GaussianReceiverNoise(2.5)
    first = sample_noise(model, 31, 3)
    second = sample_noise(model, 31, 3)
    different = sample_noise(model, 31, 4)
    @test first == second
    @test first != different
    @test length(first) == 31
    @test_throws ArgumentError GaussianReceiverNoise(-1)
end

@testset "paired experiment" begin
    primary = [2.0, 1.0, 0.5, 1.5]
    control = [1.0, 0.5, 0.25, 0.75]
    model = GaussianReceiverNoise(0.8)
    measurement = simulate_paired_measurement(
        primary,
        control;
        noise_model=model,
        signal_present=true,
        measurement_index=2,
    )
    @test measurement.primary .- primary ≈ measurement.noise
    @test measurement.control .- control ≈ measurement.noise
    @test measurement.noise_sigma == 0.8
    @test measurement.signal_present

    decisions = detect_pair(primary, control, measurement)
    @test decisions.primary.probability_kind == :Pd
    @test 0 <= decisions.primary.probability <= 1
    @test 0 <= decisions.control.probability <= 1
end

@testset "absence and zero noise" begin
    template = [1.0, 2.0, 1.0]
    absent = simulate_paired_measurement(
        template,
        template;
        noise_model=GaussianReceiverNoise(1.0),
        signal_present=false,
        measurement_index=3,
    )
    result = detect_pair(template, template, absent).primary
    @test result.probability_kind == :Pfa
    @test result.probability ≈ 0.05 atol=2e-7

    noiseless = matched_filter(template, template, 0.0)
    @test noiseless.found
    @test isinf(noiseless.score_sigma)
    @test noiseless.probability == 1
end

@testset "complex phase-sensitive traces" begin
    template = ComplexF64[1 + im, 2 - im, -1 + 0.5im]
    measurement = 0.75 .* template
    detector = MatchedFilterDetector(template, 1.0)
    result = detect(detector, measurement)
    @test result.score_sigma ≈ 0.75 * norm(template)
    @test rms_amplitude(template) ≈ norm(template) / sqrt(3)
end

include("canonical_experiments.jl")
