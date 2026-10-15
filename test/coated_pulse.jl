using Test,MetaNDT,LinearAlgebra
BLAS.set_num_threads(1)

@testset "Direct water–steel–coating–air transient" begin
    s=CoatedPulseExperiment()
    @test validate_experiment(s)===s
    mktemp() do path,io
        close(io);save_experiment(path,s)
        @test experiment_toml(load_experiment(path))==experiment_toml(s)
    end
    cases=coated_pulse_comparison(s)
    point=(.015,.006)
    traces=coated_trace.(cases,Ref(point))
    for (i,r) in enumerate(cases)
        t=traces[i]
        @test t.pressure[1]==t.scattered[1]==0
        @test t.pressure≈t.scattered+[coated_incident_pressure(r.experiment,point,ti) for ti in t.times]
        @test r.diagnostics["energy_balance_relative"]<1e-7
        println("$(r.experiment.model.scenario): ",r.diagnostics)
        if i>1
            @test r.diagnostics["interface_skew_relative"]<1e-14
            @test r.diagnostics["peak_displacement_m"]>0
            @test r.diagnostics["peak_air_pressure_Pa"]>0
            @test maximum(r.solid_energy)>0 && maximum(r.air_energy)>0
            @test maximum(abs,t.scattered)>0
            @test r.diagnostics["outer_return_lower_bound_s"]>s.duration
        end
    end
    @test all(iszero,traces[1].scattered)
    @test maximum(abs,traces[2].scattered-traces[3].scattered)>0
    @test_throws ArgumentError coated_trace(cases[3],(0.,0.))
    @test_throws ArgumentError coated_trace(cases[3],(.041,0.))
    @test coated_trace(cases[1],(0.,0.)).scattered==zeros(length(cases[1].times))
    @test all(t->coated_pulse_drive(s,t)==0,range(s.pulse_cycles/s.model.frequency,s.duration;length=50))
    for i in (1,60,141)
        @test coated_trace(cases[3],point;indices=[cases[3].frame_steps[i]]).pressure[1]≈traces[3].pressure[cases[3].frame_steps[i]]
    end
    near=coated_trace(cases[1],(-.020,0.));far=coated_trace(cases[1],(.020,0.))
    @test near.times[argmax(abs.(near.pressure))]<far.times[argmax(abs.(far.pressure))]
    println("Scattered RMS bare/coated: ",norm(traces[2].scattered)," / ",norm(traces[3].scattered))
end
