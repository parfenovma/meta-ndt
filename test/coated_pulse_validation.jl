using Test,MetaNDT,LinearAlgebra
BLAS.set_num_threads(1)
function fourier_integral(times,values,f)
    y=values.*exp.(2π*im*f.*times)
    sum(diff(times).*(y[1:end-1]+y[2:end])/2)
end
function resample_trace(trace,times)
    [begin
        j=clamp(searchsortedlast(trace.times,t),1,length(trace.times)-1)
        w=(t-trace.times[j])/(trace.times[j+1]-trace.times[j])
        (1-w)*trace.scattered[j]+w*trace.scattered[j+1]
    end for t in times]
end

@testset "Transient versus independently validated harmonic response" begin
    point=(.015,.006);s=CoatedPulseExperiment();r=simulate_pulse(s);trace=coated_trace(r,point)
    # IMPORTANT: the default 460 μs movie cuts a real ringdown tail. It is
    # valid as a finite-window experiment, not as an infinite-time Fourier
    # transfer estimate. Use 4 ms for this separate frequency-domain check.
    # At ka<0.22 modes 0–3 suffice; the reference retains the full 0–12 modes.
    longspec=CoatedPulseExperiment(model=remake_experiment(s.model;modes=3),duration=.004,frame_count=2)
    longtrace=coated_trace(simulate_pulse(longspec),point)
    input=[coated_incident_pressure(longspec,(0.,0.),t) for t in longtrace.times]
    H=fourier_integral(longtrace.times,longtrace.scattered,s.model.frequency)/fourier_integral(longtrace.times,input,s.model.frequency)
    expected=coated_pressure(simulate(s.model),point;scattered=true)
    println("10 kHz transfer: transient ",H,"; harmonic ",expected,"; relative ",abs(H-expected)/abs(expected))
    @test abs(H-expected)/abs(expected)<.02
    fine=simulate_pulse(remake_experiment(s;time_step=s.time_step/2));tf=coated_trace(fine,point)
    temporal=norm(trace.scattered-resample_trace(tf,trace.times))/norm(resample_trace(tf,trace.times))
    println("dt/2 trace change = ",temporal)
    @test temporal<.03
    m=remake_experiment(s.model;steel_elements=24,coating_elements=96,modes=16)
    refined=simulate_pulse(remake_experiment(s;model=m,water_mesh_size=s.water_mesh_size/2,water_interface_divisions=32,air_elements=96))
    tr=coated_trace(refined,point)
    spatial=norm(trace.scattered-tr.scattered)/norm(tr.scattered)
    println("Radial + modal refinement trace change = ",spatial)
    @test spatial<.05
    @test refined.diagnostics["energy_balance_relative"]<1e-7
    # The boundary return estimate is not an exact discrete causality proof.
    # Actually move it out, keeping the near mesh, and compare the microphones.
    larger=simulate_pulse(remake_experiment(s;outer_radius=1.5s.outer_radius))
    tl=coated_trace(larger,point)
    @test norm(trace.scattered-tl.scattered)/norm(trace.scattered)<1e-5
end
