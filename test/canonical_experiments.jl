using Test, MetaNDT, LinearAlgebra

@testset "Canonical membrane and causal ports" begin
    m=membrane_model()
    @test isposdef(Symmetric(m.M)) && isposdef(Symmetric(m.K))
    @test real(membrane_response(m,150.).mass)<0
    @test real(membrane_response(m,800.).mass)>0
    for f in (80.,150.,250.,800.)
        r=membrane_response(m,f)
        @test abs2(r.r)+abs2(r.t)≈1 atol=1e-12
    end
    r=membrane_record(m,192.)
    @test r.diagnostics["energy_balance_relative"]<1e-9
    @test all(isfinite,r.field)
    @test iszero(r.internal[:,1])
    @test maximum(abs,r.field[:,end])<maximum(abs,r.field)
    a=membrane_response(m,250.;back_depth=.025,loss_rate=0.)
    @test a.absorption≈0 atol=1e-10
    @test a.t==0
    wall=membrane_record(m,250.;back_depth=.025,wall=true)
    @test wall.diagnostics["source_work_J"]==0
    @test iszero(wall.internal)
    @test_throws ArgumentError membrane_model(added_mass=-1.)
end

@testset "Helmholtz network, not a painted stop band" begin
    @test real(helmholtz_modulus(29000.))<0
    @test real(helmholtz_modulus(18000.))>0
    @test real(helmholtz_modulus(29000.;open=false))>0
    closed=helmholtz_record(29000.;open=false)
    opened=helmholtz_record(29000.)
    @test opened.diagnostics["energy_balance_relative"]<1e-9
    @test closed.diagnostics["energy_balance_relative"]<1e-9
    @test maximum(abs,opened.internal)>0
    @test sum(abs2,opened.field[end,:])<.5sum(abs2,closed.field[end,:])
end

@testset "Spatial pulses are computed time histories" begin
    flat=canonical_spatial_record(:cloak,:flat;nx=65,ny=41)
    bump=canonical_spatial_record(:cloak,:bump;nx=65,ny=41)
    @test all(isfinite,bump.field)
    @test count(bump.valid)<count(flat.valid)
    @test bump.diagnostics["max_wave_energy_ratio"]<1.08
    @test_throws ArgumentError canonical_spatial_record(:unknown,:flat)
end
