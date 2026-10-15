using Test, LinearAlgebra, MetaNDT, Gridap, GridapGmsh, Gmsh
BLAS.set_num_threads(1)
ext=Base.get_extension(MetaNDT,:MetaNDTGridapExt)

@testset "solid spec and CAD contract" begin
    s=SolidLensExperiment()
    @test validate_experiment(s) === s
    @test first(solid_wave_speeds(s)) > last(solid_wave_speeds(s)) > 3000
    @test_throws ArgumentError remake_experiment(s;polarization=:SH)
    @test_throws ArgumentError remake_experiment(s;missing_cells=[1,1])
    @test_throws ArgumentError remake_experiment(s;bricks=s.bricks[1:2])
    @test_throws ArgumentError remake_experiment(s;near_mesh_size=.002)
    @test_throws ArgumentError remake_experiment(s;frequency=1e6)
    @test_throws ArgumentError remake_experiment(s;source_x=0.0)
    @test_throws ArgumentError MetaNDT.validate_shape(SinusoidalBrick(end_margin=1e-10))
    @test MetaNDT.validate_shape(SinusoidalBrick(end_margin=0.0)).end_margin == 0
    @test MetaNDT.validate_shape(SinusoidalBrick(end_margin=1e-6)).end_margin == 1e-6
    @test_throws ArgumentError MetaNDT.validate_shape(RectangularBrick(height=1e-8))
    @test remake_experiment(s).bricks !== s.bricks
    @test remake_experiment(s).missing_cells !== s.missing_cells
    @test length(solid_examples()) == 6
    @test MetaNDT.profile_width(collect(-2.:2.),[0.,.5,1.,.5,0.]) ≈ 2
    @test isnan(MetaNDT.profile_width(collect(-2.:2.),ones(5)))
end

@testset "P1 sampler independent affine check" begin
    # Both orientations, diagonal/shared edges, vertices and translated coordinates.
    xy=[-2. -1. -1. -2.; 3. 3. 4. 4.]; tris=[1 1;2 4;3 3]; dofs=reshape(1:4,1,:)
    pts=[(-2.,3.),(-1.,4.),(-1.5,3.5),(-1.75,3.8),(-1.02,3.99)]
    S=ext.sensor_sampler(xy,tris,dofs,pts,4)
    @test S*ones(4) ≈ ones(length(pts))
    @test S*vec(2xy[1,:]-3xy[2,:].+7) ≈ [2p[1]-3p[2]+7 for p in pts]
    @test_throws ArgumentError ext.sensor_sampler(xy,tris,dofs,[(-2.1,3.5)],4)
end

design=design_stg_lens(SolidLensExperiment())
s=design.experiment
@testset "real calibrated geometry and replay" begin
    @test length(design.candidates) == 27
    @test maximum(abs,design.phase_error) < .4
    @test length(unique(design.selected)) > 3
    @test all(isfinite,design.transmission)
    for id in (e.id for e in solid_examples())
        spec=solid_preset(design,id)
        mktempdir() do dir
            save_experiment(joinpath(dir,"spec.toml"),spec)
            @test experiment_toml(load_experiment(joinpath(dir,"spec.toml"))) == experiment_toml(spec)
        end
    end
    @test solid_preset(design,:focused;frequency=100e3).bricks == s.bricks
    @test sort(specimen_area.(solid_preset(design,:scrambled).bricks)) == sort(specimen_area.(s.bricks))
    @test length(solid_lens_outlines(solid_preset(design,:missing))) == 12
    @test length(split(chomp(lens_design_csv(design)),'\n')) == 14
end

results=Dict(id=>simulate(solid_preset(design,id)) for id in (e.id for e in solid_examples()))
@testset "full connected wall and controls" begin
    for (id,r) in results
        @test r.diagnostics["power_balance_relative"] < 1e-8
        @test r.diagnostics["linear_residual_relative"] < 1e-8
        @test r.diagnostics["input_power_W_m"] > 0
        @test all(isfinite,r.displacement_y)
        @test all(isfinite,r.profile_kinetic)
        expected=r.experiment.scenario==:solid ? s.width*s.height :
            (s.width-first(s.bricks).length)*s.height+
            sum(specimen_area(b) for (i,b) in enumerate(r.experiment.bricks) if !(i in r.experiment.missing_cells))
        @test r.diagnostics["solid_area_m2"] ≈ expected rtol=1e-9
        println(id," ",r.diagnostics)
    end
    @test results[:focused].diagnostics["solid_area_m2"] ≈ results[:straight].diagnostics["solid_area_m2"] rtol=1e-10
    @test results[:focused].diagnostics["solid_area_m2"] ≈ results[:scrambled].diagnostics["solid_area_m2"] rtol=1e-10
    @test results[:focused].diagnostics["profile_fwhm_m"] < results[:uniform].diagnostics["profile_fwhm_m"]*.6
    @test results[:focused].diagnostics["profile_fwhm_m"] < results[:straight].diagnostics["profile_fwhm_m"]
    @test abs(results[:focused].diagnostics["profile_peak_y_m"]) < .003
    @test abs(results[:scrambled].diagnostics["profile_peak_y_m"]) > .010
    @test results[:focused].diagnostics["profile_peak_per_input_power"] < results[:straight].diagnostics["profile_peak_per_input_power"]
    @test length(split(chomp(solid_profile_csv(results[:focused])),'\n')) == 302
    zero=simulate(remake_experiment(s;source_strength=0.0))
    @test iszero(norm(zero.displacement_y))
    doubled=simulate(remake_experiment(s;source_strength=2.0))
    @test doubled.profile_kinetic ≈ 4results[:focused].profile_kinetic rtol=1e-10
    @test doubled.diagnostics["input_power_W_m"] ≈ 4results[:focused].diagnostics["input_power_W_m"] rtol=1e-10
    # Check actual FE connectivity (not only the drawn polygons) and PSD damping.
    mktempdir() do dir
        a=ext.solid_assemble(s,joinpath(dir,"mesh.msh"))
        @test norm(a.M-a.M')/norm(a.M) < 1e-12
        @test norm(a.K-a.K')/norm(a.K) < 1e-12
        @test norm(a.C-a.C')/norm(a.C) < 1e-12
        @test dot(ones(size(a.C,1)),a.C*ones(size(a.C,1))) > 0
        neighbors=[Int[] for _ in axes(a.xy,2)]
        for tri in eachcol(a.tri), i in tri, j in tri
            i!=j && push!(neighbors[i],j)
        end
        visited=falses(length(neighbors)); queue=[1]; visited[1]=true
        for node in queue, next in neighbors[node]
            if !visited[next]; visited[next]=true; push!(queue,next); end
        end
        @test all(visited)
    end
end

@testset "refinement and remote boundary" begin
    for (label,spec) in ((:fine,remake_experiment(s;mesh_size=.001,near_mesh_size=.00035)),
                        (:domain,remake_experiment(s;width=.240)))
        r=simulate(spec); control=simulate(remake_experiment(spec;scenario=:straight))
        println(label," ",r.diagnostics," straight ",control.diagnostics)
        @test r.diagnostics["profile_fwhm_m"] < control.diagnostics["profile_fwhm_m"]
        @test abs(r.diagnostics["profile_peak_y_m"]) < .003
        @test abs(r.diagnostics["profile_fwhm_m"]/results[:focused].diagnostics["profile_fwhm_m"]-1) < .1
    end
end

@testset "real pulse, energy accounting and step refinement" begin
    p=simulate_pulse(remake_experiment(s;frame_count=11))
    fine=simulate_pulse(remake_experiment(s;time_step=.05e-6,frame_count=11))
    width=MetaNDT.profile_width(p.profile_y,p.profile_integrated_kinetic)
    println("pulse ",p.diagnostics," FWHM=",width," refined ",fine.diagnostics)
    @test p.diagnostics["energy_balance_relative"] < 1e-8
    @test fine.diagnostics["energy_balance_relative"] < 1e-8
    @test iszero(norm(p.displacement_y[:,1]))
    @test maximum(abs,p.displacement_y) > 0
    @test .014 < width < .022
    @test norm(p.profile_integrated_kinetic-fine.profile_integrated_kinetic)/norm(fine.profile_integrated_kinetic) < .03
    @test p.diagnostics["input_work_J_m"] > 0
end
