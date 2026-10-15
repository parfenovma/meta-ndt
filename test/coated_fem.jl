using MetaNDT, SpecialFunctions, LinearAlgebra, Test
BLAS.set_num_threads(1)

# Independent strong-form solution: longitudinal/transverse Bessel potentials
# in each elastic annulus, exact acoustic fluids. No FEM/mixed assembly calls.
function analytic_coated(s,n)
    ω=2π*s.frequency; k=ω/s.water_speed
    dj(f,n,z)=(f(n-1,z)-f(n+1,z))/2
    ddj(f,n,z)=-dj(f,n,z)/z-(1-n^2/z^2)*f(n,z)
    μs=s.steel_young/(2(1+s.steel_poisson))
    λs=s.steel_young*s.steel_poisson/((1+s.steel_poisson)*(1-2s.steel_poisson))
    μc=s.coating_density*s.coating_shear_speed^2*(1-im*s.coating_shear_loss)
    Kc=s.coating_density*(s.coating_pressure_speed^2-4s.coating_shear_speed^2/3)
    materials=[(s.steel_density,complex(λs),complex(μs))]
    s.scenario == :coated && push!(materials,(s.coating_density,Kc-2μc/3,μc))
    function basis(r,mat)
        ρ,λ,μ=mat; kp=ω/sqrt((λ+2μ)/ρ); ks=ω/sqrt(μ/ρ)
        B=zeros(ComplexF64,4,n == 0 ? 2 : 4)
        for (j,Z) in enumerate((besselj,bessely))
            p=Z(n,kp*r); dp=kp*dj(Z,n,kp*r); ddp=kp^2*ddj(Z,n,kp*r)
            B[:,j]=[dp,-n*p/r,-λ*kp^2*p+2μ*ddp,-2μ*n*(dp/r-p/r^2)]
            if n > 0
                t=Z(n,ks*r); dt=ks*dj(Z,n,ks*r); ddt=ks^2*ddj(Z,n,ks*r)
                B[:,j+2]=[n*t/r,-dt,2μ*n*(dt/r-t/r^2),μ*(-ddt+dt/r-n^2*t/r^2)]
            end
        end
        n == 0 ? B[[1,3],:] : B
    end
    dim=n == 0 ? 2 : 4; traction=n == 0 ? 2 : 3
    count=length(materials); R=s.outer_radius+(count == 2 ? s.coating_thickness : 0)
    ka=ω/s.air_speed; b=s.inner_radius
    Jair=besselj(n,ka*b); dJair=ka*dj(besselj,n,ka*b)
    H=hankelh1(n,k*R); dH=k*dj(hankelh1,n,k*R)
    Zw=s.water_density*ω^2*H/dH
    incident=(n == 0 ? 1 : 2)*im^n*s.incident_pressure
    dp0=incident*k*dj(besselj,n,k*R); g=incident*besselj(n,k*R)-H*dp0/dH
    A=zeros(ComplexF64,dim*count+1,dim*count+1); rhs=zeros(ComplexF64,dim*count+1)
    inner=basis(b,materials[1]); outer=basis(R,materials[end]); rows=1
    A[rows,1:dim]=inner[traction,:]; A[rows,end]=Jair; rows+=1
    A[rows,1:dim]=-s.air_density*ω^2*inner[1,:]; A[rows,end]=dJair; rows+=1
    if n > 0
        A[rows,1:dim]=inner[4,:]; rows+=1
    end
    if count == 2
        A[rows:rows+dim-1,1:dim]=basis(s.outer_radius,materials[1])
        A[rows:rows+dim-1,dim+1:2dim]=-basis(s.outer_radius,materials[2]); rows+=dim
    end
    lastlayer=((count-1)*dim+1):(count*dim)
    A[rows,lastlayer]=outer[traction,:]+Zw*outer[1,:]; rhs[rows]=-g; rows+=1
    n > 0 && (A[rows,lastlayer]=outer[4,:])
    # Equilibrate the very different displacement/traction units.
    columns=vec(maximum(abs.(A);dims=1)); A=A./transpose(columns)
    scales=vec(maximum(abs.(A);dims=2)); d=(A./scales)\(rhs./scales); d=d./columns
    U=sum(outer[1,:].*d[lastlayer])
    (s.water_density*ω^2*U-dp0)/dH
end

@testset "air cavity resonances are not hidden by coarse sampling" begin
    s=CoatedCylinderExperiment(); seeds=coated_cavity_resonances(s,5000,60000)
    @test length(seeds) == 3
    @test seeds[1] ≈ 26897.2417617 atol=1e-6
    @test_throws ArgumentError coated_cavity_resonances(s,0,60000)
    for f in (seeds[1],seeds[1]+1.5,seeds[1]+5.5), case in (:bare,:coated)
        spec=remake_experiment(s;frequency=f,scenario=case); r=simulate(spec)
        exact=[analytic_coated(spec,n) for n in 0:5]
        @test norm(r.coefficients[1:6]-exact)/norm(exact) < 1e-4
        @test r.diagnostics["power_balance_relative"] < 1e-5
    end
    scan=coated_frequency_sweep(s,5000:5000:60000;resolve_cavity=true)
    @test length(scan.frequencies) > 12
    @test all(diff(scan.frequencies).>0)
    @test scan.cavity_frequencies == seeds
    @test any(scan.ratio .> 1) # no unconditional broadband-cloak claim
    @test any(abs.(scan.frequencies.-(seeds[1]+1.5)).<1e-6)
end

@testset "coated cylinder API and invariants" begin
    s=CoatedCylinderExperiment()
    @test validate_experiment(s) === s
    @test_throws ArgumentError remake_experiment(s;inner_radius=s.outer_radius)
    @test_throws ArgumentError remake_experiment(s;observation_radius=.0049)
    @test_throws ArgumentError remake_experiment(s;frequency=200e3,coating_elements=2)
    @test_throws ArgumentError remake_experiment(s;coating_pressure_speed=20)
    @test_throws ArgumentError remake_experiment(s;incident_pressure=0)
    @test_throws ArgumentError coated_frequency_sweep(s,[20000,10000])
    mktemp() do path,io
        close(io); save_experiment(path,s)
        @test experiment_toml(load_experiment(path)) == experiment_toml(s)
    end
    empty,bare,coated=coated_comparison(s)
    @test empty.scattering_width == 0
    @test isempty(empty.radial_nodes)
    @test coated_pressure(empty,(0.,0.)) == 1
    @test coated_pressure(empty,(0.,0.);scattered=true) == 0
    @test_throws ArgumentError coated_pressure(coated,(.0045,0.))
    @test_throws ArgumentError coated_displacement(coated,(0.,0.))
    @test bare.diagnostics["steel_average_density"] ≈ s.water_density
    @test bare.scattering_width > 0
    @test 0 < coated.scattering_width/bare.scattering_width < .02
    @test coated.diagnostics["scaled_residual"] < 1e-9
    @test coated.diagnostics["modal_tail_fraction"] < 1e-8
    @test maximum(abs.(coated.displacement_tangential[:,1])) == 0
    @test maximum(abs.(coated.displacement_radial)) > 0
    double=simulate(remake_experiment(s;incident_pressure=2))
    @test double.coefficients ≈ 2coated.coefficients
    @test double.scattering_width ≈ coated.scattering_width
    @test double.scattered_pressure ≈ 2coated.scattered_pressure
    rotated=simulate(remake_experiment(s;incident_angle=π/2))
    @test coated_pressure(rotated,(0.,.020)) ≈ coated_pressure(coated,(.020,0.))
    @test rotated.scattering_width ≈ coated.scattering_width
    @test coated_displacement(rotated,(0.,s.outer_radius))[2] ≈ coated_displacement(coated,(s.outer_radius,0.))[1]
    outer=simulate(remake_experiment(s;observation_radius=.060))
    @test outer.scattering_width == coated.scattering_width
    @test outer.scattered_pressure != coated.scattered_pressure
    # Far-field flux defines a genuine width, independent of a sensor radius.
    dθ=2π/(length(coated.angles)-1)
    integral=dθ*(sum(abs2,coated.far_field_amplitude)-abs2(coated.far_field_amplitude[1]))/s.incident_pressure^2
    @test integral ≈ coated.scattering_width rtol=1e-12
    for radius in (.010,.030,.060)
        k=2π*s.frequency/s.water_speed
        flux=0.0
        for θ in coated.angles[1:end-1]
            p=sum(coated.coefficients[n+1]*hankelh1(n,k*radius)*cos(n*θ) for n in 0:s.modes)
            dp=sum(coated.coefficients[n+1]*k*(hankelh1(n-1,k*radius)-hankelh1(n+1,k*radius))/2*cos(n*θ) for n in 0:s.modes)
            flux+=imag(conj(p)*dp)*radius*dθ/k/s.incident_pressure^2
        end
        @test flux ≈ coated.scattering_width rtol=1e-10
    end
    lossless=simulate(remake_experiment(s;coating_shear_loss=0))
    @test lossless.scattering_width/bare.scattering_width < .02
    @test abs(lossless.diagnostics["absorption_width_m"]) < 1e-12
    @test lossless.diagnostics["power_balance_relative"] < 1e-5
    for case in (:bare,:coated)
        zero_poisson=simulate(remake_experiment(s;scenario=case,steel_poisson=0))
        exact=[analytic_coated(zero_poisson.experiment,n) for n in 0:5]
        @test all(isfinite,zero_poisson.coefficients)
        @test norm(zero_poisson.coefficients[1:6]-exact)/norm(exact) < 1e-4
    end
end

@testset "independent Bessel reference, radial and modal convergence" begin
    for f in (5000.,10000.,20000.,40000.,60000.), case in (:bare,:coated)
        s=CoatedCylinderExperiment(;frequency=f,scenario=case)
        r=simulate(s)
        exact=[analytic_coated(s,n) for n in 0:5]
        @test norm(r.coefficients[1:6]-exact)/norm(exact) < 1e-4
        @test r.diagnostics["absorption_width_m"] >= -1e-10
        @test r.diagnostics["power_balance_relative"] < 1e-5
        refined=simulate(remake_experiment(s;steel_elements=24,coating_elements=96))
        @test abs(refined.scattering_width/r.scattering_width-1) < 1e-4
        more=simulate(remake_experiment(s;modes=16))
        @test more.scattering_width ≈ r.scattering_width rtol=1e-8
    end
    # Changing thickness can spoil cancellation: not a hard-coded gain.
    s=CoatedCylinderExperiment()
    nominal=simulate(s).scattering_width
    @test simulate(remake_experiment(s;coating_thickness=.0005)).scattering_width > 5nominal
    @test simulate(remake_experiment(s;coating_thickness=.0015)).scattering_width > 5nominal
    sweep=coated_frequency_sweep(s,[5000,10000,20000,40000,60000])
    @test all(sweep.coated_width .> 0)
    @test all(sweep.ratio .< 1)
    @test sweep.ratio[end] > sweep.ratio[1]
    @test length(split(chomp(coated_sweep_csv(sweep)),'\n')) == 6
end
