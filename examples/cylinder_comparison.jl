using MetaNDT, Gridap, GridapGmsh, Gmsh, TOML

# Exactly the same public interface used by the Pluto client.
spec=CylinderExperiment()
destination=isempty(ARGS) ? mktempdir(prefix="meta-ndt-cylinder-") : only(ARGS)
mkpath(destination)
save_experiment(joinpath(destination,"experiment.toml"),spec)
for result in simulate_comparison(spec)
    name=String(result.experiment.scenario)
    write(joinpath(destination,"$(name)_sensors.csv"),sensor_csv(result))
    open(joinpath(destination,"$(name)_diagnostics.toml"),"w") do io
        TOML.print(io,result.diagnostics;sorted=true)
    end
    println(name,": ",result.diagnostics)
end
println("Saved experiment and pressure traces to ",destination)
