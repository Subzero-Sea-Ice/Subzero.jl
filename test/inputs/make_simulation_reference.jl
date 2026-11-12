#=
Make the reference results of test_simulation_reference.jl with the original CPU
implementation. Run it with Subzero at ae96ff2 (before the GPU port), e.g. in the
env-baseline of the benchmarks (see benchmark/setup.sh), which also has JLD2, GeoInterface
and StructArrays:

    julia --project=<work dir>/env-baseline test/inputs/make_simulation_reference.jl

It makes a new floe field, saves it as plain arrays, runs the simulation from the saved
floes, and saves the results, all in simulation_reference.jld2 next to this script. A new
floe field is different every time (the Voronoi tesselation isn't reproducible), so rerun it
only to make new reference results.
=#
using Subzero, JLD2, Logging, Random, Statistics, StructArrays
import GeoInterface as GI
include(joinpath(@__DIR__, "simulation_reference_setup.jl"))

grid = reference_grid()
generator = VoronoiTesselationFieldGenerator(; nfloes = 10, concentrations = [0.75], hmean = 0.25, Δh = 0.0)
# Before Subzero passes the rng to GO.voronoi, the Voronoi tesselation uses the global rng
Random.seed!(1)
floes = initialize_floe_field(Float64; generator, domain = reference_domain(grid), rng = Xoshiro(1),
    floe_settings = reference_floe_settings(grid))
floe_dict = reference_floes_to_dict(floes)
results = run_reference_simulation(reference_simulation(reference_floes_from_dict(floe_dict)))
commit = readchomp(`git -C $(pkgdir(Subzero)) rev-parse --short HEAD`)
jldsave(joinpath(@__DIR__, "simulation_reference.jld2");
    commit, floes = floe_dict, results = Dict(string(k) => v for (k, v) in results))
final = results[maximum(REFERENCE_TIMESTEPS)]
@info "Saved reference results" commit length(floes) length(final["id"])
