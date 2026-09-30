# Benchmark timestep_floe_properties! on main (a Threads.@threads loop over a StructArray).
# Usage: julia -t <nthreads> --project=env-main bench_main.jl
include("common.jl")

println("main, Julia threads = $(Threads.nthreads())")
fs = floe_settings()
for n in SIZES
    floes, Δt = load_floes(n)
    trial = @benchmark(
        Subzero.timestep_floe_properties!(f, 1, $Δt, $fs),
        setup = (f = deepcopy($floes)), evals = 1, seconds = SECONDS,
    )
    report("main: timestep_floe_properties!", n, trial)
end

# Save the result of one call for bench_branch.jl to compare against
floes, Δt = load_floes(SIZES[1])
Subzero.timestep_floe_properties!(floes, 1, Δt, fs)
jldsave("result_main.jld2"; n = SIZES[1], result = result_fields(floes))
