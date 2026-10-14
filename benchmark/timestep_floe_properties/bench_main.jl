# Benchmark timestep_floe_properties! on main (a Threads.@threads loop over a StructArray).
# Usage: julia -t <nthreads> --project=env-main bench_main.jl
include("common.jl")

println("main, Julia threads = $(Threads.nthreads())")
fs = floe_settings()
for n in SIZES
    floes, Δt = load_floes(n)
    # The floes contain nested arrays that the function changes in place, so every sample
    # starts from a deepcopy. The garbage collector runs before every sample (not timed), so
    # that the copies don't fill up the memory.
    trial = @benchmark(
        Subzero.timestep_floe_properties!(f, 1, $Δt, $fs),
        setup = (f = deepcopy($floes)), evals = 1, samples = SAMPLES, seconds = SECONDS,
        gcsample = true,
    )
    report("main: timestep_floe_properties!", n, trial)
end
report_host_memory()

# Save the result of one call for bench_branch.jl to compare against
floes, Δt = load_floes(CHECK_N)
Subzero.timestep_floe_properties!(floes, 1, Δt, fs)
jldsave("result_main.jld2"; n = CHECK_N, result = result_fields(floes))
