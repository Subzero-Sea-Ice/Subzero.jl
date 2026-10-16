# Run the benchmarks in steps/ with the baseline or the current Subzero, and write their
# results to results/NN_<name>.toml.
#
# Usage, from the work directory made by setup.sh:
#     julia -t <nthreads> --project=env-baseline path/to/bench.jl baseline
#     julia -t <nthreads> --project=env-current path/to/bench.jl current
# Run the baseline first: it saves the results for the correctness check.
# Environment variables: see README.md.
using Subzero, JLD2, StructArrays, BenchmarkTools, Logging, Printf, Random, SHA, Statistics, TOML, Dates
import GeoInterface as GI

const ENV_NAME = only(ARGS)
ENV_NAME in ("baseline", "current") || error("usage: bench.jl baseline|current")

include("lib/harness.jl")
include("lib/input.jl")
include("lib/results.jl")
if ENV_NAME == "current"
    include("lib/backends.jl")
else
    const BACKENDS = ["CPU"]
    device_name(_) = cpu_name()
    track_gpu_memory(_) = nothing
    report_gpu_memory() = nothing
end

# Every file in steps/ registers one step. Steps run in the order of their file names.
for file in sort(filter(endswith(".jl"), readdir(joinpath(@__DIR__, "steps"))))
    STEP_FILE[] = splitext(file)[1]
    include(joinpath(@__DIR__, "steps", file))
end
selected(step) = isempty(STEP_FILTER) ||
    any(s -> s == step.name || s == split(step.name, "_"; limit = 2)[2], STEP_FILTER)
for s in STEP_FILTER
    any(step -> s in (step.name, split(step.name, "_"; limit = 2)[2]), STEPS) ||
        error("unknown step $s in BENCH_STEPS, expected one of $(join(getfield.(STEPS, :name), ", "))")
end

# Check that the current code gives the same result as the baseline: the baseline saves its
# result, the current code compares with it.
function check(step)
    n = step.check_n
    if ENV_NAME == "baseline"
        step.baseline_result === nothing && return Dict{String, Any}[]
        mkpath(CHECK_DIR)
        jldsave(check_file(step); n, input = input_checksum(), result = step.baseline_result(n))
        return Dict{String, Any}[]
    end
    f = is_ported(step) ? step.current_result : (n, _) -> step.baseline_result(n)
    (f === nothing || check_input(step) === nothing) && return Dict{String, Any}[]
    if check_input(step) != input_checksum()
        println("check    skipped: the baseline's result was made with input $(check_input(step)), ",
            "not $(input_checksum()). Run the baseline first.")
        return Dict{String, Any}[]
    end
    reference = load(check_file(step), "result")
    checks = Dict{String, Any}[]
    for name in (is_ported(step) ? BACKENDS : ["CPU"])
        maxrel, ok = compare(f(n, get_backend(name)), reference; rtol = step.check_rtol)
        println("check    $name vs baseline, n=$n: all fields isapprox = $ok, ",
            "max relative difference = $maxrel")
        push!(checks, Dict{String, Any}("backend" => name, "n" => n, "ok" => ok,
            "max_rel_diff" => round(maxrel; sigdigits = 2)))
    end
    return checks
end

function run_step(step)
    println("\n== $(step.name): $(step.title)")
    if ENV_NAME == "baseline" && !FORCE && has_baseline(step)
        println("baseline results with $(Threads.nthreads()) threads on this input exist, skipping (BENCH_FORCE=1 reruns them)")
        # The current code still needs the baseline's result for the correctness check
        check_input(step) == input_checksum() || check(step)
        return
    end
    checks = check(step)
    # The baseline code runs on the CPU only. So does the current code until it's ported.
    cpu_code = ENV_NAME == "baseline" || !is_ported(step)
    targets = cpu_code ? filter(>(step.extrapolate_above), SIZES) : Int[]
    sizes = unique([setdiff(SIZES, targets); isempty(targets) ? Int[] : step.extrapolate_from])
    records = Dict{String, Any}[]
    for name in (cpu_code ? ["CPU"] : BACKENDS), n in sizes
        GC.gc()
        ctx = Context(ENV_NAME, n, name, records)
        if cpu_code
            step.baseline(ctx)
        else
            for FT in (step.float32 && name != "CPU" ? (Float64, Float32) : (Float64,))
                step.current(ctx, get_backend(name), FT)
            end
        end
        track_gpu_memory(name)
    end
    if !isempty(targets)
        append!(records, extrapolate(filter(r -> r["n"] in step.extrapolate_from, records), targets))
    end
    save_results!(step, ENV_NAME, records, checks)
end

println("$ENV_NAME, Subzero $(subzero_commit()), Julia threads = $(Threads.nthreads()), CPUs = $(cpu_set()), machine = $MACHINE")
foreach(name -> name == "CPU" || print_device(name), BACKENDS)
foreach(run_step, filter(selected, STEPS))
report_host_memory()
report_gpu_memory()
