# Settings, the step registry, and measuring.

# Numbers of floes to benchmark. Larger numbers than in the input file are made by copying
# the floe field next to itself in x (see load_floes).
const SIZES = parse.(Int, split(get(ENV, "BENCH_SIZES", "100000"), ","))
# Number of samples per benchmark, and the maximum time per benchmark in s
const SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "20"))
const SECONDS = parse(Float64, get(ENV, "BENCH_SECONDS", "60"))
# Steps to run: comma-separated names, with or without their number (default: all)
const STEP_FILTER = filter(!isempty, split(get(ENV, "BENCH_STEPS", ""), ","))
# Rerun the baseline even if its results exist
const FORCE = get(ENV, "BENCH_FORCE", "0") == "1"
# Default number of floes for the correctness check
const CHECK_N = 1000

# Silence the per-event @info logging of the baseline, so that the benchmark measures the
# calculation
global_logger(NullLogger())

"""
    Step

A part of `timestep_sim!` to benchmark, defined in `steps/NN_<name>.jl` with `register_step!`.

- `title`: what is benchmarked, for the overview table (Markdown).
- `baseline(ctx)`: benchmarks the baseline code with `@bench` and records the trials with
  `record!`. Also used in the current environment as long as `current` isn't given.
- `current(ctx, backend, FT)`: the same for the current code on a KernelAbstractions
  backend, in Float64 and (if `float32` is true and the backend is a GPU) Float32.
- `baseline_result(n)` and `current_result(n, backend)`: a `Dict` of arrays after running
  the step once on `n` floes, to check that the current code gives the same result as the
  baseline. `current_result` defaults to `baseline_result` while the step isn't ported.
- `extrapolate_above` and `extrapolate_from`: the baseline code is too slow for more than
  `extrapolate_above` floes. For those sizes in `BENCH_SIZES`, it runs on `extrapolate_from`
  floes instead, and the time is extrapolated with a power law.
- `check_n`: number of floes for the correctness check.
- `check_rtol`: relative tolerance of the correctness check, for steps whose result depends
  on the number of threads.
"""
Base.@kwdef struct Step
    name::String
    title::String
    baseline::Function
    current::Union{Function, Nothing} = nothing
    baseline_result::Union{Function, Nothing} = nothing
    current_result::Union{Function, Nothing} = nothing
    float32::Bool = false
    extrapolate_above::Int = typemax(Int)
    extrapolate_from::Vector{Int} = Int[]
    check_n::Int = CHECK_N
    check_rtol::Float64 = sqrt(eps())
end

const STEPS = Step[]
const STEP_FILE = Ref("")
register_step!(; kwargs...) = push!(STEPS, Step(; name = STEP_FILE[], kwargs...))

is_ported(step) = step.current !== nothing

# What a step function gets: where it runs and the list its measurements go to
struct Context
    env::String      # "baseline" or "current"
    n::Int           # number of floes
    backend::String  # "CPU", "CUDA", "AMDGPU"
    records::Vector{Dict{String, Any}}
end

# Every benchmark takes SAMPLES samples (or less if it takes longer than SECONDS). The
# garbage collector runs before every sample (not timed), so that the arrays allocated by
# earlier samples don't fill up the memory.
macro bench(ex, setup = :nothing)
    esc(:(@benchmark($ex, setup = $setup, evals = 1, samples = SAMPLES, seconds = SECONDS,
        gcsample = true)))
end

"""
    record!(ctx, label, trial; FT = Float64, primary = false)

Print the median and minimum time of `trial` and keep them for the results file. Every step
has one `primary` measurement per backend and float type: the one in the overview table.
"""
function record!(ctx, label, trial; FT = Float64, primary = false)
    t = median(trial).time  # ns
    @printf("%-8s %-6s %-7s %-40s n=%-7d median %10.3f ms  min %10.3f ms  (%.1f ns/floe, %d samples)\n",
        ctx.env, ctx.backend, FT, label, ctx.n, t / 1e6, minimum(trial).time / 1e6, t / ctx.n,
        length(trial))
    flush(stdout)
    push!(ctx.records, Dict{String, Any}(
        "label" => label, "primary" => primary, "backend" => ctx.backend, "FT" => string(FT),
        "n" => ctx.n, "median_ms" => t / 1e6, "min_ms" => minimum(trial).time / 1e6,
        "samples" => length(trial),
    ))
    return
end

# Peak memory use of this Julia process
report_host_memory() = @printf("Peak host memory (max RSS): %.2f GB\n", Sys.maxrss() / 1e9)
