# Results files (results/NN_<name>.toml), extrapolation, and the correctness check.

# Results of different machines can't be compared, so every machine has its own folder
const MACHINE = get(ENV, "BENCH_MACHINE", gethostname())
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", MACHINE)
# Results of the baseline for the correctness check, in the work directory (not committed)
const CHECK_DIR = "check"

results_file(step) = joinpath(RESULTS_DIR, step.name * ".toml")
read_results(step) = isfile(results_file(step)) ? TOML.parsefile(results_file(step)) : Dict{String, Any}()

# Commit of the Subzero that is loaded, with "-dirty" if its src/ or Project.toml has
# uncommitted changes
function subzero_commit()
    dir = pkgdir(Subzero)
    try
        commit = readchomp(`git -C $dir rev-parse --short=7 HEAD`)
        dirty = !isempty(readchomp(`git -C $dir status --porcelain -- src Project.toml`))
        dirty ? commit * "-dirty" : commit
    catch
        "unknown"
    end
end

cpu_name() = strip(replace(Sys.cpu_info()[1].model, r"\s+" => " "))
# The CPUs Julia runs on: set by run.sh with taskset, or else the CPUs the process may use
# (e.g. those of a Slurm job)
function cpu_set()
    s = get(ENV, "BENCH_CPUS", "")
    isempty(s) || return s
    line = isfile("/proc/self/status") ? filter(startswith("Cpus_allowed_list"), readlines("/proc/self/status")) : []
    isempty(line) ? "all" : split(only(line))[2]
end
# Whether hyperthreading (SMT) is on, and the number of CPU sockets. Linux only.
function smt()
    f = "/sys/devices/system/cpu/smt/active"
    isfile(f) ? readchomp(f) == "1" : "unknown"
end
# The CPUs this process may use (e.g. those of a Slurm job, or set with taskset)
function allowed_cpus()
    isfile("/proc/self/status") || return Int[]
    line = filter(startswith("Cpus_allowed_list"), readlines("/proc/self/status"))
    isempty(line) && return Int[]
    reduce(vcat, [length(r) == 1 ? r : r[1]:r[2]
        for r in (parse.(Int, split(x, "-")) for x in split(split(only(line))[2], ","))]; init = Int[])
end
# The number of sockets of those CPUs
function sockets()
    ids = Set{String}()
    for cpu in allowed_cpus()
        f = "/sys/devices/system/cpu/cpu$cpu/topology/physical_package_id"
        isfile(f) && push!(ids, readchomp(f))
    end
    isempty(ids) ? "unknown" : length(ids)
end

# Round to 4 significant digits, so the results files stay readable
round4(x) = round(x; sigdigits = 4)

"""
    extrapolate(records, sizes)

Fit `t = a nᵇ` to the times of each measurement (label, backend, float type) in `records`
and return records for `sizes`, marked as extrapolated, with the exponent `b`.
"""
function extrapolate(records, sizes)
    out = Dict{String, Any}[]
    groups = unique(r -> (r["label"], r["backend"], r["FT"]), records)
    for g in groups
        rs = filter(r -> (r["label"], r["backend"], r["FT"]) == (g["label"], g["backend"], g["FT"]), records)
        length(rs) >= 2 || continue
        logn = log.(getindex.(rs, "n"))
        fit(t) = begin  # least squares on log t = log a + b log n
            X = [ones(length(logn)) logn]
            X \ log.(t)
        end
        loga_med, b = fit(getindex.(rs, "median_ms"))
        loga_min, b_min = fit(getindex.(rs, "min_ms"))
        for n in sizes
            push!(out, merge(g, Dict{String, Any}(
                "n" => n, "median_ms" => exp(loga_med) * n^b, "min_ms" => exp(loga_min) * n^b_min,
                "samples" => 0, "extrapolated" => true, "exponent" => round(b; digits = 2),
                "fitted_sizes" => getindex.(rs, "n"),
            )))
            @printf("%-8s %-6s %-7s %-40s n=%-7d median ≈%9.3f ms  (extrapolated, t ∝ n^%.2f)\n",
                "", g["backend"], g["FT"], g["label"], n, exp(loga_med) * n^b, b)
        end
    end
    return out
end

"""
    save_results!(step, env, records, checks)

Merge `records` (timings) and `checks` into the step's results file. Earlier results of the
same environment, number of threads, backend and number of floes are replaced, everything
else is kept.
"""
function save_results!(step, env, records, checks)
    d = read_results(step)
    d["title"] = step.title
    env == "current" && (d["ported"] = is_ported(step))
    get!(d, "ported", false)
    threads = Threads.nthreads()
    common = Dict{String, Any}(
        "env" => env, "threads" => threads, "commit" => subzero_commit(),
        "date" => string(Dates.today()), "cpu" => cpu_name(), "cpus" => cpu_set(),
        "smt" => smt(), "sockets" => sockets(), "input" => input_checksum(),
    )
    new = [merge(r, common, Dict("device" => device_name(r["backend"]))) for r in records]
    for r in new
        r["median_ms"] = round4(r["median_ms"])
        r["min_ms"] = round4(r["min_ms"])
        get!(r, "extrapolated", false)
    end
    key(r) = (r["env"], r["threads"], r["backend"], r["n"])
    replaced = Set(key.(new))
    timings = filter(r -> key(r) ∉ replaced, get(d, "timing", Dict{String, Any}[]))
    append!(timings, new)
    sort!(timings; by = r -> (r["env"], r["backend"], r["FT"], r["threads"], r["n"], r["label"]))
    d["timing"] = timings
    if !isempty(checks)
        checks = [merge(c, Dict{String, Any}("commit" => common["commit"], "date" => common["date"],
            "threads" => threads)) for c in checks]
        ckey(c) = (c["backend"], get(c, "threads", 0), c["n"])
        replaced = Set(ckey.(checks))
        all_checks = filter(c -> ckey(c) ∉ replaced, get(d, "check", Dict{String, Any}[]))
        append!(all_checks, checks)
        sort!(all_checks; by = ckey)
        d["check"] = all_checks
    end
    mkpath(RESULTS_DIR)
    open(io -> TOML.print(io, d; sorted = true), results_file(step), "w")
    return
end

# True if the results file has the baseline's primary timing for every size in SIZES with
# this number of threads
function has_baseline(step)
    timings = get(read_results(step), "timing", [])
    all(SIZES) do n
        any(r -> r["env"] == "baseline" && r["threads"] == Threads.nthreads() && r["n"] == n &&
            r["primary"], timings)
    end
end

# Flatten a field (a vector of numbers or of arrays) to a vector of numbers
flat(x) = eltype(x) <: Number ? vec(collect(x)) : reduce(vcat, vec.(x); init = Float64[])

"""
    compare(result, reference; rtol)

Maximum difference between two results, relative to the largest absolute value of each
field, and whether every field is `isapprox` with `rtol`.
"""
function compare(result, reference; rtol = sqrt(eps()))
    keys(result) == keys(reference) || return (Inf, false)
    maxrel = 0.0
    ok = true
    for k in keys(reference)
        a, b = flat(result[k]), flat(reference[k])
        if length(a) != length(b)
            @printf("  %s: %d values, baseline %d\n", k, length(a), length(b))
            return (Inf, false)
        end
        isempty(b) && continue
        maxrel = max(maxrel, maximum(abs.(a .- b)) / max(maximum(abs.(b)), eps()))
        ok &= isapprox(a, b; rtol)
    end
    return (maxrel, ok)
end

check_file(step) = joinpath(CHECK_DIR, step.name * ".jld2")
