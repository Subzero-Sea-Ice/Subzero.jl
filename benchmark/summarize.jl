# Make the overview tables from results/<machine>/*.toml. Needs no packages:
#     julia benchmark/summarize.jl          # update the tables
#     julia benchmark/summarize.jl --check  # fail if a table isn't up to date
# Every machine gets results/<machine>/README.md, and the reference machine's table also goes
# in README.md. BENCH_SIZES (the first size) chooses the number of floes (default 100000).
using TOML, Printf

const DIR = @__DIR__
const README = joinpath(DIR, "README.md")
const START, STOP = "<!-- overview:start -->", "<!-- overview:end -->"
const N = parse(Int, first(split(get(ENV, "BENCH_SIZES", "100000"), ",")))
# The machine whose results are in README.md
const REFERENCE = "laptop-i7-12700H"
# Checksum of the input file of the reference results (see input_checksum in lib/input.jl)
const REFERENCE_INPUT = "d1805bc093a4"
# The step whose time is 100% in the "share of timestep" column
const TIMESTEP = "00_timestep_sim"

thousands(s) = replace(s, r"(\d)(?=(\d{3})+(\.|$))" => s"\1,")
function fmt(r)
    r === nothing && return "–"
    t = r["median_ms"]
    s = thousands(t >= 100 ? @sprintf("%.0f", t) : t >= 10 ? @sprintf("%.1f", t) : @sprintf("%.2f", t))
    # The exponent goes in the cell rather than in a numbered footnote, so that a change to
    # one step only changes its own row
    get(r, "extrapolated", false) ? "≈$s (n<sup>$(r["exponent"])</sup>)" : s
end
ratio(a, b) = a === nothing || b === nothing ? "–" : @sprintf("%.3g×", a["median_ms"] / b["median_ms"])
function share(a, b)
    (a === nothing || b === nothing) && return "–"
    p = 100a["median_ms"] / b["median_ms"]
    p < 0.5 ? "<1%" : @sprintf("%.0f%%", p)
end
function check(d)
    cs = get(d, "check", [])
    isempty(cs) ? "–" : all(c -> c["ok"], cs) ? "✓" : "✗"
end

"""
    table(machine, steps_dir)

The overview table of `machine` in Markdown, with links to the step notes in `steps_dir`
(relative to where the table goes).
"""
function table(machine, steps_dir)
    dir = joinpath(DIR, "results", machine)
    steps = [(splitext(f)[1], TOML.parsefile(joinpath(dir, f)))
        for f in sort(filter(endswith(".toml"), readdir(dir)))]
    timings(d) = filter(r -> r["primary"] && r["n"] == N, get(d, "timing", []))
    all_timings = reduce(vcat, timings.(last.(steps)); init = [])

    # Columns: the most threads that were measured, and every GPU
    threads = maximum(r -> r["threads"], all_timings; init = 1)
    gpu_fts = [(gpu, FT) for gpu in sort(unique(r["backend"] for r in all_timings if r["backend"] != "CPU"))
        for FT in ("Float64", "Float32") if any(r -> r["backend"] == gpu && r["FT"] == FT, all_timings)]
    function find(d, env, backend, FT, threads)
        rs = filter(r -> r["env"] == env && r["backend"] == backend && r["FT"] == FT &&
            (threads === nothing || r["threads"] == threads), timings(d))
        isempty(rs) ? nothing : rs[argmax([r["date"] for r in rs])]
    end

    s(t) = t == 1 ? "" : "s"
    header = ["step", "ported", "baseline 1 thread", "baseline $threads thread$(s(threads))",
        "current CPU 1 thread", "current CPU $threads thread$(s(threads))",
        ["current $gpu $FT" for (gpu, FT) in gpu_fts]..., "speedup", "share of timestep",
        "check", "measured"]
    align = ["---", ":-:", fill("--:", length(header) - 5)..., "--:", ":-:", "---"]
    lines = ["| " * join(header, " | ") * " |", "|" * join(align, "|") * "|"]
    timestep = get(Dict(steps), TIMESTEP, Dict())
    for (name, d) in steps
        link = isfile(joinpath(DIR, "steps", name * ".md")) ? "[$(d["title"])]($steps_dir/$name.md)" : d["title"]
        base = find(d, "baseline", "CPU", "Float64", threads)
        cur_cpu = find(d, "current", "CPU", "Float64", threads)
        cur_gpu = [find(d, "current", gpu, FT, nothing) for (gpu, FT) in gpu_fts]
        best = filter(!isnothing, [cur_cpu; [r for (r, (_, FT)) in zip(cur_gpu, gpu_fts) if FT == "Float64"]])
        best = isempty(best) ? nothing : best[argmin([r["median_ms"] for r in best])]
        current = filter(r -> r["env"] == "current", timings(d))
        measured = isempty(current) ? "–" :
            (r -> "$(r["date"]) `$(r["commit"])`")(current[argmax([r["date"] for r in current])])
        cells = [
            # The whole timestep is ported as far as its steps are
            link, name == TIMESTEP ? "partly" : d["ported"] ? "yes" : "no",
            fmt(find(d, "baseline", "CPU", "Float64", 1)), fmt(base),
            fmt(find(d, "current", "CPU", "Float64", 1)), fmt(cur_cpu),
            fmt.(cur_gpu)...,
            d["ported"] ? ratio(base, best) : "–",
            name == TIMESTEP ? "100%" : share(base, find(timestep, "baseline", "CPU", "Float64", threads)),
            check(d), measured,
        ]
        push!(lines, "| " * join(cells, " | ") * " |")
    end
    push!(lines, "", "Median time of one call in ms for $(thousands(string(N))) floes. ≈: extrapolated " *
        "from fewer floes as t ∝ nᵇ, with the exponent b.", "", machine_info(machine, all_timings))
    return join(lines, "\n")
end

# What the results of a machine were measured on, from the results files
function machine_info(machine, timings)
    isempty(timings) && return "Machine `$machine`: no results for $(thousands(string(N))) floes."
    values(key) = join(sort(unique(string(get(r, key, "unknown")) for r in timings)), ", ")
    cpu = filter(r -> r["backend"] == "CPU", timings)
    gpus = sort(unique(r["device"] for r in timings if r["backend"] != "CPU"))
    smt = values("smt")
    smt = smt == "true" ? "on" : smt == "false" ? "off" : smt
    input = values("input")
    input = input == REFERENCE_INPUT ? "the reference input `$input`" :
        "input `$input` (not the reference input `$REFERENCE_INPUT`, so not directly comparable)"
    "Machine `$machine`: $(values("cpu")), $(values("sockets")) socket$(values("sockets") == "1" ? "" : "s"), hyperthreading " *
        "$smt, run on CPUs $(join(sort(unique(r["cpus"] for r in cpu)), ", ")); " *
        "GPU: $(isempty(gpus) ? "none" : join(gpus, ", ")); $input."
end

machines = sort(filter(m -> isdir(joinpath(DIR, "results", m)), readdir(joinpath(DIR, "results"))))
# The new content of every file, to write or to check
files = Dict{String, String}()
for machine in machines
    files[joinpath(DIR, "results", machine, "README.md")] = """
        # Benchmark results: `$machine`

        Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
        [the benchmark README](../../README.md) for what the columns mean.

        $(table(machine, "../../steps"))
        """
end
readme = read(README, String)
i, j = findfirst(START, readme), findfirst(STOP, readme)
(i === nothing || j === nothing) && error("README.md has no $START … $STOP markers")
REFERENCE in machines || error("no results for the reference machine $REFERENCE")
files[README] = readme[1:last(i)] * "\n" * table(REFERENCE, "steps") * "\n" * readme[first(j):end]

if "--check" in ARGS
    outdated = [relpath(f, DIR) for (f, new) in files if !isfile(f) || read(f, String) != new]
    if !isempty(outdated)
        println(stderr, "Out of date: $(join(sort(outdated), ", ")). Run julia benchmark/summarize.jl")
        exit(1)
    end
    println("The overview tables are up to date")
else
    foreach(((f, new),) -> write(f, new), files)
end
