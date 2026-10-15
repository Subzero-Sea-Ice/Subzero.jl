# Make the overview table in README.md from results/*.toml. Needs no packages:
#     julia benchmark/summarize.jl          # update the table in README.md
#     julia benchmark/summarize.jl --check  # fail if the table isn't up to date
# BENCH_SIZES (the first size) chooses the number of floes of the table (default 100000).
using TOML, Printf

const DIR = @__DIR__
const README = joinpath(DIR, "README.md")
const START, STOP = "<!-- overview:start -->", "<!-- overview:end -->"
const N = parse(Int, first(split(get(ENV, "BENCH_SIZES", "100000"), ",")))

steps = [(splitext(f)[1], TOML.parsefile(joinpath(DIR, "results", f)))
    for f in sort(filter(endswith(".toml"), readdir(joinpath(DIR, "results"))))]
timings(d) = filter(r -> r["primary"] && r["n"] == N, get(d, "timing", []))
all_timings = reduce(vcat, timings.(last.(steps)); init = [])

# Columns: the most threads that were measured, and every GPU
const THREADS = maximum(r -> r["threads"], all_timings; init = 1)
const GPUS = sort(unique(r["backend"] for r in all_timings if r["backend"] != "CPU"))
const GPU_FTS = [(gpu, FT) for gpu in GPUS for FT in ("Float64", "Float32")
    if any(r -> r["backend"] == gpu && r["FT"] == FT, all_timings)]

function find(d, env, backend, FT, threads)
    rs = filter(r -> r["env"] == env && r["backend"] == backend && r["FT"] == FT &&
        (threads === nothing || r["threads"] == threads), timings(d))
    isempty(rs) ? nothing : rs[argmax([r["date"] for r in rs])]
end

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

function table()
    timestep = isempty(steps) ? Dict() : last(first(steps))
    header = ["step", "ported", "baseline 1 thread", "baseline $THREADS thread$(THREADS == 1 ? "" : "s")",
        "current CPU 1 thread", "current CPU $THREADS thread$(THREADS == 1 ? "" : "s")",
        ["current $gpu $FT" for (gpu, FT) in GPU_FTS]..., "speedup", "share of timestep",
        "check", "measured"]
    align = ["---", ":-:", fill("--:", length(header) - 5)..., "--:", ":-:", "---"]
    lines = ["| " * join(header, " | ") * " |", "|" * join(align, "|") * "|"]
    for (name, d) in steps
        link = isfile(joinpath(DIR, "steps", name * ".md")) ? "[$(d["title"])](steps/$name.md)" : d["title"]
        base = find(d, "baseline", "CPU", "Float64", THREADS)
        cur_cpu = find(d, "current", "CPU", "Float64", THREADS)
        cur_gpu = [find(d, "current", gpu, FT, nothing) for (gpu, FT) in GPU_FTS]
        best = filter(!isnothing, [cur_cpu; [r for (r, (_, FT)) in zip(cur_gpu, GPU_FTS) if FT == "Float64"]])
        best = isempty(best) ? nothing : best[argmin([r["median_ms"] for r in best])]
        current = filter(r -> r["env"] == "current", timings(d))
        measured = isempty(current) ? "–" :
            (r -> "$(r["date"]) `$(r["commit"])`")(current[argmax([r["date"] for r in current])])
        cells = [
            # The whole timestep is ported as far as its steps are
            link, name == first(first(steps)) ? "partly" : d["ported"] ? "yes" : "no",
            fmt(find(d, "baseline", "CPU", "Float64", 1)), fmt(base),
            fmt(find(d, "current", "CPU", "Float64", 1)), fmt(cur_cpu),
            fmt.(cur_gpu)...,
            d["ported"] ? ratio(base, best) : "–",
            name == first(first(steps)) ? "100%" : share(base, find(timestep, "baseline", "CPU", "Float64", THREADS)),
            check(d), measured,
        ]
        push!(lines, "| " * join(cells, " | ") * " |")
    end
    push!(lines, "", "Median time of one call in ms for $(thousands(string(N))) floes. ≈: extrapolated " *
        "from fewer floes as t ∝ nᵇ, with the exponent b.")
    return join(lines, "\n")
end

readme = read(README, String)
i, j = findfirst(START, readme), findfirst(STOP, readme)
(i === nothing || j === nothing) && error("README.md has no $START … $STOP markers")
new = readme[1:last(i)] * "\n" * table() * "\n" * readme[first(j):end]
if "--check" in ARGS
    new == readme || (println(stderr, "The overview table in README.md is out of date: run julia benchmark/summarize.jl"); exit(1))
    println("The overview table in README.md is up to date")
else
    write(README, new)
end
