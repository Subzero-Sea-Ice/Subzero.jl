# Shared by bench_main.jl and bench_branch.jl: load the floes saved by make_floes.jl.
using Subzero, JLD2, StructArrays, BenchmarkTools, Logging, Printf

const INFILE = get(ENV, "BENCH_FLOES", "floes_10k.jld2")
# Numbers of floes to benchmark. Larger numbers than in INFILE are made by tiling copies
# of the floe field next to each other in x.
const SIZES = parse.(Int, split(get(ENV, "BENCH_SIZES", "1000,10000,100000"), ","))
const SECONDS = parse(Float64, get(ENV, "BENCH_SECONDS", "10"))

function load_floes(n; file = INFILE)
    d = load(file)
    nsaved = length(d["coords"])
    Lx = d["Lx"]
    floes = Floe{Float64}[]
    for k in 0:n-1
        i = k % nsaved + 1
        shift = (k ÷ nsaved) * Lx  # x-offset of this copy of the floe field
        inters = copy(d["interactions"][i])
        inters[:, Int(Subzero.xpoint)] .+= shift
        push!(floes, Floe{Float64}(;
            poly = make_polygon([[[x + shift, y] for (x, y) in d["coords"][i]]]),
            centroid = d["centroid"][i] .+ [shift, 0],
            x_subfloe_points = Float64[], y_subfloe_points = Float64[],
            interactions = inters,
            (f => copy(d[string(f)][i]) for f in (
                :height, :area, :mass, :rmax, :moment, :angles,
                :α, :u, :v, :ξ, :fxOA, :fyOA, :trqOA, :hflx_factor, :overarea,
                :collision_force, :collision_trq, :num_inters,
                :stress_accum, :stress_instant, :strain, :damage,
                :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt,
            ))...,
        ))
    end
    return StructArray(floes), d["Δt"]
end

# Same stress calculator, max_floe_height, and maximum_ξ as used by make_floes.jl
floe_settings(FT = Float64) = FloeSettings(FT)

function report(label, n, trial)
    t = median(trial).time  # ns
    @printf("%-45s n=%-7d median %10.3f ms  min %10.3f ms  (%.1f ns/floe, %d samples)\n",
        label, n, t / 1e6, minimum(trial).time / 1e6, t / n, length(trial))
    flush(stdout)
end

# Silence the per-event @info logging of timestep_floe_properties! on main, so that the
# benchmark measures the calculation. (On the branch, logging happens outside the function.)
global_logger(NullLogger())

# Fields to compare between main and the branch after one call of timestep_floe_properties!
const CHECK_FIELDS = (:centroid, :height, :mass, :moment, :α, :u, :v, :ξ,
    :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt, :stress_accum, :stress_instant, :strain)
result_fields(floes) = Dict(string(f) => collect(getproperty(floes, f)) for f in CHECK_FIELDS)
