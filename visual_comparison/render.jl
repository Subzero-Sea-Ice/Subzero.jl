#=
Draw the runs of both versions of an example side by side, as an mp4:
    julia --project=env-current render.jl <runs dir> <videos dir> <example>...
Reads <runs dir>/original/<example>.jld2 and <runs dir>/current/<example>.jld2 (made by
run_example.jl) and writes <videos dir>/<example>.mp4.

Three panels: the original CPU implementation, the current implementation, and both on top of
each other (original outlines wide, current outlines narrow on top). Above them, the number
of floes in each run and the largest distance between the centroids of the same floe in
both runs. Fracture and ridging make new floes with new ids, and once they differ between
the runs, floes with the same new id aren't the same floe. So the distance is only taken
over the initial floes that are still in both runs, and the line says from when the runs
have different floes. Through periodic boundaries, the distance is the shortest one.
=#
using Subzero, CairoMakie, JLD2, Random, Statistics, StructArrays, Printf
const GI = Subzero.GI
include(joinpath(@__DIR__, "examples.jl"))

runs_dir, videos_dir = ARGS[1], ARGS[2]
names = ARGS[3:end]

# Colours of the documentation's animations, and two for the overlay
const OCEAN_COLOR = RGBf(4/255, 31/255, 74/255)
const FLOE_COLOR = RGBf(217/255, 226/255, 225/255)
const TOPO_COLOR = RGBf(147/255, 129/255, 102/255)
const ORIGINAL_COLOR = RGBf(1.0, 0.6, 0.1)
const CURRENT_COLOR = RGBf(0.3, 0.9, 1.0)

to_polygons(rings) = [Makie.Polygon(Point2d.(ring)) for ring in rings]

# Distance along an axis, the shortest one through the boundaries if they're periodic
axis_distance(a, b, L, periodic) = periodic ? min(abs(a - b), L - abs(a - b)) : abs(a - b)

#= Number of floes in each run, whether they have the same floes (ids), and the largest
centroid distance of the initial floes (ids in `initial_ids`) that are in both runs =#
function compare_frame(orig, cur, k, initial_ids, domain)
    a = Dict(zip(orig["ids"][k], orig["centroids"][k]))
    b = Dict(zip(cur["ids"][k], cur["centroids"][k]))
    Lx, Ly = domain.east.val - domain.west.val, domain.north.val - domain.south.val
    px, py = domain.east isa PeriodicBoundary, domain.north isa PeriodicBoundary
    dist(p, q) = hypot(axis_distance(p[1], q[1], Lx, px), axis_distance(p[2], q[2], Ly, py))
    common = intersect(keys(a), keys(b), initial_ids)
    max_dist = isempty(common) ? 0.0 : maximum(id -> dist(a[id], b[id]), common)
    return (; n_orig = length(a), n_cur = length(b), same_floes = keys(a) == keys(b),
        n_initial = length(common), max_dist)
end

function render(name)
    orig = load(joinpath(runs_dir, "original", "$name.jld2"))
    cur = load(joinpath(runs_dir, "current", "$name.jld2"))
    orig["tsteps"] == cur["tsteps"] || error("$name: the runs saved different timesteps")
    tsteps = orig["tsteps"]
    domain = example_setup(name).domain
    xmin, xmax = domain.west.val, domain.east.val
    ymin, ymax = domain.south.val, domain.north.val
    topography = to_polygons(ring_coords(domain.topography.poly))

    initial_ids = Set(orig["ids"][1])
    stats = [compare_frame(orig, cur, i, initial_ids, domain) for i in eachindex(tsteps)]
    first_different = findfirst(s -> !s.same_floes, stats)
    k = Observable(1)
    polys_orig = @lift(to_polygons(orig["polys"][$k]))
    polys_cur = @lift(to_polygons(cur["polys"][$k]))
    status = @lift begin
        s = stats[$k]
        floes = s.same_floes ? @sprintf("%d floes, the same in both runs", s.n_orig) :
            @sprintf("floes: %d original, %d current (different since %s)", s.n_orig,
                s.n_cur, Subzero.prettytime(tsteps[first_different] * Δt))
        @sprintf("%s  |  %s  |  largest centroid difference of the %d initial floes in both: %.3g m",
            Subzero.prettytime(tsteps[$k] * Δt), floes, s.n_initial, s.max_dist)
    end

    fig = Figure(size = (1500, 620))
    Label(fig[1, 1:3], EXAMPLE_TITLES[name]; fontsize = 22, font = :bold)
    Label(fig[2, 1:3], status; fontsize = 15)
    panel_titles = (
        "Original CPU implementation ($(orig["commit"]), $(orig["backend"]))",
        "Current implementation ($(cur["commit"]), $(cur["backend"]))",
        "Overlay: original (orange, wide) and current (cyan)",
    )
    axes = [Axis(fig[3, j]; title = panel_titles[j], backgroundcolor = OCEAN_COLOR,
        limits = (xmin, xmax, ymin, ymax), aspect = DataAspect(),
        xticklabelrotation = pi/4, yticklabelrotation = pi/4) for j in 1:3]
    poly!(axes[1], polys_orig; color = FLOE_COLOR, strokecolor = :black, strokewidth = 0.5)
    poly!(axes[2], polys_cur; color = FLOE_COLOR, strokecolor = :black, strokewidth = 0.5)
    poly!(axes[3], polys_orig; color = :transparent, strokecolor = ORIGINAL_COLOR,
        strokewidth = 3)
    poly!(axes[3], polys_cur; color = :transparent, strokecolor = CURRENT_COLOR,
        strokewidth = 1)
    if !isempty(topography)
        for ax in axes
            poly!(ax, topography; color = TOPO_COLOR)
        end
    end
    Label(fig[4, 1:3],
        "Both runs start from the same floes. With collisions, round-off differences grow " *
        "over time, so floes can drift apart late in a run without anything being wrong. " *
        "Look for floes that jump, vanish, overlap, or stop moving in only one run.";
        fontsize = 13, color = :gray30)

    mkpath(videos_dir)
    file = joinpath(videos_dir, "$name.mp4")
    record(fig, file, eachindex(tsteps); framerate = 10) do i
        k[] = i
    end
    s = stats[end]
    different = isnothing(first_different) ? "the same floes" :
        "different floes since timestep $(tsteps[first_different])"
    println("$name: $(length(tsteps)) frames to $file; at the end $(s.n_orig) and $(s.n_cur) ",
        "floes ($different), ", @sprintf("largest centroid difference %.3g m", s.max_dist),
        " of $(s.n_initial) initial floes")
    return file
end

foreach(render, names)
