module SubzeroMakieExt

using CairoMakie, GeoInterfaceMakie
using Subzero, JLD2
import Subzero: prettytime, plot_sim

# Midpoint of edge k of a ring (ring[k] -> ring[k+1])
@inline function _bond_edge_mid(ring, k)
    return ((ring[k][1] + ring[k+1][1]) / 2, (ring[k][2] + ring[k+1][2]) / 2)
end

#=
For one saved timestep, build the geometry of all bonds. A bond sits on the edge the two floes share.
  Active bonds (bonds[i][k] > 0):
    - `segs`: pairs of points; each consecutive pair is the segment between the midpoint of the bonded
              edge on floe i and the midpoint of the matching edge on floe j. For an unstretched bond
              the two midpoints coincide (zero length); the segment appears as the floes move apart/slip.
    - `mids`: the center of each segment, drawn as a dot, so every bond is visible from the start.
  Broken bonds (bonds[i][k] = -1 tension, -2 compression, -3 shear):
    - `broken[mode]`: midpoint of the broken edge on every floe that has one. The neighbour id is no
      longer stored once a bond is broken, so each of the two floes gets its own marker; as the floes
      separate, the two markers separate with them.
Each active bond is drawn once. Floes can have ghost copies (periodic boundaries), so every
parent/ghost pair is tried and only pairs closer than `maxlen` are drawn: a bond across a periodic
boundary shows up as short segments next to the seam and never as one line across the domain.
=#
function _bond_geometry(coords, bonds, ids, ghost_ids, maxlen)
    copies = Dict{Int, Vector{Int}}()               # id -> indices of parent + ghosts
    parent = Dict{Int, Int}()                       # id -> index of the parent (ghost_id == 0)
    for (i, id) in enumerate(ids)
        push!(get!(copies, id, Int[]), i)
        ghost_ids[i] == 0 && (parent[id] = i)
    end

    segs = Point2f[]
    mids = Point2f[]
    for (id_i, ia) in copies
        pi_ = get(parent, id_i, 0)
        pi_ == 0 && continue
        for (k, nid) in enumerate(bonds[pi_])
            (nid <= 0 || nid < id_i) && continue    # active bonds only, each bond once
            ib = get(copies, nid, nothing)
            pj = get(parent, nid, 0)
            (ib === nothing || pj == 0) && continue
            kn = findfirst(==(id_i), bonds[pj])     # edge of floe j that points back to i
            kn === nothing && continue
            for a in ia, b in ib
                mi = _bond_edge_mid(coords[a][1], k)
                mj = _bond_edge_mid(coords[b][1], kn)
                hypot(mi[1] - mj[1], mi[2] - mj[2]) > maxlen && continue
                push!(segs, Point2f(mi[1], mi[2]), Point2f(mj[1], mj[2]))
                push!(mids, Point2f((mi[1] + mj[1]) / 2, (mi[2] + mj[2]) / 2))
            end
        end
    end

    broken = (Point2f[], Point2f[], Point2f[])      # tension, compression, shear
    for i in eachindex(bonds)
        ring = coords[i][1]
        for (k, b) in enumerate(bonds[i])
            -3 <= b <= -1 || continue
            m = _bond_edge_mid(ring, k)
            push!(broken[-b], Point2f(m[1], m[2]))
        end
    end
    return segs, mids, broken
end

# Tick labels in km for an axis whose data are in meters (20000.0 -> "20", 2500.0 -> "2.5")
function _km_ticks(values)
    return map(values) do v
        km = round(v / 1000; digits = 3)
        isinteger(km) ? string(Int(km)) : string(km)
    end
end

"""
    plot_sim(floe_fn, initial_state_fn, Δt, output_fn; plot_bonds = true, plot_broken_bonds = true,
        bond_color = :red, broken_bond_colors, bond_legend = true,
        ocean_color, floe_color = :white, px_per_unit = 2)

Basic plotting of a simulation using the simulation's floe and initial state files. This
function is meant for basic plotting during testing and as an example of how to create
a video using Makie. The user can (and should!) write their own plotting code to their own
specifications and needs.

This code does not have an underlying ocean, but there are comments on how to add an ocean
heatmap and colorbar.

If the floe file contains bonds, they are drawn on top of the floes. An active bond is a segment
between the midpoints of the two bonded edges plus a dot at its center: at the start it is just a dot,
and it turns into a line as the floes move relative to each other. Broken bonds are marked with a
cross on the broken edge, colored by how they broke (tension, compression or shear).

## Arguments:
- `floe_fn::String`: $(Subzero.FLOE_FN_DEF)
- `initial_state_fn::String`: $(Subzero.INITIAL_STATE_FN_DEF)
- $(Subzero.ΔT_DEF)
- `output_fn::String`: $(Subzero.MP4_OUTPUT_FN)

## Keyword arguments:
- `plot_bonds::Bool`: draw the bonds if the floe file has them (default `true`)
- `bond_color`: color of the bonds (default `:red`)
- `plot_broken_bonds::Bool`: mark broken bonds (default `true`)
- `broken_bond_colors`: colors for broken bonds, default
  `(tension = :orange, compression = :dodgerblue, shear = :magenta)`
- `bond_legend::Bool`: add a legend for the bond colors (default `true`)
- `ocean_color`: color of the ocean (axis background), default a medium blue
- `floe_color`: color of the ice floes, default white
- `px_per_unit`: resolution multiplier of the video frames (default 2; raise for sharper output)

The axes are labelled in km (the simulation data are in meters) and the text uses the
Computer Modern family (Makie's `theme_latexfonts`), applied only while this figure is made.
"""
function plot_sim(
    floe_fn,
    initial_state_fn,
    Δt,
    output_fn;
    plot_bonds = true,
    plot_broken_bonds = true,
    bond_color = :red,
    broken_bond_colors = (tension = :orange, compression = :dodgerblue, shear = :magenta),
    bond_legend = true,
    ocean_color = RGBf(30/255, 100/255, 180/255),
    floe_color = :white,
    px_per_unit = 2,
)
    # Color Information
    topo_color = RGBf(147/255, 129/255, 102/255)
    # Domain Information
    domain = load(initial_state_fn)["sim"].model.domain
    xmax, xmin = domain.east.val, domain.west.val
    ymax, ymin = domain.north.val, domain.south.val
    Δx, Δy = xmax - xmin, ymax - ymin
    # Floe Information
    file = jldopen(floe_fn)
    sim_polys = file["poly"]
    timesteps = keys(sim_polys)
    # Set up observables for recording (updated whenever `time[]` is set to a new value)
    time = Observable(timesteps[1])  # note these are strings as we index into a JLD2 file
    time_polys = @lift(sim_polys[$time])
    title_string = @lift(Subzero.prettytime(parse(Float64, $time) * Δt))

    #=
    Note that if you had Nx by Ny by time ocean data, this would be a good place to connect
    it to the time observable! To learn more about Observables and how to use them with
    animations see here: https://docs.makie.org/stable/explanations/animation
    =#

    # Fonts and sizes (Computer Modern text; applied only while this figure is built and recorded)
    plot_theme = merge(
        Theme(
            fontsize = 14,
            Axis = (titlesize = 15, xlabelsize = 14, ylabelsize = 14, xticklabelsize = 14, yticklabelsize = 14),
        ),
        theme_latexfonts(),
    )

    with_theme(plot_theme) do
        # Set up figure
        fig = Figure()
        ax = Axis(
            fig[1, 1];
            backgroundcolor = ocean_color,
            limits = (xmin, xmax, ymin, ymax),
            title = title_string,
            xlabel = "x (km)",
            xtickformat = _km_ticks,
            ylabel = "y (km)",
            ytickformat = _km_ticks,
        )

        # Plot starting state (floes + topography)
        poly!(time_polys, color = floe_color, strokecolor = :black, strokewidth = 0.5)  # floes
        if !isempty(domain.topography)  # topography
            topo_mp = domain.topography.poly
            poly!(topo_mp, color = topo_color)
        end

        # Plot bonds on top of the floes (only if the output file has them)
        if plot_bonds && all(haskey(file, key) for key in ("bonds", "coords", "id", "ghost_id"))
            sim_coords = file["coords"]
            sim_bonds = file["bonds"]
            sim_ids = file["id"]
            sim_ghost_ids = file["ghost_id"]
            maxlen = 0.5 * min(Δx, Δy)  # longer "bonds" are periodic-boundary wraps, not real bonds
            bond_data = @lift begin
                bonds_t = sim_bonds[$time]
                _bond_geometry(sim_coords[$time], bonds_t, sim_ids[$time], sim_ghost_ids[$time], maxlen)
            end
            bond_segs = @lift($bond_data[1])
            bond_mids = @lift($bond_data[2])
            linesegments!(ax, bond_segs, color = bond_color, linewidth = 2, label = "Bond")
            scatter!(ax, bond_mids, color = bond_color, markersize = 4)
            if plot_broken_bonds
                mode_names = (:tension, :compression, :shear)
                mode_labels = ("Broken: tension", "Broken: compression", "Broken: shear")
                for m in 1:3
                    broken_pts = @lift($bond_data[3][m])
                    scatter!(
                        ax, broken_pts;
                        color = getproperty(broken_bond_colors, mode_names[m]),
                        marker = :xcross, markersize = 8, label = mode_labels[m],
                    )
                end
            end
            if bond_legend
                axislegend(ax; position = :rt, labelcolor = :white, bgcolor = (:black, 0.5), framevisible = false)
            end
        end

        #=
        Note that if you had Nx by Ny by time ocean data, this would be a good place to plot it
        with a heat map and to add a color bar on the fig[1,2] axis! Check out Makie
        documentation here: https://docs.makie.org/stable/reference/plots/heatmap 
        =#

        # Resize Figure
        colsize!(fig.layout, 1, Aspect(1, Δx / Δy))
        resize_to_layout!(fig)

        # Create movie
        record(fig, output_fn, timesteps; framerate = 20, px_per_unit) do t
            time[] = t  # updating the time auto updates the related time_polys and bonds
        end
    end

    # End movie and close file
    close(file)
end

end  # module