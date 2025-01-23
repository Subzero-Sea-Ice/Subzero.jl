module SubzeroMakieExt

using CairoMakie, GeoInterfaceMakie
using Subzero, JLD2
import Subzero: prettytime, plot_sim

"""
    plot_sim(floe_fn, initial_state_fn, Δt, output_fn; max_side_pixels = 800)

Basic plotting of a simulation using the simulation's floe and initial state files. This
function is meant for basic plotting during testing and as an example of how to create
a video using Makie. The user can (and should!) write their own plotting code to their own
specifications and needs.

This code does not have an underlying ocean, but there are comments on how to add an ocean
heatmap and colorbar.

## Arguments:
- `floe_fn::String`: $(Subzero.FLOE_FN_DEF)
- `initial_state_fn::String`: $(Subzero.INITIAL_STATE_FN_DEF)
- `Δt::Int`: $(Subzero.ΔT_DEF)
- `output_fn::String`: $(Subzero.MP4_OUTPUT_FN)
"""
function plot_sim(
    floe_fn,
    initial_state_fn,
    Δt,
    output_fn
)
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

    # Set up figure
    fig = Figure()
    Axis(
        fig[1, 1];
        limits = (xmin, xmax, ymin, ymax),
        title = title_string,
        xlabel = "Meters",
        xticklabelrotation = pi/4,
        ylabel = "Meters",
        yticklabelrotation = pi/4,
    )

    # Plot starting state (floes + topography)
    poly!(time_polys, color = :lightblue, strokecolor = :black, strokewidth = 0.5)  # floes
    if !isempty(domain.topography)  # topography
        topo_mp = domain.topography.poly
        poly!(topo_mp, color = :lightgrey)
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
    record(fig, output_fn, timesteps; framerate = 20) do t
        time[] = t  # updating the time auto updates the related time_polys
    end

    # End movie and close file
    close(file)
end

end  # module