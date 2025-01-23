using Subzero, JLD2, CairoMakie, GeoInterfaceMakie, Random, Statistics

const FT = Float64
const Lx = 8e4
const Ly = 8e4
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20
const nΔt = 15000

# Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

uvels = repeat(
    [range(0, 0.2, length = 21); range(0.2, 0, length = 20)],
    outer = (1, 41),
)
ocean = Ocean(;
    u = uvels',
    grid,
    v = 0.0,
    temp = 0.0,
)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe creation
floe_arr = initialize_floe_field(FT, 400, [0.85], domain, hmean, Δh; rng = Xoshiro(1))

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
dir = "output/title/"

# Output setup
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(floewriter)

simulation = Simulation(; model, consts, Δt, nΔt, writers, verbose = true, rng = Xoshiro(1))
run!(simulation)

function plot_logo(floe_fn, Lx, Ly, dir)
    # Floe Information
    file = jldopen(floe_fn)
    sim_polys = file["poly"]
    timesteps = keys(sim_polys)

    # Set up observables for recording (updated whenever `time[]` is set to a new value)
    time = Observable(timesteps[1])  # note these are strings as we index into a JLD2 file
    time_polys = @lift(sim_polys[$time])

    # Set up figure
    fig = Figure(; size = (800, 200))

    floe_rgb = (255, 255, 255)
    floe_color = RGBf((floe_rgb ./ 255)...)
    ocean_rgb = (0, 157, 196)
    ocean_color = RGBf((ocean_rgb ./ 255)...)

    ax = Axis(fig[1, 1]; limits = (0.0, Lx, 0.0, Ly), backgroundcolor = ocean_color)
    hidedecorations!(ax)

    # Plot starting state (floes + topography)
    poly!(time_polys, color = floe_color, strokecolor = :black, strokewidth = 0.5)  # floes
    text!(
        Lx/2, Ly/2;
        align = (:center, :center),
        text = "Subzero.jl",
        font = :bold,
        fontsize = 130,
        strokecolor = :white,
        strokewidth = 5,
    )

    # Create movie
    record(fig, dir*"title.gif", timesteps; framerate = 20) do t
        time[] = t  # updating the time auto updates the related time_polys
    end

    # End movie and close file
    close(file)
end

plot_logo(dir*"floes.jld2", Lx, Ly, dir)