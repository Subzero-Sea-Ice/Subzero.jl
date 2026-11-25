using JLD2, Random, Statistics, Subzero, CairoMakie, GeoInterfaceMakie

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 1e4
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

uvels = repeat(
    hcat(
        range(0.1 ,0.6, step = 0.1)',
        range(0.5, 0.1, step = -0.1)',
    ),
    outer = 11,
)
ocean = Ocean(; u = uvels, grid, v = 0, temp = 0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    50,
    [0.8],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Output writers
dir = "output/converge_diverge"
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
checkpointwriter = CheckpointOutputWriter(1000, dir = dir, overwrite = true)

writers = OutputWriters(initwriter, floewriter, checkpointwriter)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 3000,
    verbose = true,
    writers = writers,
)

# Run simulation
run!(simulation)

plot_sim(
    "output/converge_diverge/floes.jld2",
    "output/converge_diverge/initial_state.jld2",
    20,
    "output/converge_diverge/converge_diverge.mp4",
)