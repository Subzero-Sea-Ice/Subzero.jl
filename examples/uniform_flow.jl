using JLD2, Random, Statistics, Subzero, CairoMakie, GeoInterfaceMakie

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
ocean = Ocean(; grid, u = 0.1, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)

domain = (; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    5,
    [0.4],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) = @time run!(simulation)
dir = "output/uniform_flow"

    # Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 20,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = CouplingSettings(two_way_coupling_on = true),
)
run_time!(simulation)

Subzero.plot_sim(
    "output/uniform_flow/floes.jld2",
    "output/uniform_flow/initial_state.jld2",
    20,
    "output/uniform_flow/uniform_flow.mp4",
)