using JLD2, Random, Statistics, Subzero, CairoMakie, GeoInterfaceMakie

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.125
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = MovingBoundary(North; grid, u = 0.0, v = -0.1)
sboundary = MovingBoundary(South; grid, u = 0.0, v = 0.1)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    500,
    [1.0],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
)
nfloes = length(floe_arr)
floe_arr.u .= 0
floe_arr.v .= -0.01
# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, Cd_io = 0.0, f = 0.0, turnθ = 0.0)

# Run simulation
run_time!(simulation) = @time run!(simulation)
dir = "output/moving_bounds"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(40, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# Simulation settings 
ridgeraft_settings = RidgeRaftSettings(
    ridge_raft_on = true,
    Δt = 150,
    domain_gain_probability = 0.5
)
weld_settings = WeldSettings(
    weld_on = true,
    Δts = [150, 300, 600],
    Nxs = [2, 1, 1],
    Nys = [2, 2, 1],
)

coupling_settings = CouplingSettings(two_way_coupling_on = true)


simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 5000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
    ridgeraft_settings = ridgeraft_settings,
)
run_time!(simulation)

Subzero.plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    20,
    joinpath(dir, "moving_bounds.mp4"),
)