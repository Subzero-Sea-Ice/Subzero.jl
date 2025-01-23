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

uvels = repeat(
    [range(0, 0.5, length = 26); range(0.5, 0, length = 25)],
    outer = (1, 51),
)
ocean = Ocean(;
    u = uvels',
    grid,
    v = 0,
    temp = 0,
)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

coupling_settings = CouplingSettings(
    two_way_coupling_on = true,
)
floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)
# Floe creation
floe_arr = initialize_floe_field(
    FT,
    500,
    [0.8],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) =  @time run!(simulation)
dir = "output/shear_flow"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 5000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
)
run_time!(simulation)
plot_sim(
    "output/shear_flow/floes.jld2",
    "output/shear_flow/initial_state.jld2",
    20,
    "output/shear_flow/shear_flow.mp4",
)