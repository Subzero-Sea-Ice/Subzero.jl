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

# Set ocean u velocities
ocean_uvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))
for i in CartesianIndices(ocean_uvels)
    r, c = Tuple(i)
    if r ≤ 5
        ocean_uvels[i] = 0.2
    elseif r ≥ grid.Nx - 4
        ocean_uvels[i] = -0.2
    end
end
ocean_uvels[20:40, 20:30] .= 0.15
# Set ocean v velocities
ocean_vvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))
for i in CartesianIndices(ocean_vvels)
    r, c = Tuple(i)
    if c ≤ 5
        ocean_vvels[i] = 0.2
    elseif c ≥ grid.Ny - 4
        ocean_vvels[i] = -0.2
    end
end
# Create ocean
ocean = Ocean(
    u = ocean_uvels,
    v = ocean_vvels,
    grid, 
    temp = 0,
)

# Create atmosphere
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = OpenBoundary(North; grid)
sboundary = OpenBoundary(South; grid)
eboundary = OpenBoundary(East; grid)
wboundary = OpenBoundary(West; grid)

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)
# Floe creation - bound floes within smaller part of the domain
floe_bounds = [[[0.1Lx, 0.1Ly], [0.1Lx, 0.9Ly], [0.9Lx, 0.9Ly], [0.9Lx, 0.1Ly], [0.1Lx, 0.1Ly]]]
floe_arr = initialize_floe_field(
    FT,
    300,
    [0.4],
    domain,
    hmean,
    Δh;
    floe_bounds = floe_bounds,
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
dir = "output/contained"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 15000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
)
run_time!(simulation)

plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    Δt,
    joinpath(dir, "contained.mp4"),
)