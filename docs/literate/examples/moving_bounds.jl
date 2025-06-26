# # Moving Boundaries

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../moving_bounds/moving_bounds.mp4" type="video/mp4">
# </video>
# ```

# This simulation has two moving boundaries, on the top and the bottom of the simulation. They push
# the floes inward towards the center of the domain. As the floes are pushed inward, they ridge and raft,
# gaining height and losing area. Users can also create shear moving boundaries, rather than the compression
# boundaries seen in this example.

using Subzero, CairoMakie, GeoInterfaceMakie
using JLD2, Random, Statistics

# ## User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.125
const Δt = 20
const nΔt = 1500;

# ## Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# ## Domain creation
nboundary = MovingBoundary(North; grid, u = 0.0, v = -0.1)
sboundary = MovingBoundary(South; grid, u = 0.0, v = 0.1)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# ## Floe creation
floe_arr = initialize_floe_field(
    FT,
    100,
    [1.0],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
)
nfloes = length(floe_arr)
floe_arr.u .= 0  # set the inital floe velocities manually
floe_arr.v .= -0.01;

# ## Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# ## Output Writer Setup
dir = "moving_bounds"
init_fn, floe_fn = "moving_bounds_init_state.jld2", "moving_bounds.jld2"
initwriter = InitialStateOutputWriter(dir = dir, filename = init_fn, overwrite = true)
floewriter = FloeOutputWriter(50, dir= dir, filename = floe_fn, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# ## Simulation settings 
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, Cd_io = 0.0, f = 0.0, turnθ = 0.0)

ridgeraft_settings = RidgeRaftSettings(
    ridge_raft_on = true,
    Δt = 150,
    domain_gain_probability = 0.5
)
weld_settings = WeldSettings(
    weld_on = true,
    Δts = [150, 300, 600],  # weld at these specific timesteps
    Nxs = [2, 1, 1],  # split the domain into nx by ny sections and weld within each section
    Nys = [2, 2, 1],
)

# ## Create Simulation
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 5000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    ridgeraft_settings = ridgeraft_settings,
)

# ## Running the Simulation
run!(simulation)

# ## Plotting the Simulation
plot_sim(joinpath(dir, floe_fn), joinpath(dir, init_fn), Δt, joinpath(dir, "moving_bounds.mp4"))

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../moving_bounds/moving_bounds.mp4" type="video/mp4">
# </video>
# ```

# !!! note
#       Note that this is just using the built-in basic plotting. However, it is easy to write
#       your own plotting code. See the source code for a basic outline.