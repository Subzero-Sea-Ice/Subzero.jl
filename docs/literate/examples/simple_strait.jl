
# # Simple Strait Simulation

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../simple_strait/simple_strait.mp4" type="video/mp4">
# </video>
# ```

# This simulation creates a north to south strait that ice can flow through, pushed by the
# ocean. The north and south boundaries form a `periodic` pair, so that the ice can endlessly flow
# through the strait. The east and west boundaries are `collision` bounds, but they are
# completely covered with topography forming the edges of the domain. This is a good simulation
# to understand how to setup topography and how to turn on fractures using the fracture settings.

using Subzero, CairoMakie, GeoInterfaceMakie
using JLD2, Random, Statistics

# ## Making the Simulation

# ## User Inputs
const FT = Float64  # Float type used to run simulation
const Lx = 1e5      # grid x-length
const Ly = 1e5      # grid y-length
const Δgrid = 2e3   # grid cell edge-size
const hmean = 0.25  # mean floe height
const Δh = 0.0      # difference in floe heights - here all floes are the same height
const Δt = 20       # timestep
const nΔt = 5000;    # number of timesteps to run

# ## Grid Creation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

# ## Domain Creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = CollisionBoundary(East; grid)
wboundary = CollisionBoundary(West; grid)

island1 = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
island2 = [[[4e4, 6e4], [4e4, 6.5e4], [4.5e4, 6.5e4], [4.5e4, 6e4], [4e4, 6e4]]]
topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]
topo_arr = initialize_topography_field(FT; coords = [island1, topo1, topo2])

domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary, topography = topo_arr)

# ## Ocean Creation
ocean = Ocean(; grid, u = 0.0, v = -0.3, temp = 0.0)

# ## Atmos Creation
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)

# ## Floe Creation
floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
    stress_calculator = DecayAreaScaledCalculator(),
)

floe_arr = initialize_floe_field(
    FT,
    75,
    [0.7],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(3),
    floe_settings = floe_settings,
)

# ## Model Creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# ## Constants Creation
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# ## Settings Creation
# ### Fracture Settings
fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        deform_on = false,
)
# ### Ridge Raft Settings
ridgeraft_settings = RidgeRaftSettings(
    ridge_raft_on = true,
    Δt = 150
)

# ## Output Creation
dir = "simple_strait"
init_fn, floe_fn = "simple_strait_init_state.jld2", "simple_strait_floes.jld2"
initwriter = InitialStateOutputWriter(dir = dir, filename = init_fn, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, filename = floe_fn, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# ## Simulation Creation

simulation = Simulation(; model, consts, writers, Δt, nΔt,
    floe_settings, fracture_settings, ridgeraft_settings,
    verbose = true, rng = Xoshiro(1))
    
# ## Running the Simulation
run!(simulation)

# ## Plotting the Simulation
plot_sim(joinpath(dir, floe_fn), joinpath(dir, init_fn), Δt, joinpath(dir, "simple_strait.mp4"))

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../simple_strait/simple_strait.mp4" type="video/mp4">
# </video>
# ```

# !!! note
#       Note that this is just using the built-in basic plotting. However, it is easy to write
#       your own plotting code. See the source code for a basic outline.