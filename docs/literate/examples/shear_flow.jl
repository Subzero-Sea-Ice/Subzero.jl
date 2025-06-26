
# # Shear Flow with Periodic Boundaries

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../shear_flow/shear_flow.mp4" type="video/mp4">
# </video>
# ```

# This simulation creates a periodic simulation with an ocean shear flow current. It is very
# similar to the example created in the tutorial. The tutorial explains the
# setup in greater detail. It is a good starter simulation for understanding how the code
# works and playing around with a basic setup.

# This simulation shows how to setup a basic `domain` where all four boundaries are `periodic`.
# It also creates non-constant, shear `ocean` field. The `ocean` `u`-velocities are zero at
# the minimum and maximum y-extents (constant across x-values). The `u` values increase from
# the top and bottom of the domain towards the center to 0.5m/s, again constant across x-values,
# varying with y-values.

using Subzero, CairoMakie, GeoInterfaceMakie
using JLD2, Random, Statistics

# ## User Inputs
const FT = Float64  # Float type used to run simulation
const Lx = 1e5      # grid x-length
const Ly = 1e5      # grid y-length
const Δgrid = 2e3   # grid cell edge-size
const hmean = 0.25  # mean floe height
const Δh = 0.0      # difference in floe heights - here all floes are the same height
const Δt = 20       # timestep
const nΔt = 5000;   # number of timesteps to run

# ## Grid Creation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

# ## Domain Creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# ## Ocean Creation
half_nx = fld(grid.Nx, 2)
u_vec = [range(0, 0.5, length = half_nx); 0.5; range(0.5, 0, length = half_nx)]
u_row = transpose(u_vec)
uvels = repeat(u_row, outer = (grid.Ny + 1, 1))
ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)

# ## Atmosphere Creation
atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)

# ## Floe Creation
floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(grid, 2))
floe_arr = initialize_floe_field(
    FT,
    50,
    [0.75],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

# ## Model Creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# ## Constants Creation
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# ## Output Creation
dir = "shear_flow"
init_fn, floe_fn = "shear_flow_init_state.jld2", "shear_flow_floes.jld2"
initwriter = InitialStateOutputWriter(dir = dir, filename = init_fn, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, filename = floe_fn, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# ## Simulation Creation
simulation = Simulation(; model, consts, writers, Δt, nΔt, floe_settings,
    verbose = true, rng = Xoshiro(1))

# ## Running the Simulation
run!(simulation)

# ## Plotting the Simulation
plot_sim(joinpath(dir, floe_fn), joinpath(dir, init_fn), Δt, joinpath(dir, "shear_flow.mp4"))

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../shear_flow/shear_flow.mp4" type="video/mp4">
# </video>
# ```

# !!! note
#       Note that this is just using the built-in basic plotting. However, it is easy to write
#       your own plotting code. See the source code for a basic outline.
