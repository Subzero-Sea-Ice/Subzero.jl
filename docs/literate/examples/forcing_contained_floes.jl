# # Floes Bounded by Ocean Currents

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../forcing_contained_floes/contained_floes.mp4" type="video/mp4">
# </video>
# ```

# This simulation has four `open` boundaries. The ocean is created such that there is a current in the 
# middle of the domain that pushes the floes from left to right, and there are also currents at each of
# boundaries that push the floes back into the middle of the domain. The main point of this simulation
# is to highlight that the user can create a bounding box to control initial floe placement and that
# the initial set of floes does not need to span the entire domain.

using Subzero, CairoMakie, GeoInterfaceMakie
using JLD2, Random, Statistics

#  ## User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20
const nΔt = 10000;

# ## Grid Creation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

# ## Domain Creation
nboundary = OpenBoundary(North; grid)
sboundary = OpenBoundary(South; grid)
eboundary = OpenBoundary(East; grid)
wboundary = OpenBoundary(West; grid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# ## Ocean Creation

# ### Set ocean u and v velocities
ocean_uvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))
ocean_vvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))

for i in CartesianIndices(ocean_vvels)
    r, c = Tuple(i)
    if r ≤ 5  # set u vals
        ocean_uvels[i] = 0.6
    elseif r ≥ grid.Nx - 4
        ocean_uvels[i] = -0.6
    end
    if c ≤ 5  # set v vals
        ocean_vvels[i] = 0.6
    elseif c ≥ grid.Ny - 4
        ocean_vvels[i] = -0.6
    end
end
ocean_uvels[10:40, 20:30] .= 0.3;

# ### Instatiate Ocean
ocean = Ocean(;
    grid, 
    u = ocean_uvels,
    v = ocean_vvels,
    temp = 0,
)

# We can then plot the ocean for a better understanding of the above setup.
fig = Figure();
ax1 = Axis(fig[1, 1]; title = "Ocean U-Velocities [m/s]", xticklabelrotation = pi/4)
ax2 = Axis(fig[2, 1]; title = "Ocean V-Velocities [m/s]", xticklabelrotation = pi/4)
xs = grid.x0:grid.Δx:grid.xf
ys = grid.y0:grid.Δy:grid.yf
u_hm = heatmap!(ax1, xs, ys, ocean.u)
Colorbar(fig[1, end+1], u_hm)
v_hm = heatmap!(ax2, xs, ys, ocean.v)
Colorbar(fig[2, end+1], v_hm)
resize_to_layout!(fig)
fig

# ## Atmosphere Creation
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)
# Floe Creation - bound floes within smaller part of the domain
floe_bounds = Subzero.make_polygon([[[0.1Lx, 0.1Ly], [0.1Lx, 0.9Ly], [0.9Lx, 0.9Ly], [0.9Lx, 0.1Ly], [0.1Lx, 0.1Ly]]])
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

# ## Model Creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# ## Output Writer Creation
dir = "forcing_contained_floes"
init_fn, floe_fn = "contained_floes_init_state.jld2", "contained_floes.jld2"
initwriter = InitialStateOutputWriter(filename = init_fn, dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, filename = floe_fn, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# ## Simulation Creation
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = nΔt,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
)

# ## Running the Simulation
run!(simulation)

# ## Plotting the Simulation
plot_sim(joinpath(dir, floe_fn), joinpath(dir, init_fn), Δt, joinpath(dir, "contained_floes.mp4"))

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../forcing_contained_floes/contained_floes.mp4" type="video/mp4">
# </video>
# ```

# !!! note
#       Note that this is just using the built-in basic plotting. However, it is easy to write
#       your own plotting code. See the source code for a basic outline.
