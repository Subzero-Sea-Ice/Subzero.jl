# # Restart Simulation

# This simulation is split into two pieces, run one after the other. The purpose of this example
# is to show how to use the restart functionality. 

using Subzero, CairoMakie, GeoInterfaceMakie
using JLD2, Random, Statistics

const FT = Float64
const Δt = 20
const nΔt = 2500
const nfloes = 20
const L = 1e5
const Δgrid = 1e4
const hmean = 2
const concentration = 0.7
const uomax = 2

dirs = [joinpath("restart_sim", "run" * string(i)) for i in 1:2]

# ## Create Grid
grid = RegRectilinearGrid(; x0 = 0.0, xf = L, y0 = 0.0, yf = L, Δx = Δgrid, Δy = Δgrid)

# ## Create Domain
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary,west =  wboundary)

# ## Create Ocean
ngrid = Int(L/Δgrid) + 1
ygrid = range(0,L,ngrid)
uoprofile = @. uomax * (1 - abs(1 - 2 * ygrid/L))
uvels_ocean = repeat(
    uoprofile,
    outer = (1, ngrid),
)
ocean = Ocean(;
    u = uvels_ocean',
    grid,
    v = 0,
    temp = 0,
)

# ## Create Atmos
atmos = Atmos(FT; grid, u = 0.0, v = 0.0, temp = 0.0)

# ## Create Floes
floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(grid, 2))
floe_arr = initialize_floe_field(
    FT,
    nfloes,
    [concentration],
    domain,
    hmean,
    0;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

# ## Create Model
model = Model(grid, ocean, atmos, domain, floe_arr)

# ## Create Outout Writers
initwriter = InitialStateOutputWriter(dir = dirs[1], overwrite = true)
checkpointer = CheckpointOutputWriter(
    250,
    dir = dirs[1],
    filename = "checkpoint.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)
floewriter = FloeOutputWriter(50, dir = dirs[1], overwrite = true)
writers = OutputWriters(initwriter, floewriter, checkpointer)

# ## Create Simulation and Constants 
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, f = 0, turnθ = 0)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = nΔt,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
)

# ## Run the first part of the simulation
run!(simulation)

# ## Run second part of the simulation
new_initwriter = InitialStateOutputWriter(initwriter; dir = dirs[2])
new_floewriter = FloeOutputWriter(floewriter; dir = dirs[2])
writers = OutputWriters(new_initwriter, new_floewriter)  # didn't save checkpoint since not restarting
Subzero.restart!(
    dirs[1] * "/initial_state.jld2",
    dirs[1] * "/checkpoint.jld2",
    nΔt, writers; start_tstep = nΔt)


# ## Plot all simulation parts
for i in 1:2
    plot_sim(
        joinpath(dirs[i], "floes.jld2"),
        joinpath(dirs[i], "initial_state.jld2"),
        Δt,
        joinpath(dirs[i], "restart_floes.mp4"),
    )
end

# ### First Part of Simulation 
# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../restart_sim/run1/restart_floes.mp4" type="video/mp4">
# </video>
# ```

# ### Second Part of Simulation 
# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../restart_sim/run2/restart_floes.mp4" type="video/mp4">
# </video>
# ```