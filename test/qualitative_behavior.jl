"""
Simulations that show simple, fundamental behavior of sea ice.
Run and check these simulations to confirm behavior has not changed. 
"""

using Subzero, StructArrays, JLD2

const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const nΔt = 4000
const Δt = 10
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Setup for Simulations
grid = RegRectilinearGrid(; x0 = -2.5e4, xf = Lx, y0 = -2.5e4, yf = Ly, Δx = Δgrid, Δy = Δgrid)

zero_ocn = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
meridional_ocn = Ocean(; grid, u = 0.0, v = 1.0, temp = 0.0)

zero_atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)
zonal_atmos = Atmos(; grid, u = -15.0, v = 0.0, temp = 0.0)

open_domain_no_topo = Subzero.Domain(;
    north = OpenBoundary(North; grid),
    south = OpenBoundary(South; grid),
    east = OpenBoundary(East; grid),
    west = OpenBoundary(West; grid),
)

collision_domain_topo = Subzero.Domain(;
    north = CollisionBoundary(North; grid),
    south = CollisionBoundary(South; grid),
    east = CollisionBoundary(East; grid),
    west = CollisionBoundary(West; grid),
    topography = initialize_topography_field(; coords = [[[[2e4, 0.0], [2e4, 2e4], [2.5e4, 2e4], [2.5e4, 0.0], [2e4, 0.0]]]])
)

stationary_rect_floe = StructArray([Floe(
    [[
        [0.0, 0.0],
        [0.0, 2e4],
        [0.5e4, 2e4],
        [0.5e4, 0.0],
        [0.0, 0.0],
    ]],
    hmean,
    Δh,
)])
zonal_3rect_floes = initialize_floe_field(
    FT,
    [  # List of 3 floe coordinates
        [[
            [0.0, 0.0],
            [0.0, 2e4],
            [-0.5e4, 2e4],
            [-0.5e4, 0.0],
            [0.0, 0.0]
        ]],
        [[
            [5e4, 0.0],
            [5e4, 2e4],
            [5.5e4, 2e4],
            [5.5e4, 0.0],
            [5e4, 0.0],
        ]],
        [[
            [8.e4, 0.0],
            [8e4, 2e4],
            [8.5e4, 2e4],
            [8.5e4, 0.0],
            [8e4, 0.0],
        ]],
    ],
    collision_domain_topo,
    hmean,
    Δh;
)
zonal_3rect_floes.u .= [3.0, -3.0, 0.0]

collisions_off_settings = CollisionSettings(collisions_on = false)

sim_arr = Vector{Simulation}(undef, 0)
"""
Simulation 1:
    One floe pushed by meridional (south-to-north) 1m/s ocean flow. Floe is
    initally stationary.
Expected Behavior:
    Floe velocity quickly reaches ocean velocity and flows northward. 
"""
model1 = Model(
    grid,
    meridional_ocn,
    zero_atmos,
    open_domain_no_topo,
    deepcopy(stationary_rect_floe),
)
writers1 = OutputWriters(
    InitialStateOutputWriter(
        dir = "test/output/sim1",
        overwrite = true
    ),
    FloeOutputWriter(
        30,
        dir = "test/output/sim1",
        overwrite = true,
    ),
)
simulation1 = Simulation(
    name = "sim1",
    model = model1,
    Δt = Δt,
    nΔt = nΔt,
    collision_settings = collisions_off_settings,
    writers = writers1,
)
push!(sim_arr, simulation1)
"""
Simulation 2:
    One floe pushed by zonal (west-to-east) -15m/s atmos flow. Floe is initally
    stationary.
Expected Behavior:
    Floe should drift to the right of movement due to Coriolis force in the
    northern hemisphere. 
"""
model2 = Model(
    grid,
    zero_ocn, 
    zonal_atmos,
    open_domain_no_topo,
    deepcopy(stationary_rect_floe),
)
writers2 = OutputWriters(
    InitialStateOutputWriter(
        dir = "test/output/sim2",
        overwrite = true
    ),
    FloeOutputWriter(
        30,
        dir = "test/output/sim2",
        overwrite = true,
    ),
)
simulation2 = Simulation(
    name = "sim2",
    model = model2,
    Δt = Δt,
    nΔt = nΔt,
    collision_settings = collisions_off_settings,
    writers = writers2,
)
push!(sim_arr, simulation2)

"""
Simulation 3:
    One floe with initial velocity flows into a collision boundary and then
    bounces off a topography element. Another floe woth opposite initial
    velocity hits the other side of the collision boundary. No drag so the floes
    can move free of the ocean.
Expected Behavior:
    Floes should bounce off of the topography elements, walls, and each other. 
"""
model3 = Model(
    grid,
    zero_ocn,
    zero_atmos,
    collision_domain_topo,
    deepcopy(zonal_3rect_floes),
)
writers3 = OutputWriters(
    InitialStateOutputWriter(
        dir = "test/output/sim3",
        overwrite = true
    ),
    FloeOutputWriter(
        30,
        dir = "test/output/sim3",
        overwrite = true,
    ),
)
simulation3 = Simulation(
    name = "sim3",
    model = model3,
    Δt = Δt,
    nΔt = nΔt,
    writers = writers3,
    coupling_settings = CouplingSettings(coupling_on = false),
)
push!(sim_arr, simulation3)

"""
Simulation 4:
    Two initial floes and double periodic boundaries. One floe has (1, 1)
    initial velocity and the other has (1, 0) velocity. One topography element
    accross from second floe. No drag so the floes can move free of ocean.
Expected Behavior:
    Floe one should pass through top right the corner and you should see 3 ghost
    floes appear before it eventually passes through bottom the left corner and
    later bounces off of the topogrpahy element. Floe 2 passes eastern wall,
    populating a ghost floe which  hits the topography element before bounding
    back through the western wall. 
"""
periodic_bounds_topo = Subzero.Domain(;
    north = PeriodicBoundary(North; grid),
    south = PeriodicBoundary(South; grid),
    east = PeriodicBoundary(East; grid),
    west = PeriodicBoundary(West; grid),
    topography = StructArray([TopographyElement(; coords = [[[-1.5e4, 4.5e4], [-1.5e4, 6.5e4], [2.5e4, 6.5e4], [2.5e4, 4.5e4], [-1.5e4, 4.5e4]]])]),
)

p1_coords = [[
    [7.5e4, 7.5e4],
    [7.5e4, 9.5e4],
    [9.5e4, 9.5e4],
    [9.5e4, 7.5e4],
    [7.5e4, 7.5e4],
]]
p2_coords = [[
    [6.5e4, 4.5e4],
    [6.5e4, 6.5e4],
    [8.5e4, 6.5e4],
    [8.5e4, 4.5e4],
    [6.5e4, 4.5e4],
]]
p_floe_arr = StructArray(
    [Floe(c, hmean, Δh) for c in [p1_coords, p2_coords]]
)
p_floe_arr.u[1] = 1
p_floe_arr.v[1] = 1
p_floe_arr.u[2] = 1
model4 = Model(
    grid,
    zero_ocn,
    zero_atmos,
    periodic_bounds_topo,
    deepcopy(p_floe_arr),
)
writers4 = OutputWriters(
    InitialStateOutputWriter(
        dir = "test/output/sim4",
        overwrite = true
    ),
    FloeOutputWriter(
        30,
        dir = "test/output/sim4",
        overwrite = true,
    ),
)
simulation4 = Simulation(
    name = "sim4",
    model = model4,
    Δt = Δt,
    nΔt = nΔt,
    writers = writers4,
    coupling_settings = CouplingSettings(coupling_on = false),
)
push!(sim_arr, simulation4)


"""
Simulation 5:
    Input file full of strangley shaped floes, each given a very small initial
    velocity. Collisions are enabled.
Expected Behavior:
    Floes should all bounce off of one another, without becoming unstable.
"""
file = jldopen("test/inputs/floe_shapes.jld2", "r")
funky_floe_coords = file["floe_vertices"][1:100]
funky_floe_arr = initialize_floe_field(
    FT,
    funky_floe_coords,
    collision_domain_topo,
    hmean,
    Δh;
)
close(file)
funky_floe_arr.u .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))
funky_floe_arr.v .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))

model5 = Model(
    grid,
    zero_ocn,
    zero_atmos,
    open_domain_no_topo,
    deepcopy(funky_floe_arr),
)
writers5 = OutputWriters(
    InitialStateOutputWriter(
        dir = "test/output/sim5",
        overwrite = true
    ),
    FloeOutputWriter(
        30,
        dir = "test/output/sim5",
        overwrite = true,
    ),
)
simulation5 = Simulation(
    name = "sim5",
    model = model5,
    Δt = Δt,
    nΔt = nΔt,
    writers = writers5,
    coupling_settings = CouplingSettings(coupling_on = false),
)
push!(sim_arr, simulation5)



# Run the simulations
for sim in sim_arr
    run!(sim)
    plot_sim(
        joinpath("test/output", sim.name, "floes.jld2"),
        joinpath("test/output", sim.name, "initial_state.jld2"),
        Δt,
        joinpath("test/output", sim.name, string(sim.name, ".mp4")),
    )
end