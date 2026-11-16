using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays, SplitApplyCombine, Test
import LibGEOS as LG

function setup_floes_with_inters(coords, domain, consts,
    collision_settings, lock,  Δt, Δx = nothing, Δy = nothing,
)
    floes = initialize_floe_field(
        Float64,
        coords,
        domain,
        1.0,
        0.0;
    )
    if !isnothing(Δx)
        for i in eachindex(Δx)
            Subzero.translate!(floes.coords[i], Δx[i], Δy[i])
            floes.centroid[i][1] += Δx[i]
            floes.centroid[i][2] += Δy[i]
        end
    end
    add_ghosts!(floes, domain)
    Subzero.timestep_collisions!(  # Add interactions
        floes,
        length(floes),
        domain,
        consts,
        10,
        collision_settings, 
        lock,
    )
    return floes
end
function update_height(floes, i, new_height, consts)
    floes.height[i] = new_height
    floes.mass[i] = floes.area[i] * floes.height[i] * consts.ρi
    floes.moment[i] = Subzero._calc_moment_inertia(
        floes.coords[i],
        floes.centroid[i],
        floes.height[i],
        ρi = consts.ρi,
    )
end

grid = RegRectilinearGrid(
    (0, 1e5),
    (0, 1e5),
    1e4,
    1e4,
)
topo_coords = [[[5e4, 5e4], [5e4, 7e4], [7e4, 7e4], [7e4, 5e4], [5e4, 5e4]]]
collision_domain = Subzero.Domain(;
    north = CollisionBoundary(North; grid),
    south = CollisionBoundary(South; grid),
    east = CollisionBoundary(East; grid),
    west = CollisionBoundary(West; grid),
    topography = initialize_topography_field(; coords = [topo_coords])
)

periodic_domain = Subzero.Domain(; 
    north = PeriodicBoundary(North; grid),
    south = PeriodicBoundary(South; grid),
    east = PeriodicBoundary(East; grid),
    west = PeriodicBoundary(West; grid),
)

consts = Constants()
coupling_settings = CouplingSettings()
simp_settings = SimplificationSettings(min_floe_area = 1e7)
collision_settings = CollisionSettings(floe_floe_max_overlap = 1.0) # don't fuse
lock = Threads.SpinLock()

ridge_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 1.0,  # no ridging
    raft_probability = 0.0,  # no rafting
    min_overlap_frac = 0.0001,
)
coords = [
    [[[3e4, -0.2e4], [3e4, 0.2e4], [5e4, -0.1e4], [8e4, 0.2e4], [8e4, -0.2e4], [3e4, -0.2e4]]]
]
base_floes = setup_floes_with_inters(coords, collision_domain, consts,
    collision_settings, lock, Δt,
)
no_rr_frac_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 1.0,
    raft_probability = 1.0,
    min_overlap_frac = 1.0,  # need 100% overlap to ridge/raft
)
ridge_keep_mass_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 1.0,  # force ridging
    raft_probability = 0.0,
    min_overlap_frac = 0.001,
    domain_gain_probability = 0.0,
)

floes = deepcopy(base_floes)
#update_height(floes, 1, 1.0, consts)  # floe 1 will ridge onto floe 2
# update_height(floes, 1, 0.1, consts)
# update_height(floes, 3, 0.1, consts)
floes.status[1].tag = Subzero.active
# total_mass = floes.mass[1] + floes.mass[2]
# h1, h2 = floes.height
# area1, area2 = floes.area
# cent1, cent2 = floes.centroid
pieces_list = StructArray{Floe{Float64}}(undef, 0)

max_id = Subzero.timestep_ridging_rafting!(
    floes,
    pieces_list,
    collision_domain,
    maximum(floes.id),
    coupling_settings,
    ridge_keep_mass_settings,
    simp_settings,
    consts,
    10,
)

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 10
# Model instantiation
grid = RegRectilinearGrid(
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
zonal_ocn = Ocean(; grid, u = 0.5, v = 0.0, temp = 0.0)

zero_atmos = Atmos(; grid, u = 0.5, v = 0.0, temp = 0.0)


domain = Subzero.Domain(;
    north = PeriodicBoundary(North; grid),
    south = PeriodicBoundary(South; grid),
    east = PeriodicBoundary(East; grid),
    west = PeriodicBoundary(West; grid),
)

# Floe instantiation
coords = [[
    [0.0, 0.0],
    [5e4, 1e5],
    [1e5, 1e5],
    [1e5, 0.0],
    [0.0, 0.0],
]]

floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 1),
)
coupling_settings = CouplingSettings(
    two_way_coupling_on = true,
)

floe_arr = initialize_floe_field(
    FT,
    30,
    [0.8],
    domain,
    0.5,
    0.0;
    floe_settings = floe_settings,
    rng = Xoshiro(1),
)
floe_arr.α .= 1e-5
model = Model(
    grid,
    zonal_ocn,
    zero_atmos,
    domain,
    floe_arr,
)
dir = "output/sim"
writers = OutputWriters(
    # InitialStateOutputWriter(
    #     dir = dir,
    #     overwrite = true
    # ),
    # FloeOutputWriter(
    #     250,
    #     dir = dir,
    #     overwrite = true,
    # ),
)
simulation = Simulation(
    name = "sim",
    model = model,
    Δt = 10,
    nΔt = 200,
    writers = writers,
    verbose = true,
    consts = Constants(Cd_ia = 0, Cd_ao = 0),
    coupling_settings = coupling_settings,
    floe_settings = floe_settings,
    
)

#@benchmark timestep_sim!(simulation, 10) setup=(sim=deepcopy(simulation))
function run_with_saving(simulation, tstep_between_save)
    nsave = fld(simulation.nΔt, tstep_between_save) + 1
    τx = zeros(FT, grid.Nx + 1, grid.Ny + 1, nsave)
    Subzero.startup_sim(simulation, nothing, 1)
    tstep = 0
    counter = 1
    while tstep <= simulation.nΔt
        # Timestep the simulation forward
        timestep_sim!(simulation, tstep)
        if mod(tstep, tstep_between_save) == 0
            τx[:, :, counter] .= simulation.model.ocean.τx
            counter += 1
        end
        tstep+=1
    end
    Subzero.teardown_sim(simulation)
    return τx
end

τx = run_with_saving(simulation, 50)


# time_run(simulation) = @time run!(simulation)
# time_run(simulation)

# Subzero.create_sim_gif("output/sim/floes.jld2", 
#                        "output/sim/initial_state.jld2",
#                        "output/sim/sim.gif")
# # Run simulation
#time_run(simulation)
#Profile.Allocs.clear()
#@time run!(simulation)
#ProfileView.@profview run!(simulation)
#Profile.Allocs.@profile timestep_sim!(simulation, 1)
# Profile.Allocs.@profile sample_rate=1 Subzero.timestep_coupling!(
#     simulation.model.floes,
#     simulation.model.grid,
#     simulation.model.domain,
#     simulation.model.ocean,
#     simulation.model.atmos,
#     simulation.consts,
#     simulation.coupling_settings,
# )


# Profile.Allocs.@profile sample_rate=1 Subzero.timestep_coupling!(
#     simulation.model,
#     simulation.consts,
#     simulation.coupling_settings,
#     Threads.SpinLock(),
# )

# PProf.Allocs.pprof(from_c = false)
# last(sort(results.allocs, by=x->x.size))
# Subzero.create_sim_gif(
#     joinpath(dir, "floes.jld2"), 
#     joinpath(dir, "initial_state.jld2"),
#     joinpath(dir, "test.gif"),
# )
