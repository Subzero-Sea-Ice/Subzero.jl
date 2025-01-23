"""
Structs and functions to create and run a Subzero simulation
"""

@kwdef struct Constants{FT<:AbstractFloat}
    ρo::FT = 1027.0             # Ocean density
    ρa::FT = 1.2                # Air density
    Cd_io::FT = 3e-3            # Ice-ocean drag coefficent
    Cd_ia::FT = 1e-3            # Ice-atmosphere drag coefficent
    Cd_ao::FT = 1.25e-3         # Atmosphere-ocean momentum drag coefficient
    f::FT = 1.4e-4              # Ocean coriolis frequency
    turnθ::FT = 15π/180         # Ocean turn angle
    L::FT = 2.93e5              # Latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    ν::FT = 0.3                 # Poisson's ratio
    μ::FT = 0.2                 # Coefficent of friction
    E::FT = 6e6                 # Young's Modulus
end


"""
Constants(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Constants
constructor. A Constants of type FT will be created by passing all other
arguments to the correct constructor. 
"""
Constants(::Type{FT}, args...; kwargs...) where {FT <: AbstractFloat} =
    Constants{FT}(args...; kwargs...)

"""
    Constants(args...)

If a type isn't specified, Constants will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
Constants(args...; kwargs...) = Constants{Float64}(args...; kwargs...)

"""
    Simulation{FT<:AbstractFloat, DT<:Domain{FT}}

Simulation which holds a model and parameters needed for running the simulation.
This includes physical constants (consts), a random number generator (rng), the
number of seconds in a timestep (Δt), and the number of timesteps to run (nΔt).
We also have a flag for verbose, which will print out the number of timesteps
every 50 timesteps and a simulation name, which can be used when saving files.
The user can also define settings for each physical process.
"""
@kwdef struct Simulation{
    FT<:AbstractFloat,
    MT<:Model{FT, <:AbstractRectilinearGrid, <:Domain},
    CT<:AbstractFractureCriteria,
    PT<:AbstractSubFloePointsGenerator{FT},
    ST<:AbstractStressCalculator{FT},
    RT<:Random.AbstractRNG,
    OT<:OutputWriters{
        <:StructVector{<:InitialStateOutputWriter},
        <:StructVector{<:FloeOutputWriter},
        <:StructVector{<:GridOutputWriter},
        <:StructVector{<:CheckpointOutputWriter},
    },
}
    model::MT                           # Model to simulate
    consts::Constants{FT} = Constants() # Constants used in Simulation
    rng::RT = Xoshiro()                     # Random number generator 
    verbose::Bool = false               # String output printed during run
    name::String = "sim"                # Simulation name for printing/saving
    # Timesteps ----------------------------------------------------------------
    Δt::Int = 10                        # Simulation timestep (seconds)
    nΔt::Int = 7500                     # Total timesteps simulation runs for
    # Physical Processes -------------------------------------------------------
    floe_settings::FloeSettings{FT, PT, ST} = FloeSettings()
    coupling_settings::CouplingSettings = CouplingSettings()
    collision_settings::CollisionSettings{FT} = CollisionSettings()
    fracture_settings::FractureSettings{CT} = FractureSettings()
    simp_settings::SimplificationSettings{FT} = SimplificationSettings()
    ridgeraft_settings::RidgeRaftSettings{FT} = RidgeRaftSettings()
    weld_settings::WeldSettings{FT} = WeldSettings()
    # Output Writers -----------------------------------------------------------
    writers::OT = OutputWriters()
end

"""
timestep_sim!(sim, tstep, start_tstep)

Run one step of the simulation and write output. 
Inputs:
    sim          <Simulation> simulation to advance
    tstep        <Int> current timestep
    start_tstep  <Int> timestep simulation started on
Outputs:
    None. Simulation advances by one timestep. 
"""
function timestep_sim!(sim, tstep, start_tstep = 0)
    sim.verbose && mod(tstep, 50) == 0 && println(tstep, " timesteps")
    if !isempty(sim.model.floes)
        max_floe_id = maximum(sim.model.floes.id)
        # Need to lock some operations when multi-threading
        spinlock = Threads.SpinLock()
        # Add ghost floes through periodic boundaries
        n_init_floes = length(sim.model.floes) # number of floes before ghosts
        add_ghosts!(sim.model.floes, sim.model.domain)

        # Output at given timestep
        write_data!(sim, tstep, start_tstep)  # Horribly type unstable
        
        # Collisions
        if sim.collision_settings.collisions_on
            timestep_collisions!(
                sim.model.floes,
                n_init_floes,
                sim.model.domain,
                sim.consts,
                sim.Δt,
                sim.collision_settings,
                spinlock,
            )
        end
        pieces_buffer = StructArray{Floe{Float64}}(undef, 0)
        # Ridge and raft floes that meet overlap conditions
        if (
            sim.ridgeraft_settings.ridge_raft_on &&
            mod(tstep, sim.ridgeraft_settings.Δt) == 0
        )
            max_floe_id = timestep_ridging_rafting!(
                sim.model.floes,
                pieces_buffer,
                sim.model.domain,
                max_floe_id,
                sim.ridgeraft_settings,
                sim.floe_settings,
                sim.simp_settings,
                sim.Δt,
                sim.rng,
            )
        end
        # Remove the ghost floes - only used for collisions
        for i in reverse(n_init_floes+1:length(sim.model.floes))
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    sim.model.floes,
            )
        end
        empty!.(sim.model.floes.ghosts)

        # Add new pieces to the end of floe list
        append!(sim.model.floes, pieces_buffer)

        # Physical processes without ghost floes
        # Effects of ocean and atmosphere on ice and visa versa
        if (
            sim.coupling_settings.coupling_on &&
            mod(tstep, sim.coupling_settings.Δt) == 0
        )
            timestep_coupling!(
                sim.model,
                sim.Δt,
                sim.consts,
                sim.coupling_settings,
                sim.floe_settings,
            )
        end
        
        # Move and update floes based on collisions and ocean/atmosphere forcing
        timestep_floe_properties!(
            sim.model.floes,
            tstep,
            sim.Δt,
            sim.floe_settings,
        )
        # Fracture floes
        if sim.fracture_settings.fractures_on && mod(tstep, sim.fracture_settings.Δt) == 0
            max_floe_id =
                fracture_floes!(
                    sim.model.floes,
                    max_floe_id,
                    sim.rng,
                    sim.fracture_settings,
                    sim.floe_settings,
                    sim.Δt,
                )
        end

        # Weld floes
        if sim.weld_settings.weld_on
            weld_setting_idx = findfirst(
                x -> mod(tstep, x) == 0,
                sim.weld_settings.Δts
            )
            if !isnothing(weld_setting_idx)
                max_floe_id = Subzero.timestep_welding!(
                    sim.model.floes,
                    max_floe_id,
                    sim.model.grid,
                    sim.model.domain,
                    sim.weld_settings,
                    sim.floe_settings,
                    weld_setting_idx,
                    sim.Δt,
                )
            end
        end

        # What happens if floe tried to fuse with ghost floe?? 
        max_floe_id = 
            simplify_floes!(
                sim.model,
                max_floe_id,
                sim.simp_settings,
                sim.collision_settings,
                sim.floe_settings,
                sim.Δt,
                sim.rng,
            )
    end

    # h0 = real(sqrt.(Complex.((-2Δt * newfloe_Δt) .* hflx)))
    # mean(h0)
    return 
end

"""
    startup_sim(sim)

Required actions to setup simulation. For example, setting up the simulation
logger.
Inputs:
    sim                 <Simulation>
    logger              <AbstractLogger> logger for simulation - default is
                            Subzero logger
    messages_per_tstep  <Int> number of messages to print per timestep if using
                            default SubzeroLogger, else not needed
Outputs:
    None.
"""
function startup_sim(sim, logger = nothing, messages_per_tstep = 1)
    # Set up logger
    if isnothing(logger)
        logger = SubzeroLogger(sim, messages_per_tstep)
    end
    global_logger(logger)
    # Start sim notice
    sim.verbose && println(sim.name * " is running!")
    return
end

"""
    teardown_sim(sim)

Required actions to tear down simulation. For example, flushing the simulation's
logger and closing the stream.
Inputs:
    sim <Simulation>
Outputs:
    None.
"""
function teardown_sim(sim)
    # Finish logging
    logger = current_logger()
    if hasfield(typeof(logger), :stream)
        io = logger.stream
        flush(io)
        close(io)
    end
    # End sim notice
    sim.verbose && println(sim.name * " done running!")
    return
end

"""
    run!(sim; logger = nothing, messages_per_tstep = 1, start_tstep = 0)

Run given simulation and generate output for given writers.
Simulation calculations will be done with Floats of type T (Float64 of Float32).

Inputs:
    sim                 <Simulation> simulation to Run
    logger              <AbstractLogger> logger for simulation - default is
                            Subzero logger
    messages_per_tstep  <Int> number of messages to print per timestep if using
                            default SubzeroLogger, else not needed
    start_tstep         <Int> which timestep to start the simulation on
Outputs:
    None. The simulation will be run and outputs will be saved in the output
    folder. 
"""
function run!(sim; logger = nothing, messages_per_tstep = 1, start_tstep = 0)
    startup_sim(sim, logger, messages_per_tstep)
    tstep = start_tstep
    while tstep <= (start_tstep + sim.nΔt)
        # Timestep the simulation forward
        timestep_sim!(sim, tstep, start_tstep)
        tstep+=1
    end
    teardown_sim(sim)
    return
end

"""
    restart!(initial_state_fn, checkpointer_fn, new_nΔt, new_output_writers; start_tstep = 0)

Continue the simulation run started with the given initial state and floe file for an
additional `new_nΔt` timesteps and with the new output_writers provided. The simulation will
restart with a recorded timestep of `start_tstep`.

Note that this `restart!` function may not fit your needs and you may need to write your
own. This function is meant to act as a simplest case and as a template for users to write
their own restart functions. 
"""
function restart!(initial_state_fn, checkpointer_fn, new_nΔt, new_output_writers; start_tstep = 0)
    is = jldopen(initial_state_fn)
    cp = jldopen(checkpointer_fn)
    last_tstep = maximum(parse.(Int, keys(cp["ocean"])))

    # Remove any ghost floes from floe list
    new_floes = cp["floes"][string(last_tstep)]
    filter!(f -> f.ghost_id == 0, new_floes)
    empty!.(new_floes.ghosts)

    new_model = Model(
        is["sim"].model.grid, 
        cp["ocean"][string(last_tstep)], 
        cp["atmos"][string(last_tstep)], 
        is["sim"].model.domain, 
        new_floes,
    )

    new_simulation = Simulation(
        model = new_model,
        consts = is["sim"].consts,
        Δt = is["sim"].Δt,
        nΔt = new_nΔt,
        verbose = is["sim"].verbose,
        writers = new_output_writers,
        coupling_settings = is["sim"].coupling_settings,
        simp_settings = is["sim"].simp_settings,
    )
    run!(new_simulation; start_tstep = start_tstep)
    return
end