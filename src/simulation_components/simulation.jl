export Simulation, timestep_sim!, run!, restart!

const ΔT_DEF = "`Δt::Int`:length of timestep in integer seconds"

"""
    Simulation{FT, MT, CT, PT, ST, RT, OT}

Simulation which holds a model and the parameters, settings, and output writers needed for running the simulation.

Only keyword arguments are used! 

## _Fields_ / _Keyword Arguments_
### _General_
- $MODEL_DEF
- `consts::Constants{FT}`: Constants used in Simulation (default = Constants())
- `rng::RT`: Random number generator (default = Xoshiro())
- `verbose::Bool`: String output printed during run (Default = false)
- `name::String`: Simulation name for printing/saving (Default = "sim")
### _Timesteping Information_
- `Δt::Int`: Simulation timestep in seconds
- `nΔt::Int`: Total timesteps simulation runs for
### _Physical Processes_
- `floe_settings::FloeSettings{FT, PT, ST}`: Settings that control floe size/mass/etc - no default!
- `coupling_settings::CouplingSettings`: Settings that control coupling between floes/ocean/atmosphere (Default = CouplingSettings())
- `collision_settings::CollisionSettings{FT}`: Settings that control floe collisions with other floes and the domain (Default = CollisionSettings())
- `fracture_settings::FractureSettings{CT}`: Settings that control floe fracturing (Default = FractureSettings())
- `simp_settings::SimplificationSettings{FT}`: Settings that control the simplification of floes (Default = SimplificationSettings())
- `ridgeraft_settings::RidgeRaftSettings{FT}`: Settings that control floe ridging and rafting (Default = RidgeRaftSettings())
- `weld_settings::WeldSettings{FT}`: Settings that control floe welding (Default = WeldSettings())
### _Output Writers_
- `writers::OT`: Simulation output writers (Default = OutputWriters())

!!! note
    Unlike almost all of the other constructors that are used to create the structs that are fed into the `Simulation`
    constructor, the `Simulation` constructor doesn't have a `FT` argument to set the Float-type of the simulation. This
    is because it will simply use the Float-type of all of the fields, which _must_ all match! 

!!! note
    If a `FloeSettings` object is required since it is also an input to the [`_initialize_floe_field!`](@ref) functions.
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
    model::MT                               # Model to simulate
    consts::Constants{FT} = Constants()     # Constants used in Simulation
    rng::RT = Xoshiro()                     # Random number generator 
    verbose::Bool = false                   # String output printed during run
    name::String = "sim"                    # Simulation name for printing/saving
    # Timesteps ----------------------------------------------------------------
    Δt::Int                     # Simulation timestep (seconds)
    nΔt::Int                    # Total timesteps simulation runs for
    # Physical Processes -------------------------------------------------------
    floe_settings::FloeSettings{FT, PT, ST}
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

## _Positional arguments_
- $SIM_DEF
- `tstep::Int`: simulation's current timestep
- `start_tstep::Int`: timestep simulation started on (Default = 0)

## _Returns_
-  None. Simulation advances by one timestep. 

!!! note
    The order of the function calls within `timestep_sim!` matter quite a lot! Swapping the order can break things as certian fields are cleared at the end of function calls. This migth be worth debugging at some point.

!!! note
    This function must be updated to add any new functionalitites. It might be worth modularizing somehow to smooth that process over and make it easier to add new science functionality.
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

# Required actions to setup simulation. Right now, this only entails printing a notice.
function startup_sim(sim)
    # Start sim notice
    sim.verbose && println(sim.name * " is running!")
    return
end

# Required actions to tear down simulation. Right now, this only entails printing a notice.
function teardown_sim(sim)
    # End sim notice
    sim.verbose && println(sim.name * " done running!")
    return
end

"""
    run!(sim; logger, messages_per_tstep, start_tstep)

Run given simulation and generate output for given output writers.
Simulation calculations will be done with Floats of type FT (Float64 of Float32).

## _Positional arguments_
- $SIM_DEF

## _Keyword arguments_
- `logger::AbstractLogger`: logger for simulation (Default = Nothing, which triggers use of [`SubzeroLogger`](@ref)).
    The logger is only used while the simulation runs; the global logger is not changed.
    A default `SubzeroLogger` is closed at the end of the run, but a given logger is not.
- `messages_per_tstep::Int`"` number of messages to print per timestep if using default SubzeroLogger, else not needed (Default = 1)
- `start_tstep::Int`: which timestep to start the simulation on (Default = 0)

## _Returns_
- None. The simulation will be run and outputs will be saved in the output folder. 
"""
function run!(sim; logger = nothing, messages_per_tstep = 1, start_tstep = 0)
    # Only use the logger within this function. Replacing the global logger breaks
    # packages that log from an older world age, such as GPUCompiler.
    own_logger = isnothing(logger)
    if own_logger
        logger = SubzeroLogger(; sim, messages_per_tstep)
    end
    try
        with_logger(logger) do
            startup_sim(sim)
            tstep = start_tstep
            while tstep <= (start_tstep + sim.nΔt)
                # Timestep the simulation forward
                timestep_sim!(sim, tstep, start_tstep)
                tstep+=1
            end
            teardown_sim(sim)
        end
    finally
        if own_logger
            close(logger.stream)
        end
    end
    return
end

"""
    restart!(initial_state_fn, checkpointer_fn, new_nΔt, new_output_writers; start_tstep)

Continue the simulation run started with the given initial state and floe file for an
additional `new_nΔt` timesteps and with the new output_writers provided. The simulation will
restart with a recorded timestep of `start_tstep`.

Note that this `restart!` function may not fit your needs and you may need to write your
own. This function is meant to act as a simplest case and as a template for users to write
their own restart functions. 

## _Positional arguments_
    - `initial_state_fn::String`: file path to previously run simulation's initial state file
    - `checkpointer_fn::String`: file path to previously run simulation's checkpointer file.
    The simulation will restart right after the checkpoint captured in this file
    - `new_output_writers::OutputWriters`: new output writers for the new simulation - new ones are required to the output writes to new, unique files.

## _Keyword arguments_
    - `start_tstep::Int`: which timestep to start the simulation on (Default = 0)

## _Returns_
    - None. The simulation will be run and outputs will be saved in the output folder. 
"""
function restart!(initial_state_fn, checkpointer_fn, new_nΔt, new_output_writers; start_tstep = 0)
    is = jldopen(initial_state_fn)
    cp = jldopen(checkpointer_fn)
    last_tstep = maximum(parse.(Int, keys(cp["ocean"])))

    # Remove any ghost floes from floe list
    new_floes = cp["floes"][string(last_tstep)]
    filter!(f -> f.ghost_id == 0, new_floes)
    empty!.(new_floes.ghosts)

    new_model = Model(;
        grid = is["sim"].model.grid, 
        ocean = cp["ocean"][string(last_tstep)], 
        atmos = cp["atmos"][string(last_tstep)], 
        domain = is["sim"].model.domain, 
        floes = new_floes,
    )

    new_simulation = Simulation(
        model = new_model,
        consts = is["sim"].consts,
        Δt = is["sim"].Δt,
        nΔt = new_nΔt,
        verbose = is["sim"].verbose,
        writers = new_output_writers,
        floe_settings = is["sim"].floe_settings,
        coupling_settings = is["sim"].coupling_settings,
        collision_settings = is["sim"].collision_settings,
        fracture_settings = is["sim"].fracture_settings,
        simp_settings = is["sim"].simp_settings,
        ridgeraft_settings = is["sim"].ridgeraft_settings,
        weld_settings = is["sim"].weld_settings,
    )
    run!(new_simulation; start_tstep = start_tstep)
    return
end

# Pretty printing for Simulation showing key dimensions
function Base.show(io::IO, sim::Simulation)
    overall_summary = "Simulation"
    timestep_summary = "Timestep: $(sim.Δt) seconds"
    runtime_summary = "Runtime: $(sim.nΔt) timesteps"
    print(io, overall_summary, "\n",
        "  ⊢", timestep_summary, "\n",
        "  ⊢", runtime_summary, "\n",
        "  ⊢RNG: ", sim.rng, "\n",
        "  ⊢verbose: ", sim.verbose, "\n",
        "  ⊢model\n",
        "  ⊢consts\n",
        "  ⊢floe_settings\n",
        "  ⊢collision_settings\n",
        "  ∟ ...")
end
