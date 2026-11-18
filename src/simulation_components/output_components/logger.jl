import Logging: shouldlog, min_enabled_level, catch_exceptions, handle_message
export SubzeroLogger

# Logger for Subzero - see documentation below
struct SubzeroLogger <: Logging.AbstractLogger
    stream::IO
    min_level::Logging.LogLevel
    message_limits::Dict{Union{String, Symbol},Int}
    messages_per_tstep::Int
end

"""
    SubzeroLogger(; sim = nothing, filename = " ", messages_per_tstep = 1)

Logger for Subzero. Logs unique messages `messages_per_tstep` times per timestep
to prevent overwhelming number of messages timesteps from multiple floes
triggering the same log event.

## _Fields_
- `stream::IO`: logs are written to this IO
- `min_level::Logging.LogLevel`: minimum log event level to write
- `message_limits::Dict{Any, Int}`: dictionary with message IDs for key whose values are the
    number of times that message can still be written in current timestep - 
    current timestep is stored in same dictionary under key 'tstep'
- `messages_per_tstep::Int`: maximum number of times a given message should be written per timestep

## _Keyword arguments_
- `sim::Simulation`: Subzero simulation - used to get simulation name (Default = nothing), which is then used to name log file
- `filename::String`: if simulation isn't provided, a filename (+ path) must be provided to save log file
- `messages_per_tstep::Int`: maximum number of times a given message should be written per timestep
"""
function SubzeroLogger(; sim = nothing, min_log_level = nothing, messages_per_tstep = 1)
    filename = if !isnothing(sim)
        "./log/$(sim.name).log"
    else
        "./log/sim.log"
    end
    if isnothing(min_log_level)
        min_log_level = Logging.Info
    end
    # Create folder and file
    logfolder = dirname(filename)
    mkpath(logfolder)
    isfile(filename) && rm(filename, force=true)
    # Create logger
    return SubzeroLogger(
        open(filename, "w+"),
        min_log_level,
        Dict("tstep" => 0),
        messages_per_tstep,  # number of messages per timestep (per message)
    )
end

# Required logging functions
shouldlog(logger::SubzeroLogger, level, _module, group, id) = logger.min_level ≤ level
min_enabled_level(logger::SubzeroLogger) = logger.min_level
catch_exceptions(logger::SubzeroLogger) = false

# Returns string with log event name given log event level
function level_to_string(level)
    level == Logging.Error && return "ERROR"
    level == Logging.Warn  && return "WARN "
    level == Logging.Info  && return "INFO "
    level == Logging.Debug && return "DEBUG"
    return string(level)
end

#= 
Function that determines if log event should be written to file depending on how
many times that event has been written to file in current timestep.


This is called when a log macro is called (e.g. @warn), not explicitly by
the user. Additionally, it is not threadsafe so a message may be written
more times than `messages_per_tstep`, but it should be in the ballpark.
Putting a lock would slow down logging and isn't worth it given that this
problem only records a few extra log events. 
=#
function handle_message(
    logger::SubzeroLogger,
    level,
    message,
    _module,
    group,
    id,
    filepath,
    line;
    tstep = nothing,  # current timestep
    kwargs...,
)
    # Each message should only write `messages_per_tstep` times per timestep
    if !isnothing(tstep)
        if tstep > logger.message_limits["tstep"]  # new timestep
            for id in keys(logger.message_limits)  # reset log event count
                logger.message_limits[id] = logger.messages_per_tstep
            end
            logger.message_limits["tstep"] = tstep
        end
        remaining = get!(
            logger.message_limits,
            id,
            logger.messages_per_tstep,
        )
        remaining > 0 || return nothing
        logger.message_limits[id] = remaining - 1  # decrease log event count
    end

    buf = IOBuffer()
    iob = IOContext(buf, logger.stream)

    # Log type and message
    level_name = level_to_string(level)
    formatted_message = "$level_name $message"
    # Add simulation timestep to message
    if !isnothing(tstep)
        formatted_message *= " --> timestep $tstep"
    end
    # Add wall clock time to message
    msg_timestamp = Dates.format(Dates.now(), "[yyyy/mm/dd HH:MM:SS.sss]")
    formatted_message *= " $msg_timestamp"
    # Add log location in code to message
    file_name   = something(filepath, "nothing")
    line_number = something(line, "nothing")
    formatted_message *= " --> $file_name:$line_number"
    # Write message
    println(iob, formatted_message)
    write(logger.stream, take!(buf))
    return
end