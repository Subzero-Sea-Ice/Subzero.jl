import Logging: shouldlog, min_enabled_level, catch_exceptions, handle_message

"""
    SubzeroLogger

Logger for Subzero. Logs unique messages `messages_per_tstep` times per timestep
to prevent overwhelming number of messages timesteps from multiple floes
triggering the same log event.
Fields:
- stream: ogs are written to this IO
- min_level: minimum log event level to write
- message_limits: dictionary with message IDs for key whose values are the
    number of times that message can still be written in current timestep - 
    current timestep is stored in same dictionary under key 'tstep'
- messages_per_tstep: maximum number of times a given message should be written
    per timestep
"""
struct SubzeroLogger <: Logging.AbstractLogger
    stream::IO
    min_level::Logging.LogLevel
    message_limits::Dict{Any,Int}
    messages_per_tstep::Int

    function SubzeroLogger(
        stream,
        min_level,
        message_limits,
        messages_per_tstep,
    )
        if !haskey(message_limits, "tstep")
            message_limits["tstep"] = 0  # timestep key/value pair must exist
        end
        new(
            stream,
            min_level,
            message_limits,
            messages_per_tstep,
        )
    end
end

"""
    SubzeroLogger(sim, messages_per_tstep = 1)

Constructor from Subzero logger.
Inputs:
    filename            <String> file path to file to write log events into
    messages_per_tstep  <Int> maximum number of times a given message should be
                            written per timestep
Outputs:
    Subzero logger that saves log to file with the same name as the simulation's
    name field and optional keyword arguments set. 
"""
function SubzeroLogger(filename::String, messages_per_tstep = 1)
    # Create folder and file
    logfolder = dirname(filename)
    mkpath(logfolder)
    isfile(filename) && rm(filename, force=true)
    # Create logger
    return SubzeroLogger(
        open(filename, "w+"),
        Logging.Info,
        Dict{Any,Int}(),
        messages_per_tstep,  # number of messages per timestep (per message)
    )
end

"""
    SubzeroLogger(sim; messages_per_tstep = 1)

Created Subzero logger and writes log events to log file in current directory to
file with the same name as the simulation's name field.
Inputs:
    sim                 <Simulation>
    messages_per_tstep  <Int> maximum number of times a given message should be
                        written per timestep
Outputs:
    Subzero logger that saves log to file with the same name as the simulation's
    name field and optional keyword arguments set. 
"""
SubzeroLogger(sim, messages_per_tstep = 1) =
    SubzeroLogger(
        "./log/$(sim.name).log",
        messages_per_tstep,
    )

# Required logging functions
shouldlog(logger::SubzeroLogger, level, _module, group, id) = true
min_enabled_level(logger::SubzeroLogger) = Logging.Info
catch_exceptions(logger::SubzeroLogger) = false

"""
    level_to_string(level)

Returns string with log event name given log event level
"""
function level_to_string(level)
    level == Logging.Error && return "ERROR"
    level == Logging.Warn  && return "WARN "
    level == Logging.Info  && return "INFO "
    level == Logging.Debug && return "DEBUG"
    return string(level)
end

"""
    handle_message(
        logger::SubzeroLogger,
        level,
        message,
        _module,
        group,
        id,
        filepath,
        line;
        tstep = nothing,
        kwargs...,
    )

Function that determines if log event should be written to file depending on how
many times that event has been written to file in current timestep.

Note:
    This is called when a log macro is called (e.g. @warn), not explicitly by
    the user. Additionally, it is not threadsafe so a message may be written
    more times than `messages_per_tstep`, but it should be in the ballpark.
    Putting a lock would slow down logging and isn't worth it given that this
    problem only records a few extra log events. 
"""
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