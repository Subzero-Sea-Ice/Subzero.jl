# # Basic plotting functions (and stub functions) for Subzero simulations
export plot_sim, prettytime

# constants used for documentation
const FLOE_FN_DEF = "floe outputwriter output file path and name"
const INITIAL_STATE_FN_DEF = "initial state outputwriter output file path and name"
const MP4_OUTPUT_FN = "output video file path and name (should end with .mp4)"

#= 
## What plotting functionality is availible for Subzero simulations?

We provide very basic plotting functionalities for Subzero simulations. Given plotting is so
customizable, we just wanted to provide a roadmap for how a user might plot a simulation and
provide a few utility functions. The user should use our plotting code for early
experimentation and as an exampe of how to create animations with `Makie`.

Below are the stubs of functions that require `CairoMakie` and `GeoInterfaceMakie`. The
actual code for these functions is in `ext/SubzeroMakieExt`. These functions are only loaded
when the user loads `CairoMakie` and `GeoInterfaceMakie`. This prevents unneccesarily long
compiling times when the user just wants to run simulations. They can then load `CairoMakie`
and `GeoInterfaceMakie` once post-simulation runs and create all of the needed videos at
once, rather than loading `CairoMakie` prior to each simulation.

We also provide basic plotting utility functions, like `prettytime`, here as they do
not depend on any of the `Makie packages` and therefore don't need to be abstracted away
into the extension.
=#

# Stub functions that depend on CairoMakie implemented in ext/SubzeroMakieExt.jl
"""
    plot_sim(floe_fn, initial_state_fn, Δt, output_fn; max_side_pixels = 800)

Basic plotting of a simulation using the simulation's floe and initial state files. This
function is meant for basic plotting during testing and as an example of how to create
a video using Makie. The user can (and should!) write their own plotting code to their own
specifications and needs.

This code does not have an underlying ocean, but there are comments on how to add an ocean
heatmap and colorbar within the source code.

## Arguments:
- `floe_fn::String`: $FLOE_FN_DEF
- `initial_state_fn::String`: $INITIAL_STATE_FN_DEF
- $ΔT_DEF
- `output_fn::String`: $MP4_OUTPUT_FN
"""
function plot_sim end

# Utility functions that don't depend on CairoMakie
"""
    prettytime(t)

Turn time in seconds into units of minutes, hours, days, or years as appropriate
## _Positional Arguments_:
- `t::Real`: number of seconds
## Outputs:
- `::String`: number of seconds converted to a string value in minutes, hours, days, or years with units
## Note:
    This code was modified from the this [source code](https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl).
"""
function prettytime(t)
    minute = 60
    hour = 60 * minute
    day = 24 * hour
    year = 365 * day

    iszero(t) && return "0 seconds"
    if t < minute
        value = floor(Int, t)
        units = value == 1 ? "second" : "seconds"
    elseif t < hour
        value = floor(Int, t / minute)
        units = value == 1 ? "minute" : "minutes"
    elseif t < day
        value = floor(Int, t / hour)
        units = value == 1 ? "hour" : "hours"
    elseif t < year
        value = floor(Int, t / day)
        units = value == 1 ? "day" : "days"
    else
        value = floor(Int, t / year)
        units = value == 1 ? "year" : "years"
    end
    return @sprintf("%d %s", value, units)
end