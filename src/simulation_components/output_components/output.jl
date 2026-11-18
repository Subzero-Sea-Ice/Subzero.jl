export AbstractOutputWriter, CheckpointOutputWriter, GridOutputWriter,  FloeOutputWriter, InitialStateOutputWriter, OutputWriters

const TSTEP_DEF = "`tstep::Int`: current simulation timestep"

"""
    AbstractOutputWriter

An abstract type for output writers that provide data from simulation runs.

Right now, there are four types of output writers: [`InitialStateOutputWriter`](@ref), [`CheckpointOutputWriter`](@ref), [`FloeOutputWriter`](@ref), and [`GridOutputWriter`](@ref).

They currently don't dispatch off of anything.
"""
abstract type AbstractOutputWriter end

const FILEPATH_DEF = "`filepath::String`: filepath to the output files"
const OVERWRITE_DEF = "`overwrite::Bool `: if true, existing files with the same `filepath` will be overwritten (Default = false)"
const DIR_DEF = "`dir::String`: directory the files should be saved in"
const FILENAME_DEF = "`filename::String`: name of output file"
const JLD2_KW_DEF = "`jld2_kw::Dict`: dictionary of keywords for the `jldopen` function"


# Concrete subtype of AbstractOutputWriter - see documentation below
struct InitialStateOutputWriter<:AbstractOutputWriter
    filepath::String     # Filename for output file
    overwrite::Bool      # Remove existing files if their filenames conflict
end

"""
    InitialStateOuputWriter <: AbstractOutputWriter

Concrete subtype of AbstractOutputWriter that records the intial state of the
simulation so that you can re-load it later. Writes JLD2 file with initial simulation
state to the filepath specified. If overwrite is true, and there is a file of the given
name at the filepath, that file will be overwritten.

## _Fields_
- $FILEPATH_DEF
- $OVERWRITE_DEF

## _Keyword arguments_
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "initial_state.jld2")
- $OVERWRITE_DEF
- $JLD2_KW_DEF

Here is how to construct an `InitialStateOutputWriter`:

    InitialStateOutputWriter(; kwargs...)
"""
function InitialStateOutputWriter(
    ;
    dir = ".",
    filename = "initial_state.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}()
)
    filepath = _initialize_jld2_file!(dir, filename, overwrite, [], jld2_kw)
    return InitialStateOutputWriter(filepath, overwrite)
end

"""
    InitialStateOutputWriter <: AbstractOutputWriter

`InitialStateOutputWriter` can also be created using an existing `InitialStateOutputWriter`,
copying all fields unless the new field values are explicity provided through keyword arguments.

## _Positional arguments_
- `writer::InitialStateOutputWriter`: existing `InitialStateOutputWriter`

## _Keyword arguments_
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "initial_state.jld2")
- $OVERWRITE_DEF
- kwargs...


"""
InitialStateOutputWriter(writer::InitialStateOutputWriter;
    dir = dirname(writer.filepath),
    filename = basename(writer.filepath),
    overwrite = writer.overwrite,
    kwargs...,
) = InitialStateOutputWriter(;
    dir = dir,
    filename = filename,
    overwrite = overwrite,
    kwargs...,
)

# Concrete subtype of AbstractOutputWriter - see documentation below
struct CheckpointOutputWriter<:AbstractOutputWriter
    Δtout::Int              # Number of timesteps between checkpoint outputs
    filepath::String        # Filename for output file
    overwrite::Bool         # Remove existing files if their filenames conflict
end

"""
    CheckpointOutputWriter <: AbstractOutputWriter

Checkpoint subtype of AbstractOutputWriter that holds information for outputting
checkpoint information used for restarting the simulation from a point where the
writer saved data. Checkpoint data is saved every `Δtout` timesteps to filepath.
If the given file doesn't end in ".jld2", the extension will be appended. If
overwrite is true then if there is already a file of the given name, it will be
overwriten. Else it will thrown an error. 

## _Fields_
- `Δtout::Int`: number of timesteps between output
- $FILEPATH_DEF
- $OVERWRITE_DEF

## _Positional arguments_
- `Δtout::Int`: number of timesteps between output

## _Keyword arguments_
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "checkpoint.jld2")
- $OVERWRITE_DEF
- $JLD2_KW_DEF

Here is how to construct an `CheckpointOutputWriter`:

    CheckpointOutputWriter(Δtout; kwargs...)
"""
function CheckpointOutputWriter(
    Δtout;
    dir = ".",
    filename = "checkpoint.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
    filepath = _initialize_jld2_file!(
        dir,
        filename,
        overwrite,
        [:ocean, :atmos, :floes],
        jld2_kw,
    )
    return CheckpointOutputWriter(Δtout, filepath, overwrite)
end

"""
    CheckpointOutputWriter <: AbstractOutputWriter

A `CheckpointOutputWriter` can also be created from an existing `CheckpointOutputWriter`, copying all fields unless new field
values are explicity provided either as the optional argument Δtout or through keyword
arguments.

## _Positional arguments_
- `writer::InitialStateOutputWriter`: OPTIONAL argument - if present, the `dir`, `filename`, and `overwrite` of the provided
    writer will be used for the new writer, rather than the above argument
- `Δtout::Int`: number of timesteps between output

## _Keyword arguments_
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "checkpoint.jld2")
- $OVERWRITE_DEF
- kwargs...

Here is how to construct a new `CheckpointOutputWriter` using an existing `CheckpointOutputWriter`:

    CheckpointOutputWriter(existing_checkpointer, Δtout; kwargs...)
"""
CheckpointOutputWriter(writer::CheckpointOutputWriter, Δtout = writer.Δtout;
    dir = dirname(writer.filepath),
    filename = basename(writer.filepath),
    overwrite = writer.overwrite,
    kwargs...,
) = CheckpointOutputWriter(Δtout;
    dir = dir,
    filename = filename,
    overwrite = overwrite,
    kwargs...,
)

# Concrete subtype of AbstractOutputWriter - see documentation below
struct FloeOutputWriter<:AbstractOutputWriter
    Δtout::Int                  # Number of timesteps between floe outputs
    outputs::Vector{Symbol}     # Floe fields to output
    filepath::String            # Filename for output file
    overwrite::Bool             # Remove existing files if filenames conflict
end

"""
    FloeOutputWriter <: AbstractOutputWriter

Floe subtype of AbstractOutputWriter that holds information for outputting floe
information from model throughout simulation. Output will be saved to the file
defined by `filename` every `Δtout` timesteps. Only outputs within the outputs list will
be saved. The `outputs` field takes in a list of symbols corresponding to floe fields.
For example, if you want the floe output writer to output the floes centroid and polygons
then `outputs = [:centroid, :poly]`. If you want all floe fields then you can simply omit
the outputs field all together and all floe fields will be output.

File will be saved as a JLD2 file to filepath. If the given file
doesn't end in ".jld2", the extension will be appended. If `overwrite` is true then
if there is already a file of the given name, it will be overwriten. Else it
will thrown an error. 

## _Fields_
- `Δtout::Int`: number of timesteps between output
- `outputs::Vector{Symbol}`: list of Floe field names (as symbols) to output
- $FILEPATH_DEF
- $OVERWRITE_DEF

## _Positional arguments_
- `Δtout::Int`: number of timesteps between output

## _Keyword arguments_
- `outputs::Vector{Symbol}`: list of Floe field names (as symbols) to output (Default outputs all `Floe` fields)
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "floes.jld2")
- $OVERWRITE_DEF
- $JLD2_KW_DEF

!!! note
    If you have Periodic walls, and thus ghost floes in your simulation, these will also be saved by the `FloeOutputWriter`.
    If you want to exclude these floes from your analysis or when otherwise using the `FloeOutputWriter` output, you can do
    so by only including floes with a `ghost_id = 0` when post-processing. 

Here is how to construct an `FloeOutputWriter`:

    FloeOutputWriter(Δtout; kwargs...)
"""
function FloeOutputWriter(
    Δtout;
    outputs = collect(fieldnames(Floe)),
    dir = ".",
    filename = "floes.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
    filepath = _initialize_jld2_file!(dir, filename, overwrite, outputs, jld2_kw)
    return FloeOutputWriter(Δtout, outputs, filepath, overwrite)
end

"""
    FloeOutputWriter <: AbstractOutputWriter

A `FloeOutputWriter` can also be created using an existing `FloeOutputWriter`, copying all fields unless the new field values are explicity provided through keyword arguments.

## _Positional arguments_
- `writer::InitialStateOutputWriter`: OPTIONAL argument - if present, the `dir`, `filename`, and `overwrite` of the provided 
    writer will be used for the new writer, rather than the above argument
- `Δtout::Int`: number of timesteps between output

## _Keyword arguments_
- `outputs::Vector{Symbol}`: list of Floe field names (as symbols) to output
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "floes.jld2")
- $OVERWRITE_DEF
- kwargs...

Here is how to construct a new `FloeOutputWriter` using an existing `FloeOutputWriter`:

    FloeOutputWriter(existing_floewriter, Δtout; kwargs...)
"""
FloeOutputWriter(writer::FloeOutputWriter, Δtout = writer.Δtout;
    outputs = writer.outputs,
    dir = dirname(writer.filepath),
    filename = basename(writer.filepath),
    overwrite = writer.overwrite,
    kwargs...,
) = FloeOutputWriter(Δtout;
    outputs = outputs,
    dir = dir,
    filename = filename,
    overwrite = overwrite,
    kwargs...,
)

# Concrete subtype of AbstractOutputWriter - see documentation below
struct GridOutputWriter{FT<:AbstractFloat}<:AbstractOutputWriter
    outputs::Vector{Symbol}
    Δtout::Int                # Number of timesteps between grid outputs
    filepath::String          # Filename for output file
    overwrite::Bool           # Remove existing files if filenames conflict
    xg::Vector{FT}            # Grid lines in x-direction for ouput calculations
    yg::Vector{FT}            # Grid lines in y-direction for ouput calculations
    data::Array{FT, 3}        # Data to write to file stored here
    average::Bool
end


# Returns list of symbols that represent calculations available in `calc_eularian_grid` to average floe data.
function get_known_grid_outputs()
    return Set([
        :u_grid,
        :v_grid,
        :dudt_grid,
        :dvdt_grid,
        :overarea_grid,
        :mass_grid,
        :area_grid,
        :height_grid,
        :si_frac_grid,
        :stress_xx_grid,
        :stress_yx_grid, 
        :stress_xy_grid,
        :stress_yy_grid,
        :stress_eig_grid,
        :strain_ux_grid,
        :strain_vx_grid,
        :strain_uy_grid,
        :strain_vy_grid,
        # Adding bond stress too
        :stress_bond_xx_grid,
        :stress_bond_yx_grid, 
        :stress_bond_xy_grid,
        :stress_bond_yy_grid,
    ])
end

"""
    GridOutputWriter{FT} <: AbstractOutputWriter

Grid subtype of AbstractOutputWriter that holds information for outputting floe
data on the grid. This output does not need to be on the grid defined for the
model. This grid can be coarser, or more fine, as defined by the xg and yg
fields. Output on this scale will be saved to the file defined by filepath every
Δtout timesteps. Data will be collected in the data field during calculation for
easier writing to the NetCDF file. Only outputs within the outputs list will be
saved. There is a limited number of floe outputs that can be calculated by the
`GridOutputWriter`.

## _Fields_
- `outputs::Vector{Symbol}`: list of field names (as symbols) to output - call `get_known_grid_outputs()` to see options
- `Δtout::Int`: number of timesteps between output
- $FILEPATH_DEF
- $OVERWRITE_DEF
- `xg::Vector{FT}`: x-grid points on output grid
- `yg::Vector{FT}`: y-grid points on output grid
-  `data::Array{FT, 3}`: 3D output grid of size `nx`, `ny`, and `nt` where nt is the number of timesteps saved
- `average::Bool`: if true, average gridded data over timesteps between outputs, else just calculate at output timestep

## _Positional arguments_
- `Δtout::Int`: number of timesteps between output
- $GRID_DEF
- `dims::Tuple{Int, Int}`: output new grid dimensions for these calculations (rows -> ny, cols -> nx)

## _Keyword arguments_
- `outputs::Vector{Symbol}`: list of field names (as symbols) to output - call `get_known_grid_outputs()` to see options
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "gridded_data.nc")
- $OVERWRITE_DEF
- `average::Bool`: if true, average gridded data over timesteps between outputs, else just calculate at output timestep

Here is how to construct an `GridOutputWriter`:

    GridOutputWriter(Δtout, grid, dims; kwargs...)
"""
function GridOutputWriter{FT}(
    Δtout,
    grid::AbstractRectilinearGrid,
    dims;
    outputs::Vector{Symbol} = collect(get_known_grid_outputs()),
    dir = ".",
    filename = "gridded_data.nc",
    overwrite = false,
    average = false,
) where {FT <: AbstractFloat}
    # Check for known outputs - need routine to calculate in calc_eulerian_data
    known_grid_outputs = get_known_grid_outputs()
    remove_idx = []
    for i in eachindex(outputs)
        if !(outputs[i] in known_grid_outputs)
            @warn "$(outputs[i]) is not a known grid output. It will not be \
                output with the GridOutputWriter."
            push!(remove_idx, i)
        end
    end
    deleteat!(outputs, remove_idx)

    # Define output grid
    xg = collect(grid.x0:(grid.xf-grid.x0)/dims[1]:grid.xf)
    yg = collect(grid.y0:(grid.yf-grid.y0)/dims[2]:grid.yf)

    # Create file path
    filepath = _initialize_netcdf_file!(FT, dir, filename, overwrite, outputs, xg, yg)

    # Data output container
    data = zeros(FT, length(xg) - 1, length(yg) - 1, length(outputs))
    return GridOutputWriter{FT}(outputs, Δtout, filepath, overwrite, xg, yg, data, average)
end

# syntactic sugar for GridOutputWriter
GridOutputWriter(::Type{FT}, args...; kwargs...) where {FT <: AbstractFloat} =
    GridOutputWriter{FT}(args...; kwargs...)

# syntactic sugar for GridOutputWriter - default to Float64
GridOutputWriter(args...; kwargs...) =
    GridOutputWriter{Float64}(args...; kwargs...)

"""
    GridOutputWriter <: AbstractOutputWriter

`GridOutputWriter` can also be created using existing `GridOutputWriter`, copying all fields unless the new field values are explicity provided through keyword arguments.

## _Positional arguments_
- `writer::InitialStateOutputWriter`: OPTIONAL argument - if present, the `dir`, `filename`, and `overwrite` of the provided 
    writer will be used for the new writer, rather than the above argument
- `Δtout::Int`: number of timesteps between output

## _Keyword arguments_
- $GRID_DEF
- `dims::Tuple{Int, Int}`: output new grid dimensions for these calculations (rows -> ny, cols -> nx)
- `outputs::Vector{Symbol}`: list of field names (as symbols) to output - call `get_known_grid_outputs()` to see options
- $DIR_DEF (Default = ".")
- $FILENAME_DEF (Default = "gridded_data.nc")
- $OVERWRITE_DEF
- `average::Bool`: if true, average gridded data over timesteps between outputs, else just calculate at output timestep

!!! note
    The argument/field `average` currently doesn't do anything! Only instantaneous values can be saved. The argument was added
    so that this would be an easy change in the future. 

Here is how to construct a new `GridOutputWriter` using an existing `GridOutputWriter`:

    GridOutputWriter(existing_gridwriter, Δtout; kwargs...)
"""
GridOutputWriter(writer::GridOutputWriter, Δtout = writer.Δtout;
    grid = writer.grid,
    dims = writer.dims,
    outputs = writer.outputs,
    dir = dirname(writer.filepath),
    filename = basename(writer.filepath),
    overwrite = writer.overwrite,
    average = writer.average
) = GridOutputWriter(
    Δtout,
    grid,
    dims;
    outputs = outputs,
    dir = dir,
    filename = filename,
    overwrite = overwrite,
    average = average,
)

# Struct to hold lists of output writer types - see documentation below
struct OutputWriters{
    IW<:StructVector{<:InitialStateOutputWriter},
    FW<:StructVector{<:FloeOutputWriter},
    GW<:StructVector{<:GridOutputWriter},
    CW<:StructVector{<:CheckpointOutputWriter},
} 
    initialwriters::IW
    floewriters::FW
    gridwriters::GW
    checkpointwriters::CW
end

"""
    OutputWriters{FT}

Structure to hold all types of output writers a user might want. Fields are
vectors so that more than one of each type of output writer can be defined, and
so that a default OutputWriter object doesn't create default output writer
fields, which would create files, but rather empty lists of output writers.
If any of the fields is not provided, the default is just an empty list. 

## _Fields_
-  `initialwriters::StructVector{<:InitialStateOutputWriter}`: vector of `InitialStateOutputWriter`s
- `floewriters::StructVector{<:FloeOutputWriter}`: vector of `FloeOutputWriter`s
- `gridwriters::StructVector{<:GridOutputWriter}`: vector of `GridOutputWriter`s
- `checkpointwriters::StructVector{<:CheckpointOutputWriter}`: vector of `CheckpointOutputWriter`s

## _Positional arguments_
- All output writers (one per argument) - they will be aggregated and put into lists by type
"""
function OutputWriters(args...)
    initialwriters = Vector{InitialStateOutputWriter}()
    floewriters = Vector{FloeOutputWriter}()
    gridwriters = Vector{GridOutputWriter}()
    checkpointwriters = Vector{CheckpointOutputWriter}()
    for writer in args
        wt = typeof(writer)
        if wt <: InitialStateOutputWriter
            push!(initialwriters, writer)
        elseif wt <: FloeOutputWriter
            push!(floewriters, writer)
        elseif wt <: GridOutputWriter
            push!(gridwriters, writer)
        elseif wt <: CheckpointOutputWriter
            push!(checkpointwriters, writer)
        end
    end
    return OutputWriters(
        StructVector(initialwriters),
        StructVector(floewriters),
        StructVector(gridwriters),
        StructVector(checkpointwriters)
    )
end

# Pretty printing for OutputWriters showing key dimensions
function Base.show(io::IO, writers::OutputWriters)
    overall_summary = "OutputWriters"
    init_writer_summary = "$(length(writers.initialwriters)) InitialStateOuputWriter(s)"
    checkpoint_writer_summary = "$(length(writers.checkpointwriters)) CheckpointOutputWriter(s)"
    floe_writer_summary = "$(length(writers.floewriters)) FloeOutputWriter(s)"
    grid_writer_summary = "$(length(writers.gridwriters)) GridOutputWriter(s)"
    print(io, overall_summary, "\n",
        "  ⊢", init_writer_summary, "\n",
        "  ⊢", checkpoint_writer_summary, "\n",
        "  ⊢", floe_writer_summary, "\n",
        "  ∟", grid_writer_summary)
end

#----------------------- Write Data -----------------------#

"""
    write_data!(sim, tstep, start_tstep)

Writes data for the simulation's writers that are due to write at given tstep.

## _Positional arguments_
- $SIM_DEF
- $TSTEP_DEF
- `start_tstep::Int`: starting timestep of the simulation
"""
function write_data!(sim, tstep, start_tstep)
    # write initial state on first timestep
    if tstep == start_tstep
        write_init_state_data!(sim)
    end
    # Write checkpoint data
    write_checkpoint_data!(sim, tstep)
    # Write floe data
    write_floe_data!(
        sim.writers.floewriters,
        sim.model.floes,
        tstep,
    )
    # Write grid data
    write_grid_data!(
        sim.writers.gridwriters,
        sim.model.floes,
        sim.model.domain.topography,
        tstep,
    )
    return
end

"""
    write_init_state_data!(sim, tstep)

Save initial simulation state.

## _Positional arguments_
- $SIM_DEF
- $TSTEP_DEF
"""
function write_init_state_data!(sim)
    for filepath in sim.writers.initialwriters.filepath
        jldopen(filepath, "a") do file
            file["sim"] = sim
        end
    end
    return
end

"""
    write_checkpoint_data!(sim, tstep)

Writes model's floe, ocean, and atmosphere data to JLD2 file. Data can be used
to restart simulation run. Writes floes, ocean, and atmosphere to JLD2 file with name writer.fn for
current timestep, which will be the group in the JLD2 file. 

## _Positional arguments_
- $SIM_DEF
- $TSTEP_DEF
"""
function write_checkpoint_data!(sim, tstep)
    for filepath in sim.writers.checkpointwriters.filepath[
        mod.(tstep, sim.writers.checkpointwriters.Δtout) .== 0
    ]
        jldopen(filepath, "a+") do file
            file[string("floes/", tstep)] = sim.model.floes
            file[string("ocean/", tstep)] = sim.model.ocean
            file[string("atmos/", tstep)] = sim.model.atmos
        end
    end
    return
end

"""
    write_floe_data!(sim, tstep)

Writes desired `FloeOutputWriter` data to JLD2 file. Writes desired fields writer.outputs
to JLD2 file with name writer.fn for current timestep, which will be the group in the JLD2 file. 

## _Positional arguments_
- `writers::StructArray{FloeWriter}`: list of floe writers 
- $FLOES_DEF
- $TSTEP_DEF
"""
function write_floe_data!(writers, floes, tstep)
    for w in writers[
        mod.(tstep, writers.Δtout) .== 0
    ]
        file = jldopen(w.filepath, "a+")
        for output in w.outputs
            vals = StructArrays.component(floes, output)
            write_field(file, string(output, "/", tstep), vals)
        end
        close(file)
    end
    # once parents are recorded, they should be deleted
    empty!.(floes.parent_ids) 
    return
end

write_field(file, group, vals) = write(file, group, vals)

"""
    write_grid_data!(sim, tstep)

Writes desired GridOutputWriter data to NetCDF file. Writes desired fields writer.outputs
to file with name writer.fn for current timestep.

## _Positional arguments_
- `writers::StructArray{GridWriter}`: list of grid writers 
- $FLOES_DEF
- $TOPO_FIELD
- $TSTEP_DEF
"""
function write_grid_data!(writers, floes, topography, tstep)
    w_idx = mod.(tstep, writers.Δtout) .== 0
    if !isempty(floes) && sum(w_idx) != 0
        @views for w in writers[w_idx]
            calc_eulerian_data!(floes, topography, w)
            istep = div(tstep, w.Δtout) + 1  # Julia indicies start at 1
            
            # Open file and write data from grid writer
            ds = NCDataset(w.filepath, "a")
            for i in eachindex(w.outputs)
                name = string(w.outputs[i])
                ds[name][istep, :, :] = w.data[:, :, i]
            end
            close(ds)
        end
    end
    return
end

#----------------------- File Setup -----------------------#

#=
Initializes a JLD2 file in the given directory with the given filename. Setup
file to write given outputs. Each output is a group within the file.
=#
function _initialize_jld2_file!(dir, filename, overwrite, outputs, jld2_kw)
    # create path to file
    mkpath(dir)
    filename = auto_extension(filename, ".jld2")
    filepath = joinpath(dir, filename)
    # check if file already exists
    if !overwrite && isfile(filepath)
        throw(ErrorException("File $filepath already exists and overwrite is \
            false."))
    end
    overwrite && isfile(filepath) && rm(filepath, force=true)
    try
        # Create the file and the needed groups
        jldopen(filepath, "a+"; jld2_kw...) do file
            # Add groups for each output and each output's metadata
            for o in outputs
                JLD2.Group(file, string(o))
                # Add metadata
                attrs = getattrs(o)
                file[string("metadata/", o)] = attrs
            end
        end
    catch err
        @warn """Initialization of $filepath failed because of \
            $(sprint(showerror, err))"""
    end
    return filepath
end

#=
Initializes a NetCDF file in the given directory with the given filename. Setup
file to write given outputs. Create NetCDF file dir/filename with each output
added as a variable and with the dimensions time, x, and y. 
=#
function _initialize_netcdf_file!(
    ::Type{FT},
    dir,
    filename,
    overwrite,
    outputs,
    xg,
    yg,
) where {FT}
    mkpath(dir)
    filename = auto_extension(filename, ".nc")
    filepath = joinpath(dir, filename)
    if !overwrite && isfile(filepath)
        throw(ErrorException("File $filepath already exists and overwrite is \
            false."))
    end
    overwrite && isfile(filepath) && rm(filepath, force=true)
    try
        # Create the file and the needed groups
        dataset = NCDataset(filepath, "c")
        dataset.attrib["type"] = "Floe data averaged on the grid. The grid is \
            broken down into user provided dimensions."

        # Define dimensions time, x, and y
        defDim(dataset, "time", Inf)
        defVar(
            dataset,
            "time",
            FT,
            ("time",),
            attrib = Dict("units" => "10 seconds"),
        )

        defDim(dataset, "x", length(xg) - 1)
        x = defVar(
            dataset,
            "x",
            FT,
            ("x",),
            attrib = Dict("units" => "meters"),
        )
        x[:] = xg[1:end-1] .+ 0.5(xg[2] - xg[1])

        defDim(dataset, "y", length(yg) - 1)
        y = defVar(
            dataset,
            "y",
            FT,
            ("y",),
            attrib = Dict("units" => "meters"),
        )
        y[:] = yg[1:end-1] .+ 0.5(yg[2] - yg[1])

        # Define variables
        for o in outputs
            unit, comment = getattrs(o)
            defVar(
                dataset,
                string(o),
                FT,
                ("time", "x", "y"),
                attrib = Dict("units" => unit, "comments" => comment),
            )
        end
        # Write to file and close
        close(dataset)
    catch err
        @warn """Initialization of $filepath failed because of \
            $(sprint(showerror, err))"""
    end
    return filepath
end

#----------------------- Utils -----------------------#

# If `filename` ends in `ext`, return `filename`. Otherwise return `filename * ext`.
function auto_extension(filename, ext) 
    Next = length(ext)
    filename[end-Next+1:end] == ext || (filename *= ext)
    return filename
end

#=
Creates x-grid and y-grid. Assume xlines has length n and ylines has length m.
xgrid is the grid's xline vector repeated m times as rows in a mxn array and
ygrid is the yline vector repeated n times as columns in a mxn vector. xlines
and ylines are typically either xg and yg or xc and yc.
=#
function grids_from_lines(xlines, ylines)
    xgrid = repeat(reshape(xlines, 1, :), inner=(length(ylines),1))
    ygrid = repeat(ylines, outer = (1, length(xlines)))
    return xgrid, ygrid
end

"""
    calc_eulerian_data!(floes, topography, writer, istep)

Calculate floe data averaged on grid defined by `GridOutputWriter` for current
timestep (istep). Saved in writer.data field 

## _Positional arguments_
- $FLOES_DEF
- $TOPO_FIELD
- `writer::GridWriter`: grid output writer
"""
function calc_eulerian_data!(floes::FLT, topography, writer) where {FT <: AbstractFloat, FLT <: StructArray{<:Floe{FT}}}
    # Calculate/collect needed values
    Δx = writer.xg[2] - writer.xg[1]
    Δy = writer.yg[2] - writer.yg[1]
    cell_rmax = sqrt(Δx^2 + Δy^2)
    xgrid, ygrid = grids_from_lines(
        writer.xg[1:end-1] .+ 0.5Δx,
        writer.yg[1:end-1] .+ 0.5Δy,
    )
    dims = size(xgrid)
    floe_centroids = floes.centroid
    floe_rmax = floes.rmax

    # Identify floes that potentially overlap each grid square and create mask
    _potential_interactions = zeros(dims[1], dims[2], length(floes))
    for i in eachindex(floes)
        pint = sqrt.(
            (xgrid .- floe_centroids[i][1]).^2 .+
            (ygrid .- floe_centroids[i][2]).^2,
        )
        pint = pint .- (floe_rmax[i] + cell_rmax)
        pint[pint .> 0] .= 0
        pint[pint .< 0] .= 1
        _potential_interactions[:,:,i] = pint
    end
    FT_area_poly(p) = area_poly(p, FT)
    # Loop over each grid square
    for j in 1:dims[2]
        for i in 1:dims[1]
            pint = _potential_interactions[i,j,:]
            # If there are any potential interactions
            if sum(pint) > 0
                cell_poly_list = [_make_bounding_box_polygon(FT, writer.xg[j], writer.xg[j+1], writer.yg[i], writer.yg[i+1])]
                if length(topography) > 0
                    cell_poly_list = diff_polys(make_multipolygon(cell_poly_list), make_multipolygon(topography.poly), FT)
                end
                
                if length(cell_poly_list) == 0
                    writer.data[j, i, :] .= 0.0
                    continue
                end

                floeidx = collect(1:length(floes))[pint .== 1]
                #=  fic -> floes in cell - entirety of floe that is partially
                        within grid cell
                    pic -> partially in cell - only includes pieces of floes
                        that are within grid bounds
                =#
                pic_area = zeros(length(floeidx))
                for (i, idx) in enumerate(floeidx)
                    floe_poly = floes.poly[idx]
                    pic_area[i] = mapreduce(x -> sum(FT_area_poly, Subzero.intersect_polys(floe_poly, x, FT); init = 0.0), +, cell_poly_list; init = 0.0)
                end
                
                floeidx = floeidx[pic_area .> 0]
                pic_area = pic_area[pic_area .> 0]
                fic = floes[floeidx]
                fic_area = floes.area[floeidx]
                fic_mass = floes.mass[floeidx]

                area_ratios = pic_area ./ fic_area
                area_tot = sum(pic_area)
                mass_tot = sum(fic_mass .* area_ratios)

                if mass_tot > 0
                    # mass and area ratios
                    ma_ratios = area_ratios .* (fic_mass ./ mass_tot)
                    outputs = writer.outputs
                    for k in eachindex(outputs)
                        data = if outputs[k] == :u_grid
                            sum(floes.u[floeidx] .* ma_ratios)
                        elseif outputs[k] == :v_grid
                            sum(floes.v[floeidx] .* ma_ratios)
                        elseif outputs[k] == :dudt_grid
                            sum(floes.p_dudt[floeidx] .* ma_ratios)
                        elseif outputs[k] == :dvdt_grid
                            sum(floes.p_dvdt[floeidx] .* ma_ratios)
                        elseif outputs[k] == :si_frac_grid
                            area_tot/sum(FT_area_poly, cell_poly_list; init = 0.0)
                        elseif outputs[k] == :overarea_grid
                            sum(floes.overarea[floeidx])/length(floeidx)
                        elseif outputs[k] == :mass_grid
                            mass_tot
                        elseif outputs[k] == :area_grid
                            area_tot
                        elseif outputs[k] == :height_grid
                            sum(floes.height[floeidx] .* ma_ratios)
                        elseif outputs[k] == :stress_xx_grid
                            sum([s[1, 1] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_yx_grid
                            sum([s[1, 2] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_xy_grid
                            sum([s[2, 1] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_yy_grid
                            sum([s[2, 2] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_eig_grid
                            xx = sum([s[1, 1] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                            yx = sum([s[1, 2] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                            xy = sum([s[2, 1] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                            yy = sum([s[2, 2] for s in floes.stress_accum[floeidx]] .* ma_ratios)
                            stress = maximum(eigvals([xx yx; xy yy]))
                            if abs(stress) > 1e8
                                stress = 0.0
                            end
                            stress
                        elseif outputs[k] == :strain_ux_grid
                            sum([s[1, 1] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_vx_grid
                            sum([s[1, 2] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_uy_grid
                            sum([s[2, 1] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_vy_grid
                            sum([s[2, 2] for s in fic.strain] .* ma_ratios)
                        # Write bond_stress data too
                        ##-------------------------##  
                        elseif outputs[k] == :stress_bond_xx_grid
                            sum([s[1, 1] for s in floes.stress_bond[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_bond_yx_grid
                            sum([s[1, 2] for s in floes.stress_bond[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_bond_xy_grid
                            sum([s[2, 1] for s in floes.stress_bond[floeidx]] .* ma_ratios)
                        elseif outputs[k] == :stress_bond_yy_grid
                            sum([s[2, 2] for s in floes.stress_bond[floeidx]] .* ma_ratios)
                        ##-------------------------##    
                        end
                        writer.data[j, i, k] = data
                    end
                else
                    writer.data[j, i, :] .= 0.0
                end
            else
                writer.data[j, i, :] .= 0.0
            end
        end
    end
    return
end

#----------------------- Metadata -----------------------#

#=
Returns unit and comment attributes for each output type to be saved within
output NetCDF file
=#
function getattrs(output::Symbol)
    unit, comment =
        output == :sim ? ("Struct", "Starting state of simulation after initialization") :
        # Model structs/struct arrays
        output == :ocean ? ("Struct", "Model's ocean struct at given simulation timesteps") :
        output == :atmos ? ("Struct", "Model's atmosphere struct at given simulation timesteps") :
        output == :floes ? ("StructArray", "Model's array of floe structs at given simulation timesteps") :
        # Floe fields
        output == :centroid ? ("location", "Coordinates of centorid as [x y]") :
        output == :height ? ("m", "Floe height (uniform over floe)") :
        output == :area ? ("m^2", "Floe area") :
        output == :mass ? ("kg", "Floe mass") :
        output == :rmax ? ("m", "Maximum radius within floe") :
        output == :moment ? ("kg m^2", "Floe mass moment of intertia") :
        output == :angles ? ("kg m^2", "Floe interior angles") :
        output == :hex_bonds ? ("Floe Bonds - Neighbpuring floe ids connected with a bond") :
        output == :x_subfloe_points ? ("location", "x-coordinates for floe monte carlo points") :
        output == :y_subfloe_points ? ("location", "y-coordinates for floe monte carlo points") :
        output == :α ? ("radians", "Floe rotation since starting position") :
        output == :u ? ("m/s", "Floe x-direction velocity") :
        output == :v ? ("m/s", "Floe y-direction velocity") :
        output == :ξ ? ("rad/s", "Floe angular velocity") :
        output == :status ? ("unitless", "Flag if floe is still active in simulation") :
        output == :id ? ("unitless", "Floe ID - represents unique floe (ghosts have the same as parents since they are copies)") :
        output == :ghost_id ? ("unitless", "Ghost ID - represents individual ghosts so every ghost of 1 parent has unique ID (parents have ghost ID of 0)") :
        output == :ghosts ? ("unitless", "Index of floe's ghost floes in list") :
        output == :fxOA ? ("N", "X-directional forces on floe from the ocean and atmosphere") :
        output == :fyOA ? ("N", "Y-directional forces on floe from the ocean and atmosphere") :
        output == :trqOA ? ("N m", "Torque on floes from ocean and atmosphere") :
        output == :hflx ? ("W/m^2", "Average heatflux acting on the floe") : 
        output == :overarea ? ("m^2", "Overlap area of floe with other floes") :
        output == :collision_force ? ("N", "Total collision force from floes/topography/border in form [xf yf]") :
        output == :collision_torque ? ("N", "Total collision torque from floes/topography/border") :
        output == :bond_force ? ("N", "Total bond force") :
        output == :bond_torque ? ("N", "Total bond torque") :
        output == :stress_bond ? ("N/m^2", "Total bond stress") :
        output == :interactions ? ("[unitless, N, N, location, location, N m, m^2]",
                                    "Matrix of floe's interactions with following columns: ID of floe interacted with,
                                    collision x-force on floe, collision y-force on floe, collision x-point, collision y-point,
                                    collision torque on floe, overlap with other floe") :
        output == :stress_accum ? ("N/m^2", "Stress on the floe in the form [xx yx; xy, yy]") :
        output == :stress_instant ? ("N/m^2", "List of last stresses felt on the floe from previous timesteps") :
        output == :strain ? ("unitless", "Strain on the floe in the form [ux vx; uy vy]") :
        output == :p_dxdt ? ("m/s", "Floe x-velocity from previous time step") :
        output == :p_dydt ? ("m/s", "Floe y-velocity from previous time step") :
        output == :p_dudt ? ("m/s^2", "Floe x-acceleration from previous time step") :
        output == :p_dvdt ? ("m/s^2", "Floe y-acceleration from previous time step") :
        output == :p_dξdt ? ("rad/s^2", "Floe angular acceleration from the previous timestep") :
        output == :p_dαdt ? ("rad/s", "Floe angular velocity from previous timestep") :
        # Gridded data fields
        output == :u_grid ? ("m/s", "Average floe x-direction velocity in grid cell") :
        output == :v_grid ? ("m/s", "Average floe y-direction velocity in grid cell") :
        output == :dudt_grid ? ("m/s^2", "Average floe x-direction acceleration in grid cell") :
        output == :dvdt_grid ? ("m/s^2", "Average floe y-direction acceleration in grid cell") :
        output == :overarea_grid ? ("m", "Average overlap area of floe with other floes in grid cell") :
        output == :mass_grid ? ("kg", "Average mass of floes in grid cell") :
        output == :area_grid ? ("m^2", "Average area of floes in grid cell") :
        output == :height_grid ? ("m", "Average height of floes in grid cell") :
        output == :si_frac_grid ? ("unitless", "Fraction of grid cell covered by floes") :
        output == :stress_xx_grid ? ("N/m^2", "Average xx stress on floes in a given grid cell") :
        output == :stress_yx_grid ? ("N/m^2", "Average yx stress on floes in a given grid cell") :
        output == :stress_xy_grid ? ("N/m^2", "Average xy stress on floes in a given grid cell") :
        output == :stress_yy_grid ? ("N/m^2", "Average yy stress on floes in a given grid cell") :
        output == :stress_eig_grid ? ("N/m^2", "Maximum eigenvalue of the stress matricies [xx yx; xy yy]") :
        output == :strain_ux_grid ? ("unitless", "Average ux strain on floes in a grid cell") :
        output == :strain_vx_grid ? ("unitless", "Average vx strain on floes in a grid cell") :
        output == :strain_uy_grid ? ("unitless", "Average uy strain on floes in a grid cell") :
        output == :strain_vy_grid ? ("unitless", "Average vy strain on floes in a grid cell") :
        # Adding bond stresses too
        output == :stress_bond_xx_grid ? ("N/m^2", "Average xx stress on floes in a given grid cell for bonds") :
        output == :stress_bond_yx_grid ? ("N/m^2", "Average yx stress on floes in a given grid cell for bonds") :
        output == :stress_bond_xy_grid ? ("N/m^2", "Average xy stress on floes in a given grid cell for bonds") :
        output == :stress_bond_yy_grid ? ("N/m^2", "Average yy stress on floes in a given grid cell for bonds") :
        ("", "") # if symbol isn't found, return empty attributes
    return unit, comment
end
