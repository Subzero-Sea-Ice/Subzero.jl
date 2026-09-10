module Subzero

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Subzero

import Base.@kwdef # this is being exported as of version 1.9
import Base.show
import GeometryOps as GO
import GeometryOps.GeoInterface as GI
import StaticArrays as SA
using CoordinateTransformations, Dates, Extents, 
    Interpolations, JLD2, LinearAlgebra, Logging, Measures, NCDatasets,
    Printf, Random, Rotations, SplitApplyCombine, Statistics, StructArrays
using Adapt, CUDA, KernelAbstractions


const SIM_DEF = "`sim::Simulation`: simulation to be run"
const FT_DEF = "`FT::Type{<:AbstractFloat}`: Float type used to run the simulation, either \
`Float64` (default) or `Float32`."
const POLY_DEF = "`poly::Polys{FT}`: Polygon used to represent the shape of a floe or topography"
const POLY_LIST_DEF = "`polys::Vector{<:Polygon}`: list of polygons meant to represent a field \
of floes or topography elements. Polygons can be any polygon type that supports GeoInterface."
const COORDS_LIST_DEF = "`coords::Vector{<:PolyVec}`: list of polygon coordinates meant to \
represent a field of floes or topography elements. PolyVec refers to a Vector{Vector{<:Points}} \
where the points can be tuples, vectors, or static vectors and the innermost vector refers to each ring of the polygon."
const CENTROID_DEF = "`centroid::Vector{FT}`: Two-element vector meant to represent the (x, y) \
point that is the centroid of either a floe or topography"
const RMAX_DEF = "`rmax::FT`: Float length representing the maximum radius of a floe or topography \
from the centroid to any given vertex"


# Types
include("tools/geom_utils.jl")
# Model
include("simulation_components/grids.jl")
include("simulation_components/domain_components/abstract_domains.jl")
include("simulation_components/domain_components/boundaries.jl")
include("simulation_components/domain_components/topography.jl")
include("simulation_components/domain_components/domains.jl")
include("simulation_components/floe_components/floe_status.jl")
include("simulation_components/floe_components/floe_interaction.jl")
include("simulation_components/floe_components/floe.jl")
include("simulation_components/floe_components/floe_utils.jl")
include("simulation_components/floe_components/floe_field.jl")
include("simulation_components/floe_components/stress_calculators.jl")
include("simulation_components/floe_components/subfloe_points_generators.jl")
include("simulation_components/floe_components/floe_settings.jl")
include("simulation_components/oceans.jl")
include("simulation_components/atmos.jl")
include("simulation_components/model.jl")
# Outputs
include("simulation_components/output_components/logger.jl")
# Process Settings
include("simulation_components/process_settings/process_settings.jl")
include("simulation_components/process_settings/fracture_settings.jl")
# Simulation
include("simulation_components/output_components/output.jl")
include("simulation_components/constants.jl")
include("simulation_components/simulation.jl")
# Physical Processes
include("physical_processes/fractures.jl")
include("physical_processes/update_floe.jl")
include("physical_processes/coupling.jl")
include("physical_processes/collisions.jl")
include("physical_processes/simplification.jl")
include("physical_processes/ridge_raft.jl")
include("physical_processes/welding.jl")
# Tools
include("tools/plotting.jl")
include("tools/conservation_em.jl")
include("tools/compare_files.jl")
end