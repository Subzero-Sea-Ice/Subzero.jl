module Subzero

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Subzero

export
    Floe,
    Constants,
    Model,
    Simulation,
    timestep_sim!,
    run!,
    AbstractOutputWriter,
    CheckpointOutputWriter,
    GridOutputWriter, 
    FloeOutputWriter,
    InitialStateOutputWriter,
    GridOutput,
    FloeOutput,
    add_ghosts!,
    # Interaction field enum and elements
    InteractionFields,
    floeidx,
    xforce,
    yforce,
    xpoint,
    ypoint,
    torque,
    overlap, 
    initialize_floe_field,
    NoFracture,
    HiblerYieldCurve,
    MohrsCone,
    CollisionSettings,
    FractureSettings,
    CouplingSettings,
    SimplificationSettings,
    RidgeRaftSettings,
    WeldSettings,
    FloeSettings,
    PolyVec,
    OutputWriters,
    check_energy_momentum_conservation_julia,
    MonteCarloPointsGenerator,
    SubGridPointsGenerator

import Base.@kwdef # this is being exported as of version 1.9
import Base.show
import GeometryOps as GO
import GeometryOps.GeoInterface as GI
import GeometryOps.GeometryBasics as GB
import StaticArrays as SA
using CoordinateTransformations, Dates, Extents,
    Interpolations, JLD2, LinearAlgebra, Logging, Measures, NCDatasets,
    Printf, Random, Rotations, SplitApplyCombine, Statistics, StructArrays,
    VoronoiCells

"""
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates
 are the exterior border of the floe and the wz coordinates, or any other
 following sets of coordinates, describe holes within the floe.
 This form is for easy conversion to polygons.
"""
const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:Real

"""
Coordinates are vector of vector of points of the form:
[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]] where the xy coordinates form a
closed ring. PolyVec objects can be made out RingVec objects.
This form is for each conversion to LinearRings, which can also be made into Polygons.
"""
const RingVec{T} = R where {
    T<:Real,
    V<:AbstractArray{T},
    R <: AbstractArray{V},
}

const Polys{T, V} = GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, Nothing} where T

Base.convert(::Type{Polys{Float32}}, p::Polys{<:Real}) = GO.tuples(p, Float32)
Base.convert(::Type{Polys{Float64}}, p::Polys{<:Real}) = GO.tuples(p, Float64)

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
include("simulation_components/stress_calculators.jl")
# Model
include("simulation_components/grids.jl")
include("simulation_components/domain_components/abstract_domains.jl")
include("simulation_components/domain_components/boundaries.jl")
include("simulation_components/domain_components/topography.jl")
include("simulation_components/domain_components/domains.jl")
include("simulation_components/floe.jl")
include("floe_utils.jl")
include("simulation_components/oceans.jl")
include("simulation_components/atmos.jl")
include("simulation_components/model.jl")
# Physical Processes
include("physical_processes/fractures.jl")
include("physical_processes/update_floe.jl")
include("physical_processes/coupling.jl")
include("physical_processes/collisions.jl")
include("physical_processes/process_settings.jl")
include("physical_processes/simplification.jl")
include("physical_processes/ridge_raft.jl")
include("physical_processes/welding.jl")
# Tools
include("plotting.jl")
include("tools/conservation_em.jl")
include("tools/compare_files.jl")
include("logger.jl")
# Simulation
include("output.jl")
include("simulation_components/simulation.jl")
end