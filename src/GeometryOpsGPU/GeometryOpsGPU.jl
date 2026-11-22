"""
    GeometryOpsGPU

Geometry functions that run in GPU kernels (and on the CPU), for polygons stored in plain
arrays ([`ArrayPolygon`](@ref GeometryOpsGPU.ArrayPolygon)). They don't allocate and don't
build GeometryOps or GeoInterface objects.

The module mirrors [GeometryOps](https://github.com/JuliaGeo/GeometryOps.jl), so that
someone who knows GeometryOps finds the same functions here, and so that the code could be
moved into GeometryOps one day:
- Functions with an equivalent in GeometryOps have its name and argument order
  (`transform!(f, poly)` for `GO.transform(f, geom)`), with a `!` if they write into their
  arguments instead of returning a new geometry. Internal helpers have the names of the
  GeometryOps helpers they follow. Functions without an equivalent say so.
- Files are laid out as in GeometryOps' `src/` (`methods/`, `transformations/`, ...).
- With Float64 polygons, results are the same as GeometryOps', bit for bit. With Float32
  polygons, the maths is in Float32 throughout, even where GeometryOps converts to Float64,
  so results can differ by round-off.

The differences from GeometryOps:
- Polygons are `ArrayPolygon`s, not any GeoInterface geometry, and have no holes.
- Results are written into preallocated arrays rather than returned as new geometries.
- There's no manifold or algorithm argument: everything is planar.

This module doesn't use Subzero, so that it can become a package of its own. Subzero
imports it as `GOG`, as it imports GeometryOps as `GO`.
"""
module GeometryOpsGPU

import Extents

include("types.jl")
include("transformations/transform.jl")
include("methods/geom_relations/interiors_intersect.jl")

end
