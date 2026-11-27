#=
GeoInterface methods for `ArrayPolygon`, so that GeometryOps (`GO.area`, `GO.intersects`,
`GO.transform`, ...) and plotting accept it on the host. Kernels don't use them: they call
`getpoint` and `npoint` directly. They only work for arrays on the host, e.g. floes copied
back with `adapt(Array, floes)`, since GeoInterface reads the points one at a time.
=#

# The exterior ring of an ArrayPolygon: the same row of the same array
struct ArrayRing{A<:AbstractArray{<:Any, 3}}
    coords::A
    i::Int
    n::Int
end

GI.isgeometry(::Type{<:ArrayPolygon}) = true
GI.geomtrait(::ArrayPolygon) = GI.PolygonTrait()
GI.ngeom(::GI.PolygonTrait, ::ArrayPolygon) = 1  # no holes
GI.getgeom(::GI.PolygonTrait, p::ArrayPolygon, k::Integer) = ArrayRing(p.coords, p.i, p.n)
GI.ncoord(::GI.PolygonTrait, ::ArrayPolygon) = 2
GI.extent(::GI.PolygonTrait, p::ArrayPolygon) = extent(p)

GI.isgeometry(::Type{<:ArrayRing}) = true
GI.geomtrait(::ArrayRing) = GI.LinearRingTrait()
GI.ngeom(::GI.LinearRingTrait, r::ArrayRing) = r.n
GI.getgeom(::GI.LinearRingTrait, r::ArrayRing, j::Integer) =
    (r.coords[r.i, j, 1], r.coords[r.i, j, 2])  # tuples are GeoInterface points
GI.ncoord(::GI.LinearRingTrait, ::ArrayRing) = 2
