#=
Geometric predicates, as GeometryOps' `Predicates` module. `exact` is `True()` or `False()`.
=#
module Predicates

import AdaptivePredicates
import ..GeometryOpsGPU: True, False

#= Orientation of c relative to the line from a to b, as `GO.Predicates.orient(a, b, c;
exact)`: positive if c is to the left, negative if it's to the right, and 0 if c is on the
line. With `exact = True()`, AdaptivePredicates' `orient2p` decides exactly. With
`exact = False()`, the cross product is computed directly, and is 0 if it's approximately
zero (`isapprox`). Unlike GeometryOps, which converts the points to Float64, it computes in
the points' type, so Float32 points use Float32 maths. The points must all have the same
type. =#
orient(a::NTuple{2, T}, b::NTuple{2, T}, c::NTuple{2, T}; exact) where T =
    _orient(exact, a, b, c)

_orient(::True, a, b, c) = AdaptivePredicates.orient2p(a, b, c)

function _orient(::False, a, b, c)
    ax, ay = a[1] - c[1], a[2] - c[2]
    bx, by = b[1] - c[1], b[2] - c[2]
    c_t1 = ax * by
    c_t2 = ay * bx
    return isapprox(c_t1, c_t2) ? zero(c_t1) : sign(c_t1 - c_t2)
end

end

import .Predicates
