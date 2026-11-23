#=
Geometric predicates, as GeometryOps' `Predicates` module. `exact` is `True()` or `False()`.
=#
module Predicates

import AdaptivePredicates
import ..GeometryOpsGPU: True, False

#= Orientation of c relative to the line from a to b, as `GO.Predicates.orient(a, b, c;
exact)`: positive if c is to the left, negative if it's to the right, and 0 if c is on the
line. With `exact = True()`, AdaptivePredicates' `orient2p` decides exactly. With
`exact = False()`, the cross product is computed directly, and 0 if it's approximately zero.

TODO: like GeometryOps, both convert the points to Float64, also for Float32 points. =#
orient(a, b, c; exact) = _orient(exact, a, b, c)

_orient(::True, a, b, c) = AdaptivePredicates.orient2p(
    (Float64(a[1]), Float64(a[2])),
    (Float64(b[1]), Float64(b[2])),
    (Float64(c[1]), Float64(c[2])),
)

function _orient(::False, a, b, c)
    ax, ay = Float64(a[1]) - Float64(c[1]), Float64(a[2]) - Float64(c[2])
    bx, by = Float64(b[1]) - Float64(c[1]), Float64(b[2]) - Float64(c[2])
    c_t1 = ax * by
    c_t2 = ay * bx
    return isapprox(c_t1, c_t2) ? 0.0 : sign(c_t1 - c_t2)
end

end

import .Predicates
