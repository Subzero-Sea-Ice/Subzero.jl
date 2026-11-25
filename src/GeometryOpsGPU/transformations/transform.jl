"""
    transform!(f, poly::ArrayPolygon)

Apply `f` to every point of the polygon, in place, as `GO.transform(f, poly)`. `f` takes and
returns a point as a tuple `(x, y)`. Returns `nothing`.

To run in a kernel, `f` must be isbits and must not allocate, for example a
[`Translate`](@ref GeometryOpsGPU.Translate) or a
[`RotateTranslate`](@ref GeometryOpsGPU.RotateTranslate). GeometryOps is usually given
CoordinateTransformations' types (`Translation`, `LinearMap`, ...) instead, but they
compute in another order, so they wouldn't give the same results as the code that this
replaces.
"""
@inline function transform!(f, p::ArrayPolygon)
    for j in 1:npoint(p)
        _setpoint!(p, j, f(getpoint(p, j)))
    end
    return
end

"""
    Translate(Δx, Δy)

A point function for [`transform!`](@ref GeometryOpsGPU.transform!) that moves a point by
`Δx` and `Δy`: `(x + Δx, y + Δy)`, as CoordinateTransformations' `Translation(Δx, Δy)`.
"""
struct Translate{T}
    Δx::T
    Δy::T
end

@inline (t::Translate)((x, y)) = (x + t.Δx, y + t.Δy)

"""
    RotateTranslate(Δα, cx, cy, Δx, Δy)

A point function for [`transform!`](@ref GeometryOpsGPU.transform!) that rotates a point by
`Δα` radians anticlockwise around `(cx, cy)`, and then moves it by `Δx` and `Δy`. That's
`Translation(Δx, Δy) ∘ recenter(LinearMap(Angle2d(Δα)), (cx, cy))` in
CoordinateTransformations and Rotations, up to round-off. All arguments must have the same
type, that of the polygon's coordinates. `sincos(Δα)` is computed once, here.
"""
struct RotateTranslate{T}
    sinα::T
    cosα::T
    cx::T
    cy::T
    Δx::T
    Δy::T
end

function RotateTranslate(Δα, cx, cy, Δx, Δy)
    sinα, cosα = sincos(Δα)
    return RotateTranslate(sinα, cosα, cx, cy, Δx, Δy)
end

@inline function (t::RotateTranslate)((x, y))
    x, y = x - t.cx, y - t.cy
    return (t.cosα * x - t.sinα * y + t.cx + t.Δx, t.sinα * x + t.cosα * y + t.cy + t.Δy)
end
