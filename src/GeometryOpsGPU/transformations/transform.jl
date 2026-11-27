"""
    transform!(f, poly::ArrayPolygon)

Apply `f` to every point of the polygon, in place, as `GO.transform(f, poly)`. `f` takes and
returns a point as a tuple `(x, y)`. Returns `nothing`.

To run in a kernel, `f` must be isbits and must not allocate, for example a
[`Translate`](@ref GeometryOpsGPU.Translate) or a
[`RotateTranslate`](@ref GeometryOpsGPU.RotateTranslate). They give the same results as
the CoordinateTransformations types that GeometryOps is usually given (`Translation`,
`LinearMap`, ...), but `RotateTranslate` computes the sine and cosine once, not for every
point.
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
`Δα` radians anticlockwise around `(cx, cy)`, and then moves it by `Δx` and `Δy`: the same,
bit for bit, as `Translation(Δx, Δy) ∘ recenter(LinearMap(Angle2d(Δα)), (cx, cy))` in
CoordinateTransformations and Rotations. All arguments must have the same type, that of the
polygon's coordinates. `sincos(Δα)` and the translation are computed once, here.
"""
struct RotateTranslate{T}
    sinα::T
    cosα::T
    tx::T
    ty::T
end

#=
CoordinateTransformations composes the transformation into one affine map, `p ↦ R*p + t`
with `t = Δ + (c + R*(-c))`, where R is the rotation and c the centre. This computes t, and
then each point, in the same order, so that the results are the same.
=#
function RotateTranslate(Δα, cx, cy, Δx, Δy)
    sinα, cosα = sincos(Δα)
    rx, ry = _rotate(sinα, cosα, -cx, -cy)
    return RotateTranslate(sinα, cosα, Δx + (cx + rx), Δy + (cy + ry))
end

@inline function (t::RotateTranslate)((x, y))
    rx, ry = _rotate(t.sinα, t.cosα, x, y)
    return (rx + t.tx, ry + t.ty)
end

# (x, y) rotated by the angle with sine s and cosine c, as Rotations' `Angle2d(α) * [x, y]`
@inline _rotate(s, c, x, y) = (c * x - s * y, s * x + c * y)
