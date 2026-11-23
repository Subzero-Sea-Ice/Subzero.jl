# Whether the vectors (Δx1, Δy1) and (Δx2, Δy2) are parallel, exactly, as GeometryOps'
# `_isparallel(Δx1, Δy1, Δx2, Δy2)`
_isparallel(Δx1, Δy1, Δx2, Δy2) = (Δx1 * Δy2 == Δy1 * Δx2)
