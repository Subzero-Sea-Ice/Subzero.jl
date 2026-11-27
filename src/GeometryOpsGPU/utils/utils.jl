# Whether the point is in the closed box `ext`, as `GO._point_in_extent`
function _point_in_extent(p, ext::Extents.Extent)
    (x1, x2), (y1, y2) = ext.X, ext.Y
    return x1 ≤ p[1] ≤ x2 && y1 ≤ p[2] ≤ y2
end
