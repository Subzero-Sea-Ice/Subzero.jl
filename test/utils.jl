#=
    find_poly_coords(poly)

Syntactic sugar for to find a polygon's coordinates
Input:
    poly    <Polygon>
Output:
    <PolyVec> representing the floe's coordinates xy plane
=#
find_poly_coords(poly) = GI.coordinates(poly)

# Make a copy of given coordinates and translate by given deltas, returning new coordiantes
function translate_coords(coords, Δx, Δy)
    new_coords = [[Vector{Float64}(undef, 2) for _ in eachindex(coords[1])]]
    for i in eachindex(coords[1])
        new_coords[1][i][1] = coords[1][i][1] + Δx
        new_coords[1][i][2] = coords[1][i][2] + Δy 
    end
    return new_coords
end

# Translate each of the given coodinates by given deltas in place
function translate_coords!(coords, Δx, Δy)
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            coords[i][j][1] += Δx
            coords[i][j][2] += Δy
        end
    end
    return
end