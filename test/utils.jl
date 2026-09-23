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

#=
Deterministic set of floes for testing `timestep_floe_properties!`. Floes have different
numbers of vertices and non-zero velocities, forces, and interactions. Some floes trigger
the height limit, collision force reduction, velocity adjustment, and ξ limit.
=#
function _make_timestep_test_floes(::Type{FT}, floe_settings; n = 10) where FT
    rng = Xoshiro(42)
    floes = StructArray{Floe{FT}}(undef, 0)
    for k in 1:n
        # convex polygon with k + 3 vertices around a random center
        nverts = k + 3
        cx, cy = 1e4 .* rand(rng, 2)
        r = 500 + 1500rand(rng)
        θs = sort(2π .* rand(rng, nverts))
        ring = [[cx + r * cos(θ), cy + r * sin(θ)] for θ in θs]
        push!(ring, ring[1])
        height = 0.5 + rand(rng)
        floe = Floe{FT}([ring], height; floe_settings, rng)
        floe.u, floe.v = 0.2 .* (rand(rng, 2) .- 0.5)
        floe.ξ = 1e-6 * (rand(rng) - 0.5)
        floe.p_dxdt, floe.p_dydt = 0.2 .* (rand(rng, 2) .- 0.5)
        floe.p_dαdt = 1e-6 * (rand(rng) - 0.5)
        floe.p_dudt, floe.p_dvdt = 1e-4 .* (rand(rng, 2) .- 0.5)
        floe.p_dξdt = 1e-9 * (rand(rng) - 0.5)
        floe.fxOA, floe.fyOA = 1e6 .* (rand(rng, 2) .- 0.5)
        floe.trqOA = 1e8 * (rand(rng) - 0.5)
        floe.hflx_factor = 1e-3 * (rand(rng) - 0.5)
        floe.collision_force = 1e6 .* (rand(rng, 1, 2) .- 0.5)
        floe.collision_trq = 1e8 * (rand(rng) - 0.5)
        floe.stress_accum = 1e3 .* (rand(rng, 2, 2) .- 0.5)
        floe.stress_instant = 1e3 .* (rand(rng, 2, 2) .- 0.5)
        floe.strain = 1e-6 .* (rand(rng, 2, 2) .- 0.5)
        floe.damage = rand(rng)
        if isodd(k)  # interactions, with unused (padding) rows for some floes
            ninters = k ÷ 2 + 1
            inters = zeros(FT, ninters + (k % 3), 7)
            inters .= 1e3 .* rand(rng, size(inters)...)
            inters[:, xpoint] .+= cx
            inters[:, ypoint] .+= cy
            floe.interactions = inters
            floe.num_inters = ninters
        end
        push!(floes, floe)
    end
    # edge cases
    floes.height[1] = 2floe_settings.max_floe_height  # height limit
    floes.collision_force[2] .= [1e3 -1e3] .* floes.mass[2]  # force reduction
    floes.fxOA[3] = 10 * floes.mass[3]  # velocity adjustment (frac != 1)
    floes.fyOA[4] = -10 * floes.mass[4]
    floes.fxOA[4] = 10 * floes.mass[4]
    floes.trqOA[5] = 1e3 * floes.moment[5]  # ξ limit
    return floes
end
