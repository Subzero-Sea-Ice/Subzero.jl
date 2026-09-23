#=
Functions to update floe's shape and fields. Function typically called from
other physical processes, other than timestep_floe_properties!, which is called
from the timestep_simulation! function.
=#

#=
Updates existing floe shape and related physical properties based of the polygon
defining the floe.

Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> floe to update
    new_poly    <Polygon> polygon representing new outline of floe
    new_mass    <AbstractFloat> mass of floe
    floe_settings   <FloeSettings> simulation's settings for making floes
    rng         <RNG> random number generator

Ouputs:
    Updates a given floe's physical properties given new shape and total mass.
=#
function replace_floe!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    new_poly,
    new_mass,
    floe_settings,
    rng,
) where {FT}
    # Floe shape
    floe.poly = new_poly
    floe.centroid = collect(centroid_poly(new_poly, FT))
    floe.area = area_poly(new_poly, FT)
    floe.height = new_mass/(floe.area * floe_settings.ρi)
    floe.mass = new_mass
    floe.moment = _calc_moment_inertia(
        FT,
        new_poly,
        floe.centroid,
        floe.height;
        ρi = floe_settings.ρi,
    )
    floe.angles = angles_poly(floe.poly, FT)
    floe.α = FT(0)
    floe.rmax = _calc_max_radius(floe.poly, floe.centroid, FT)
    # Floe monte carlo points
    x_subfloe_points, y_subfloe_points, status = generate_subfloe_points(
        floe_settings.subfloe_point_generator,
        floe.poly,
        floe.centroid,
        floe.area,
        floe.status,
        rng,
    )
    floe.x_subfloe_points = x_subfloe_points
    floe.y_subfloe_points = y_subfloe_points
    # Floe status / identification
    floe.status = status
    return
end

#=
Update current and previous velocity/acceleration fields to conserve momentum of
a floe whose shape has been changed, given the previous mass, momentum, and
centroid.

Inputs:
    mass_tmp    <AbstractFloat> original mass of floe before shape change
    moment_tmp  <AbstractFloat> original moment of intertia of floe before shape
                    change
    x_tmp       <AbstractFloat> original x-coordinate of centroid of floe before
                    shape change
    y_tmp       <AbstractFloat> original y-coordinate of centroid of floe before
                    shape change
    Δt          <Int> timestep of simulation in seconds
    keep_floe   <Union{Floe, LazyRow{Floe}}> floe whose shape has been changed
    combine_floe <Union{Floe, LazyRow{Floe}}> if keep_floe's shape has been
                    changed due to an interaction with another floe, combine_floe
                    is that floe - optional parameter

Output:
    None. keep_floe's u, v, ξ, p_dxdt, p_dydt, p_dαdt, p_dudt, p_dvdt, and
    p_dξdt fields all updated to preserve momentum. 

Note:
    Function does not depend on conservation of mass
=#
function conserve_momentum_change_floe_shape!(
    mass_tmp,
    moment_tmp,
    x_tmp,
    y_tmp,
    Δt,
    keep_floe,
    combine_floe = nothing,
)
    # Calculate linear velocities to conserve linear momentum
    new_u = keep_floe.u * mass_tmp
    new_v = keep_floe.v * mass_tmp
    new_dxdt = keep_floe.p_dxdt * mass_tmp
    new_dydt = keep_floe.p_dydt * mass_tmp
    if !isnothing(combine_floe)
        new_u += combine_floe.u * combine_floe.mass
        new_v += combine_floe.v * combine_floe.mass
        new_dxdt += combine_floe.p_dxdt * combine_floe.mass
        new_dydt += combine_floe.p_dydt * combine_floe.mass
    end
    new_u /= keep_floe.mass
    new_v /= keep_floe.mass
    new_dxdt /= keep_floe.mass
    new_dydt /= keep_floe.mass
    # Calculate angular velocities to conserve rotational angular momentum
    p_x = x_tmp - Δt * keep_floe.p_dxdt
    p_y = y_tmp - Δt * keep_floe.p_dydt
    new_ξ = keep_floe.ξ * moment_tmp +  # spin momentum + orbital momentum
        mass_tmp * (x_tmp * keep_floe.v - y_tmp * keep_floe.u)
    new_dαdt = keep_floe.p_dαdt * moment_tmp +
        mass_tmp * (p_x * keep_floe.p_dydt - p_y * keep_floe.p_dxdt)
    if !isnothing(combine_floe)
        p_x = combine_floe.centroid[1] - Δt * combine_floe.p_dxdt
        p_y = combine_floe.centroid[2] - Δt * combine_floe.p_dydt
        new_ξ += combine_floe.ξ * combine_floe.moment +
            combine_floe.mass * (
                combine_floe.centroid[1] * combine_floe.v -
                combine_floe.centroid[2] * combine_floe.u
            )
        new_dαdt += combine_floe.p_dαdt * combine_floe.moment +
            combine_floe.mass * (
                p_x * combine_floe.p_dydt -
                p_y * combine_floe.p_dxdt
            )
    end
    # Subtract new orbital velocity from total previous momentum
    p_x = keep_floe.centroid[1] - Δt * new_dxdt
    p_y = keep_floe.centroid[2] - Δt * new_dydt
    new_ξ -= keep_floe.mass * (
        keep_floe.centroid[1] * new_v -
        keep_floe.centroid[2] * new_u
    )
    new_dαdt -= keep_floe.mass * (p_x * new_dydt - p_y * new_dxdt)
    new_ξ /= keep_floe.moment
    new_dαdt /= keep_floe.moment

    # Set new values
    keep_floe.u = new_u
    keep_floe.v = new_v
    keep_floe.ξ = new_ξ
    keep_floe.p_dxdt = new_dxdt
    keep_floe.p_dydt = new_dydt
    keep_floe.p_dαdt = new_dαdt
    # Calculate previous accelerations
    keep_floe.p_dudt = (keep_floe.u - keep_floe.p_dxdt) / Δt
    keep_floe.p_dvdt = (keep_floe.v - keep_floe.p_dydt) / Δt
    keep_floe.p_dξdt = (keep_floe.ξ - keep_floe.p_dαdt) / Δt
    return
end

#=
Update rotational velocities of two floes whose shapes may have changed to
maintain conservation of momentum given the floes' intial roational momentums
and the change in orbital momentum between their original shape and new shapes.
Additionally, both floes will have the same velocity at a given (x,y) point.

Inputs:
    x                   <Float> x-coordinate of point where floes share same
                            rotational velocity
    y                   <Float> y-coordinate of point where floes share same
                            rotational velocity
    floe1               <Union{Floe, LazyRow{Floe}}> first floe in pair
    floe2               <Union{Floe, LazyRow{Floe}}> second floe in pair
    init_rot_momentum   <Float> initial rotational momentum of floe 1 before
                            shape change
    init_p_rot_momentum <Float> initial rotational momentum of floe 2 before
                            shape change
    diff_orbital        <Float> change in floe 1's orbital velocity after shape
                            change
    diff_p_orbital      <Float> change in floe 1's orbital velocity after shape
                            change
    Δt                  <Int> length of timestep in seconds

Output:
    Nothing. Update both floes' rotational velocity, previous rotational
    velocity, and rotational acceleration to conserve angular momentum.
=#
function update_new_rotation_conserve!(
    x, y,
    floe1, floe2,
    init_rot_momentum, init_p_rot_momentum,
    diff_orbital, diff_p_orbital,
    Δt,
)
    # Find radius of each polygon to shared midpoint
    rad1 = sqrt(
        (floe1.centroid[1] - x)^2 +
        (floe1.centroid[2] - y)^2
    )
    rad2 = sqrt(
        (floe2.centroid[1] - x)^2 +
        (floe2.centroid[2] - y)^2
    )
    rad_ratio = rad1 / rad2
    # Determine ξ values so they are stationary at midpoint of shared edges
    floe1.ξ = (diff_orbital + init_rot_momentum) /
        (floe1.moment - (floe2.moment * rad_ratio))
    floe2.ξ = -floe1.ξ * rad_ratio
    # Determine p_dαdt values so they are stationary at intersection point
    floe1.p_dαdt = (diff_p_orbital + init_p_rot_momentum) /
        (floe1.moment - (floe2.moment * rad_ratio))
    floe2.p_dαdt = -floe1.p_dαdt * rad_ratio
    # Calculate previous rotational accelerations
    floe1.p_dξdt = (floe1.ξ - floe1.p_dαdt) / Δt
    floe2.p_dξdt = (floe2.ξ - floe2.p_dαdt) / Δt
end

#=
Update new_floes's current and previous velocity/acceleration fields to conserve
momentum when a floe has been fractured into several new floes, given the
previous mass, momentum, and centroid. The assumption is made that each new floe
has the same velocities/accelerations.

Inputs:
    init_floe   <Union{Floe, LazyRow{Floe}}> original floe
    new_floes   <StructArray{Floe}> fractured pieces of original floe
    Δt          <Int> simulation's timestep in seconds

Output:
    None. new_floes velocities and accelerations are updated for current and
    previous timestep to conserve momentum.

Note: Depends on conservation of mass.
=#
function conserve_momentum_fracture_floe!(
    init_floe,
    new_floes::StructArray{<:Floe{FT}},
    Δt,
) where {FT}
    if !isempty(new_floes)
        # conserve linear momentum by keeping linear velocities the same
        new_floes.u .= init_floe.u 
        new_floes.v .= init_floe.v
        new_floes.p_dxdt .= init_floe.p_dxdt
        new_floes.p_dydt .= init_floe.p_dydt
        new_floes.p_dudt .= init_floe.p_dudt
        new_floes.p_dvdt .= init_floe.p_dvdt
        # Reset α with new coordinates
        new_floes.α .= 0
        # MATLAB uses assumptions -> doesn't conserve rotational momentum
        new_floes.ξ .= init_floe.ξ
        new_floes.p_dαdt .= 0
        new_floes.p_dξdt .= init_floe.p_dξdt
    end
end

#=
Conserve linear momentum when mass is transfered from one floe to another floe.
Inputs:
    floes       <StructArray{Floes}> list of floes
    idx1        <Int> index of first floe in floes list
    idx2        <Int> index of second floe in floes list
    m1          <Float> initial mass of first floe before mass transfer
    m2          <Float> initial mass of second floe before mass transfer
    Δt          <Int> timestep of simulation in seconds
    pieces_list <Union{nothing, Vector{Floe}}> list of floes created from floes
                    breaking during ridging and rafting
    pieces_idx  <Int> start index of pieces included in conservation
                    calculations within the pieces_list
Outputs:
    Nothing. Update floes' velocities and accelerations to conserve linear
    momentum.
=#
function conserve_momentum_transfer_mass!(
    floes,
    idx1, idx2,
    m1, m2,
    Δt,
    pieces_list = nothing,
    pieces_idx = 0,
)
    # Conserve linear - assume resultant floes have same linear velocity
    tot_mass = m1 + m2
    new_u = (m1 * floes.u[idx1] + m2 * floes.u[idx2]) / tot_mass
    new_v = (m1 * floes.v[idx1] + m2 * floes.v[idx2]) / tot_mass
    new_p_dxdt = (m1 * floes.p_dxdt[idx1] + m2 * floes.p_dxdt[idx2]) / tot_mass
    new_p_dydt = (m1 * floes.p_dydt[idx1] + m2 * floes.p_dydt[idx2]) / tot_mass
    new_p_dudt = (floes.u[idx1] - floes.p_dxdt[idx1]) / Δt
    new_p_dvdt = (floes.v[idx1] - floes.p_dydt[idx1]) / Δt
    # MATLAB uses assumptions -> doesn't conserve rotational momentum
    floes.p_dαdt[idx1] = 0
    floes.p_dαdt[idx2] = 0
    # Update pieces list of they exist
    if !isnothing(pieces_list) && length(pieces_list) > 0
        pieces_list.p_dαdt[pieces_idx:end] .= 0
        pieces_list.u[pieces_idx:end] .= new_u
        pieces_list.v[pieces_idx:end] .= new_v
        pieces_list.p_dxdt[pieces_idx:end] .= new_p_dxdt
        pieces_list.p_dydt[pieces_idx:end] .= new_p_dydt
        pieces_list.p_dudt[pieces_idx:end] .= new_p_dudt
        pieces_list.p_dvdt[pieces_idx:end] .= new_p_dvdt
    end
    # Update floes linear velocities and accelerations
    floes.u[idx1], floes.u[idx2] = new_u, new_u
    floes.v[idx1], floes.v[idx2] = new_v, new_v
    floes.p_dxdt[idx1], floes.p_dxdt[idx2] = new_p_dxdt, new_p_dxdt
    floes.p_dydt[idx1], floes.p_dydt[idx2] = new_p_dydt, new_p_dydt
    floes.p_dudt[idx1], floes.p_dudt[idx2] = new_p_dudt, new_p_dudt
    floes.p_dvdt[idx1], floes.p_dvdt[idx2] = new_p_dvdt, new_p_dvdt
    return
end

#=
Update a parent floes and its ghosts to match velocities and accelerations at
given index.
Inputs:
    floes       <StructArray{Floe}> list of floes in the simulation
    idx         <Int> index of floe within floes list to copy to floe and ghosts
    parent_idx  <Int> index of parent floe to update, along with its ghosts
=#
function update_ghost_timestep_vals!(floes, idx, parent_idx)
    floes.u[parent_idx] = floes.u[idx]
    floes.v[parent_idx] = floes.v[idx]
    floes.ξ[parent_idx] = floes.ξ[idx]
    floes.p_dxdt[parent_idx] = floes.p_dxdt[idx]
    floes.p_dydt[parent_idx] = floes.p_dydt[idx]
    floes.p_dudt[parent_idx] = floes.p_dudt[idx]
    floes.p_dvdt[parent_idx] = floes.p_dvdt[idx]
    floes.p_dαdt[parent_idx] = floes.p_dαdt[idx]
    floes.p_dξdt[parent_idx] = floes.p_dξdt[idx]
    for gidx in floes.ghosts[parent_idx]
        floes.u[gidx] = floes.u[parent_idx]
        floes.v[gidx] = floes.v[parent_idx]
        floes.ξ[gidx] = floes.ξ[parent_idx]
        floes.p_dxdt[gidx] = floes.p_dxdt[parent_idx]
        floes.p_dydt[gidx] = floes.p_dydt[parent_idx]
        floes.p_dudt[gidx] = floes.p_dudt[parent_idx]
        floes.p_dvdt[gidx] = floes.p_dvdt[parent_idx]
        floes.p_dαdt[gidx] = floes.p_dαdt[parent_idx]
        floes.p_dξdt[gidx] = floes.p_dξdt[parent_idx]
    end
    return
end

"""
    calc_stress!(floe)

Calculates the stress on a floe for current collisions given interactions and
floe properties.

## _Positional arguments_
- `floe::Floe`: single floe
- $FLOE_SETTINGS_DEF

##  _Returns_
- Nothing. Updates floe.stress_accum and floe.stress_instant in-place.
"""
function calc_stress!(floe::FloeType{FT}, floe_settings) where {FT}
    # Stress calcultions
    xi, yi = floe.centroid
    inters = floe.interactions
    # Calculates timestep stress
    stress = zeros(FT, 2, 2)
    if floe.num_inters > 0
        for i in 1:floe.num_inters
            stress[1, 1] += (inters[i, xpoint] - xi) * inters[i, xforce]
            stress[1, 2] += (inters[i, ypoint] - yi) * inters[i, xforce] +
                (inters[i, xpoint] - xi) * inters[i, yforce]
            stress[2, 2] += (inters[i, ypoint] - yi) * inters[i, yforce]
        end
        stress[1, 2] *= FT(0.5)
        stress[2, 1] = stress[1, 2]
        stress .*= 1/(floe.area * floe.height)
    end
    # Updates accumulated stress
    _update_stress_accum!(floe_settings.stress_calculator, stress, floe)
    # Updates instantanious stress
    floe.stress_instant .= stress
    return
end

"""
    calc_strain!(...)

Calculates the strain on a floe given the velocity at each vertex.

## _Positional arguments_
- `floe::Floe`: single floe

##  _Returns_
- `strain::Matrix{AbstractFloat}`: 2x2 matrix for floe strain 
"""
function calc_strain!(floe::FloeType{FT}) where {FT}
    fill!(floe.strain, zero(FT))
    # coordinates of floe centered at centroid
    trans_poly = _translate_poly(FT, floe.poly, -floe.centroid[1], -floe.centroid[2])::Polys{FT}
    local x1, y1
    for (i, p2) in enumerate(GI.getpoint(GI.getexterior(trans_poly)))
        x2, y2 = get_tuple_point(p2, FT)
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        xdiff, ydiff = x2 - x1, y2 - y1
        rad1, rad2 = sqrt(x1^2 + y1^2), sqrt(x2^2 + y2^2)
        θ1, θ2 = atan(y1, x1), atan(y2, x2)
        u1 = floe.u - floe.ξ * rad1 * sin(θ1)
        u2 = floe.u - floe.ξ * rad2 * sin(θ2)
        v1 = floe.u + floe.ξ * rad1 * cos(θ1)
        v2 = floe.u + floe.ξ * rad2 * cos(θ2)
        udiff, vdiff = u2 - u1, v2 - v1
        floe.strain[1, 1] += udiff * ydiff
        floe.strain[1, 2] += udiff * xdiff + vdiff * ydiff
        floe.strain[2, 2] += vdiff * xdiff
        x1, y1 = x2, y2
    end
    floe.strain[1, 2] *= FT(0.5)
    floe.strain[2, 1] = floe.strain[1, 2]
    floe.strain ./= 2floe.area
    return
end

#=
Calculates the strain on floe i of a `FixedWidthFloes`. Same calculation as
`calc_strain!(floe)`, written as scalar operations so it can be called from a GPU kernel.
=#
function calc_strain!(floes::FixedWidthFloes{FT}, i) where FT
    cx = floes.centroid[i, 1]
    cy = floes.centroid[i, 2]
    u = floes.u[i]
    ξ = floes.ξ[i]
    ε11, ε12, ε22 = zero(FT), zero(FT), zero(FT)
    # coordinates of floe centered at centroid
    x1 = floes.poly[i, 1, 1] - cx
    y1 = floes.poly[i, 1, 2] - cy
    for j in 2:floes.n_points[i]
        x2 = floes.poly[i, j, 1] - cx
        y2 = floes.poly[i, j, 2] - cy
        xdiff, ydiff = x2 - x1, y2 - y1
        rad1, rad2 = sqrt(x1^2 + y1^2), sqrt(x2^2 + y2^2)
        θ1, θ2 = atan(y1, x1), atan(y2, x2)
        u1 = u - ξ * rad1 * sin(θ1)
        u2 = u - ξ * rad2 * sin(θ2)
        v1 = u + ξ * rad1 * cos(θ1)
        v2 = u + ξ * rad2 * cos(θ2)
        udiff, vdiff = u2 - u1, v2 - v1
        ε11 += udiff * ydiff
        ε12 += udiff * xdiff + vdiff * ydiff
        ε22 += vdiff * xdiff
        x1, y1 = x2, y2
    end
    ε12 *= FT(0.5)
    area2 = 2floes.area[i]
    floes.strain[i, 1, 1] = ε11 / area2
    floes.strain[i, 1, 2] = ε12 / area2
    floes.strain[i, 2, 1] = ε12 / area2
    floes.strain[i, 2, 2] = ε22 / area2
    return
end

@kernel function calc_strain_kernel!(floes)
    i = @index(Global)
    calc_strain!(floes, i)
end

# Bits in FixedWidthFloes.flags, set by kernels to report events to the host
const FLAG_HEIGHT_LIMITED = 0x01
const FLAG_FORCE_REDUCED = 0x02
const FLAG_VELOCITY_ADJUSTED = 0x04
const FLAG_ξ_SHRUNK = 0x08

@kernel function limit_height_kernel!(floes, max_height)
    # Ensure no extreme height values due to model instability
    i = @index(Global)
    if floes.height[i] > max_height
        floes.height[i] = max_height
        floes.flags[i] |= FLAG_HEIGHT_LIMITED
    end
end

#=
Ensure no extreme collision forces due to model instability. Reduces the collision force and
torque of floe i in place, so it must run before the floe's mass changes.
=#
function _limit_collision_force!(floes::FixedWidthFloes, i, Δt)
    cfx = floes.collision_force[i, 1]
    cfy = floes.collision_force[i, 2]
    ctrq = floes.collision_trq[i]
    while max(abs(cfx), abs(cfy)) > floes.mass[i]/(5Δt)
        floes.flags[i] |= FLAG_FORCE_REDUCED
        cfx /= 10
        cfy /= 10
        ctrq /= 10
    end
    floes.collision_force[i, 1] = cfx
    floes.collision_force[i, 2] = cfy
    floes.collision_trq[i] = ctrq
    return
end

@kernel function limit_collision_force_kernel!(floes, Δt)
    i = @index(Global)
    _limit_collision_force!(floes, i, Δt)
end

@kernel function thermodynamic_growth_kernel!(floes)
    # Update floe based on thermodynamic growth
    i = @index(Global)
    h = floes.height[i]
    Δh = floes.hflx_factor[i] / h
    hfrac = (h + Δh) / h
    floes.mass[i] *= hfrac
    floes.moment[i] *= hfrac
    floes.height[i] -= Δh
end

@kernel function update_ice_coordinates_kernel!(floes, Δt)
    # Update ice coordinates with velocities and rotation
    i = @index(Global)
    Δx = 1.5Δt*floes.u[i] - 0.5Δt*floes.p_dxdt[i]
    Δy = 1.5Δt*floes.v[i] - 0.5Δt*floes.p_dydt[i]
    Δα = 1.5Δt*floes.ξ[i] - 0.5Δt*floes.p_dαdt[i]
    floes.α[i] += Δα

    _move_floe!(floes, i, Δx, Δy, Δα)
    floes.p_dxdt[i] = floes.u[i]
    floes.p_dydt[i] = floes.v[i]
    floes.p_dαdt[i] = floes.ξ[i]
end

# Update ice velocities of floe i with forces and torques
function _update_velocities!(floes::FixedWidthFloes{FT}, i, Δt, maximum_ξ) where FT
    h = floes.height[i]
    dudt = (floes.fxOA[i] + floes.collision_force[i, 1])/floes.mass[i]
    dvdt = (floes.fyOA[i] + floes.collision_force[i, 2])/floes.mass[i]
    frac = if abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) > (h/2)
        frac1 = (sign(dudt)*h/2Δt)/dudt
        frac2 = (sign(dvdt)*h/2Δt)/dvdt
        min(frac1, frac2)
    elseif abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) < (h/2)
        (sign(dudt)*h/2Δt)/dudt
    elseif abs(Δt*dudt) < (h/2) && abs(Δt*dvdt) > (h/2)
        (sign(dvdt)*h/2Δt)/dvdt
    else
        one(FT)
    end
    if frac != 1
        floes.flags[i] |= FLAG_VELOCITY_ADJUSTED
        dudt = frac*dudt
        dvdt = frac*dvdt
    end
    floes.u[i] += 1.5Δt*dudt-0.5Δt*floes.p_dudt[i]
    floes.v[i] += 1.5Δt*dvdt-0.5Δt*floes.p_dvdt[i]
    floes.p_dudt[i] = dudt
    floes.p_dvdt[i] = dvdt

    dξdt = (floes.trqOA[i] + floes.collision_trq[i])/floes.moment[i]
    dξdt = frac*dξdt
    ξ = floes.ξ[i] + 1.5Δt*dξdt-0.5Δt*floes.p_dξdt[i]
    if abs(ξ) > maximum_ξ
        floes.flags[i] |= FLAG_ξ_SHRUNK
        ξ = sign(ξ) * maximum_ξ
    end
    floes.ξ[i] = ξ
    floes.p_dξdt[i] = dξdt
    return
end

@kernel function update_velocities_kernel!(floes, Δt, maximum_ξ)
    i = @index(Global)
    _update_velocities!(floes, i, Δt, maximum_ξ)
end

# Log the events that kernels reported in flags, once per type of event
function _log_flags(flags, tstep, floe_settings)
    nfloes(flag) = count(f -> f & flag != 0, flags)
    n = nfloes(FLAG_HEIGHT_LIMITED)
    n > 0 && @info "Reducing height to $(floe_settings.max_floe_height) m" tstep nfloes = n
    n = nfloes(FLAG_FORCE_REDUCED)
    n > 0 && @info "Decreasing collision forces by a factor of 10" tstep nfloes = n
    n = nfloes(FLAG_VELOCITY_ADJUSTED)
    n > 0 && @info "Adjusting u and v velocities to prevent too high" tstep nfloes = n
    n = nfloes(FLAG_ξ_SHRUNK)
    n > 0 && @info "Shrinking ξ" tstep nfloes = n
    return
end

"""
    timestep_floe_properties!(...)

Update floe position and velocities using second-order time stepping with
tendencies calculated at previous timesteps. Height, mass, stress, and strain
also updated based on previous timestep thermodynamics and interactions with
other floes. 

## _Positional arguments_
- $FLOES_DEF
- `tstep::Int`: simulation timestep
- $ΔT_DEF
- $FLOE_SETTINGS_DEF
"""
function timestep_floe_properties!(
    floes::StructArray{<:Floe{FT}},
    tstep,
    Δt,
    floe_settings,
) where FT
    Threads.@threads for i in eachindex(floes)
        # Update stress
        calc_stress!(get_floe(floes, i), floe_settings)
    end

    dev_floes = adapt(CuArray, FixedWidthFloes(floes))
    dev = get_backend(dev_floes.height)

    limit_height_kernel!(dev, 512)(dev_floes, floe_settings.max_floe_height, ndrange=size(floes.height))
    limit_collision_force_kernel!(dev, 512)(dev_floes, Δt, ndrange=size(floes.height))
    thermodynamic_growth_kernel!(dev, 512)(dev_floes, ndrange=size(floes.height))
    update_ice_coordinates_kernel!(dev, 512)(dev_floes, Δt, ndrange=size(floes.height))
    update_velocities_kernel!(dev, 512)(dev_floes, Δt, floe_settings.maximum_ξ, ndrange=size(floes.height))
    calc_strain_kernel!(dev, 512)(dev_floes, ndrange=size(floes.height))

    KernelAbstractions.synchronize(dev)
    host_floes = adapt(Array, dev_floes)
    update_floes!(floes, host_floes)
    _log_flags(host_floes.flags, tstep, floe_settings)
    return
end