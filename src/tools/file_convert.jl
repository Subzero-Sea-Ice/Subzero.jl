"""
Code for converting between .MAT and .jld2 filetypes
"""
# Uses JLD2.jl and SplitApplyCombine.jl

# Given .mat file with a struct that holds an array of polyshapes.
# This doesn't open well with Julia because polyshape doesn't exist. 
# Do a version of the following commented-out code in MATLAB.
# It opens the file and resaves just the vertices in a cell array:

# F = load('FloeShapes.mat');
# polyVertices = cell(size(F, 2), 1)
# for i = 1:size(F.poly, 2)
#    polyVertices(i, 1) = {F.poly(i).Vertices}
# end
# save('FloeVertices.mat', 'polyVertices')

# Then switch to Julia, open the new file, transform the coordiantes into
# PolyVecs using SplitApplyCombine and re-save as a JLD2 file for use

function matfloe2julfloes(filename)
    filename = auto_extension(filename, ".mat")
    vars = matread(filename)
    coords = vars["floe"]["c_alpha"]
    coords = splitdims(coords)
    coords = Subzero.translate([coords], vars["floe"]["Xi"], vars["floe"]["Yi"])
    floe = Floe(coords, vars["floe"]["h"], 0.0)
    floe.centroid = [vars["floe"]["Xi"], vars["floe"]["Yi"]]
    floe.x_subfloe_points = vec(vars["floe"]["X"])[vec(vars["floe"]["A"]) .== 1]
    floe.y_subfloe_points = vec(vars["floe"]["Y"])[vec(vars["floe"]["A"]) .== 1]
    floe.u =vars["floe"]["Ui"]
    floe.v =vars["floe"]["Vi"]
    floe.ξ = vars["floe"]["ksi_ice"]
    floe.α = vars["floe"]["alpha_i"]
    floe.interactions = vars["floe"]["interactions"]
    floe.fxOA = vars["floe"]["FxOA"]
    floe.fyOA = vars["floe"]["FyOA"]
    floe.trqOA = vars["floe"]["torqueOA"]
    floe.collision_force = vars["floe"]["collision_force"]
    floe.collision_trq = vars["floe"]["collision_torque"]
    floe.stress_accum = vars["floe"]["Stress"]
    floe.strain = vars["floe"]["strain"]
    for x in eachslice(vars["floe"]["StressH"], dims = 3)
        if sum(x) != 0
            push!(floe.stress_instant, x)
        end
    end
    return floe
end



# JLD2 to MATLAB

function julfloe2matfloe(floes, Δg, out_fn)
    xcoords = first.(floes.coords[1])
    ycoords = last.(floes.coords[1])
    reshaped_x = Vector{Matrix{Float64}}()
    reshaped_y = Vector{Matrix{Float64}}()
    for i in range(1, length(xcoords))
        push!(reshaped_x, reshape(xcoords[i], (1, length(xcoords[i]))) .- Δg)
        push!(reshaped_y, reshape(ycoords[i], (1, length(ycoords[i]))) .- Δg)
    end

    floe.stress_accum = vars["floe"]["Stress"]
    floe.strain = vars["floe"]["strain"]
    for x in eachslice(vars["floe"]["StressH"], dims = 3)
        if sum(x) != 0
            push!(floe.stress_instant, x)
        end
    end
    out_fn = auto_extension(out_fn, ".mat")
    matwrite(out_fn, Dict(
        "x_subfloe_points" => floes.x_subfloe_points,
        "y_subfloe_points" =>  floes.y_subfloe_points,
        "xcoords" => reshaped_x,
        "ycoords" => reshaped_y,
        "u" => floes.u,
        "v" => floes.v,
        "ksi_ice" => floes.ξ,
        "alpha_i" => floes.α,
        "interactions" => floes.interactions,
        "fxOA" => floes.fxOA,
        "fyOA" => floes.fyOA,
        "torqueOA" => floes.trqOA,
        "collision_force" => floes.collision_force,
        "collision_toque" => floes.collision_trq,
        "stress" => floes.stress_accum,
        "strain" => floes.strain,
        "stress_instant" => floes.stress_instant
    );)
end

    matwrite(out_fn, Dict(
        "x_subfloe_points" => floes.x_subfloe_points,
        "y_subfloe_points" =>  floes.y_subfloe_points,
        "xcoords" => reshaped_x,
        "ycoords" => reshaped_y,
    );)