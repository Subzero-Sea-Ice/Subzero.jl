module SubzeroMatExt

using Subzero, MAT, JLD2, GeoInterface
import Subzero: polygons_to_matlab, compare_matlab_julia

function polygons_to_matlab(jl_file, mat_file)
    # get julia polygons where jl_file is a floes file output by a simulation
    f = jldopen(jl_file)
    polys = f["poly"]["0"]
    close(f)

    n = length(polys)
    # MATLAB-style cell array: Vector of matrices
    poly_cells = Vector{Matrix{Float64}}(undef, n)

    for i in range(1, n)
        poly = polys[i]

        # Get exterior ring
        ring = GeoInterface.getexterior(poly)

        # Extract coordinates
        coords = GeoInterface.getpoint(ring)

        # Split into x/y
        x = first.(coords)
        y = last.(coords)

        # MATLAB expects Nx2 for polyshape
        poly_cells[i] = hcat(x, y)
    end

    # Write MAT file
    matwrite(mat_file, Dict("polygons" => poly_cells))
end

function compare_matlab_julia(mat_dir, jl_floe_file, nout, Nb = 0)
    # get MATLAB floes
    files = readdir(mat_dir)
    pattern = r"^Floe0*(\d+)\.mat$"
    idx = Int[]
    mat_files = String[]
    for f in files
        m = match(pattern, f)
        if m !== nothing
            push!(idx, parse(Int, m.captures[1]))
            push!(mat_files, f)
        end
    end
    sort!(idx)
    sort!(mat_files)

    nfloes_mat = Int[]
    ninters_mat = Int[]
    for fname in mat_files
        path = joinpath(mat_dir, fname)
        data = matread(path)
        Floe = data["Floe"]
        # Example field access
        area = Floe["area"]
        push!(nfloes_mat, length(area))
        n = sum(size.(Floe["interactions"], 1))
        push!(ninters_mat, n) # all are repeated since both floes involved record
    end

    # get Julia floes
    out_times = nout * idx
    jl_file = jldopen(jl_floe_file)

    nfloes_julia = Int[]
    ninters_julia = Int[]
    for time in ["0"; out_times[1:end-1]]
        time_str = string(time)
        
        ghost_id = jl_file["ghost_id"][time_str]
        poly     = jl_file["poly"][time_str]
        num_inters = jl_file["num_inters"][time_str]
        
        non_ghost_idx = ghost_id .== 0
        
        # Filter arrays to get non-ghost floes
        nfloes_julia_push = length(poly[non_ghost_idx])
        ninters_julia_push = sum(num_inters[non_ghost_idx])
        
        push!(nfloes_julia, nfloes_julia_push)
        push!(ninters_julia, ninters_julia_push)
    end

    return nfloes_mat, ninters_mat, nfloes_julia, ninters_julia
end

end  # module