using JLD2, NCDatasets

"""
    compare_floe_data(filename1, filename2)

Compare two files output by the floe output writer. Prints out first instances
of not matching per timestep and field. 
Inputs:
    filename1   <String> filename and path of first file
    filename2   <String> filename and path of second file
Outputs:
    If there are instances of differences, function will print time and index of
    floes that don't match
"""
function compare_floe_data(filename1, filename2)
    file1 = jldopen(filename1, "r")
    file2 = jldopen(filename2, "r")
    @assert length(keys(file1)) == length(keys(file2))
    for k in keys(file1)
        if k!= "metadata"
            println(string("Comparing key: ", k))
            for t in keys(file1[k])
                f1_data = file1[k][t]
                f2_data = file2[k][t]
                same = true
                for i in eachindex(f1_data)
                    if k == "stress_instant"
                        same = same && f1_data[i].total == f2_data[i].total
                        same = same && all(f1_data[i].cb .== f2_data[i].cb)
                    else
                        same = same && all(f1_data[i] .== f2_data[i])
                    end
                    if !same
                        println(string(
                            "Differences starting at time ",
                            t,
                            " and index ",
                            i,
                        ))
                        break
                    end
                end
            end
        end
    end
    close(file1)
    close(file2)
    return
end

"""
    compare_grid_data(filename1, filename2)
Compare two files output by the grid output writer. Prints out first instances
of not matching per field. All timesteps are compared at once. 
Inputs:
    filename1   <String> filename and path of first file
    filename2   <String> filename and path of second file
Outputs:
    If there are instances of differences, function will print the field that
    has discrepancies. 
"""
function compare_grid_data(filename1, filename2)
    file1 = Dataset(filename1, "r")
    file2 = Dataset(filename2, "r")
    @assert length(keys(file1)) == length(keys(file2))
    for k in keys(file1)
        if k!= "time" && k!= "x" && k!= "y"
            println(string("Comparing key: ", k))
            f1_data = file1[k]
            f2_data = file2[k]
            same = true
            same = same && all(f1_data[:, :, :] .== f2_data[:, :, :])
            if !same
                println(string("Differences in ", k))
                break
            end
        end
    end
    close(file1)
    close(file2)
end

"""
    compare_checkpointer_data(filename1, filename2)
Compare two files output by the checkpointer output writer. Compares ocean and
atmosphere. If there are discrepancies between the files, it will
timesteps and field. 
Inputs:
    filename1   <String> filename and path of first file
    filename2   <String> filename and path of second file
Outputs:
    If there are instances of differences, function will print the field and
    timesteps that have discrepancies. 
"""
function compare_oa_checkpointer_data(filename1, filename2)
    file1 = jldopen(filename1, "r")
    file2 = jldopen(filename2, "r")

    ocean1 = file1["ocean"]
    ocean2 = file2["ocean"]
    fields =fieldnames(Ocean)
    println("Comparing Ocean")
    for t in keys(ocean1)
        for f in fields
            same = true
            same = same && all(getfield(ocean1[t], f) .== getfield(ocean2[t], f))
            if !same
                println(string("Differences in field ", f, " and timestep ", t))
                break
            end
        end
    end

    atmos1 = file1["atmos"]
    atmos2 = file2["atmos"]
    fields =fieldnames(Atmos)
    println("Comparing Atmosphere")
    for t in keys(ocean1)
        for f in fields
            same = true
            same = same && all(getfield(atmos1[t], f) .== getfield(atmos2[t], f))
            if !same
                println(string("Differences in field ", f, "and timestep ", t))
                break
            end
        end
    end

    close(file1)
    close(file2)
end


