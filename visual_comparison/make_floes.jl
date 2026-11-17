#=
Make the floes of an example once and save them, so that both versions start from the same
floes: julia --project=env-original make_floes.jl <example> <file.jld2>
=#
using Subzero, JLD2, Random, Statistics, StructArrays
const GI = Subzero.GI
include(joinpath(@__DIR__, "examples.jl"))

name, file = ARGS
floes = make_example_floes(name)
jldsave(file; floes = floes_to_dict(floes))
println("$name: $(length(floes)) floes saved to $file")
