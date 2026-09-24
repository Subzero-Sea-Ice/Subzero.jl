using Documenter, Literate
using Subzero
using Random
import GeoInterface as GI
import GeometryOps as GO
import LibGEOS as LG

function jl_to_md(input, output)
    # turns .jl files in input to .md files in output
    for ipath in readdir(input, join = true)
        # ignore non julia files
        splitext(ipath)[2] == ".jl" || continue
        # # full path to a literate script
        # ipath = joinpath(root, file)
        # generated output path
        opath = splitdir(replace(ipath, input=>output))[1]
        # generate the markdown file calling Literate
        Literate.markdown(ipath, opath)
    end
end

# Copy over contributor guidlines and adjust EditURL to point to the true source file
open("docs/src/contribute.md", "w") do io
    editurl = """
    ```@meta
    EditURL = "../../CONTRIBUTING.md"
    ```
    """
    println(io, editurl)
    write(io, read(joinpath(@__DIR__, "..", "CONTRIBUTING.md")))
end

# Converting any files in the literate folder to markdown
println("Building tutorial...")
tutorial_input = joinpath(@__DIR__, "literate")
tutorial_output = joinpath(@__DIR__, "src")
jl_to_md(tutorial_input, tutorial_output)

println("Building examples...")
examples_input = joinpath(@__DIR__, "literate", "examples")
examples_output = joinpath(@__DIR__, "src", "examples")
jl_to_md(examples_input, examples_output)

# Documentation formatting
format = Documenter.HTML(;
    repolink = "https://github.com/Subzero-Sea-Ice/Subzero.jl",
    canonical = "https://Subzero-Sea-Ice.github.io/SubzeroDocumentation/stable",
    mathengine = MathJax3(),
    size_threshold = 5*10^6,  # 500 KiB
    collapselevel = 3,
)

# Metadata from doc tests
DocMeta.setdocmeta!(Subzero, :DocTestSetup, :(using Subzero); recursive=true)

println("Making docs...")

makedocs(;
    modules=[Subzero],
    authors="Skylar Gering and contributers",
    sitename="Subzero.jl",
    build = "build",
    format,
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "Examples" => [
            "examples/shear_flow.md",
            "examples/simple_strait.md",
            "examples/forcing_contained_floes.md",
            "examples/moving_bounds.md",
            "examples/restart_sim.md"
        ],
        "API Reference" => "api.md",
        "Contributing to Subzero" => "contribute.md",
    ],
    warnonly = Documenter.except(:doctest),
)

println("Adding .gitignore file...")
write("docs/build/.gitignore", "examples/**/*.jld2");

println("Deploying docs...")
deploydocs(;
    target = "build",
    repo="https://github.com/Subzero-Sea-Ice/Subzero.jl",
    push_preview = true,
)
