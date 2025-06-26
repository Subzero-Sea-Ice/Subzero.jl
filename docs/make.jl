using Documenter, Literate
using Subzero
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
    repolink = "https://github.com/Caltech-OCTO/Subzero.jl",
    canonical = "https://Caltech-OCTO.github.io/SubzeroDocumentation/stable",
    mathengine = MathJax3(),
    size_threshold = 5*10^5,  # 500 KiB
)

# Metadata from doc tests
DocMeta.setdocmeta!(Subzero, :DocTestSetup, :(using Subzero); recursive=true)

makedocs(;
    modules=[Subzero],
    authors="Skylar Gering and contributers",
    sitename="Subzero.jl",
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
        "Improving Subzero" => [
            "Contributing" => "contribute.md",
            "Developer Documentation" => "devdocs.md",
        ],
    ],
    warnonly = true,
)

deploydocs(;
    repo="https://github.com/Caltech-OCTO/Subzero.jl",
)
