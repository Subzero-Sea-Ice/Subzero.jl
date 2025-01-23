using Documenter, Literate
using Subzero
import GeoInterface as GI
import GeometryOps as GO
import LibGEOS as LG


# Converting any files in the literate folder to markdown
LITERATE_INPUT = joinpath(@__DIR__, "literate")
LITERATE_OUTPUT = joinpath(@__DIR__, "src")

for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    ipath = joinpath(root, file)
    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    # generate the markdown file calling Literate
    Literate.markdown(ipath, opath)
end

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
