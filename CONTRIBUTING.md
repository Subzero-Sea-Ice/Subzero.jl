# Subzero.jl contributor guide

Welcome to the Subzero.jl contributor guide.

If you are new to open source development, [here is a good guide to get started](https://github.com/firstcontributions/first-contributions). If you want something Julia specific, [check out this video](https://youtu.be/cquJ9kPkwR8).

## Local Repository
If you want to work on developing Subzero (other than for small documentation changes), you probably want to work locally on your computer. For this, you will want to create a fork of the repository, and then eventually open a pull request to the original code base. Here is a good [step-by-step guide of this process](https://www.matecdev.com/posts/julia-package-collaboration.html).

## Documentation
Contributing to the documentation is a great way to get involved in Subzero development. If something in the documents is confusing, please feel free to fix it! We always welcome improved documentation. 

Small changes can be done easily in GitHub's web interface (see [Editing
files](https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files#editing-files-in-another-users-repository)). Every page in the documentation has an `Edit on GitHub` button at the top, which takes you to the correct source file. The video [Making Julia documentation
better](https://youtu.be/ZpH1ry8qqfw) guides you through these steps.

If you want to edit larger sections of the documentation, you probably want to use the above instructions to work on a local repository. You can then make changes there and open a pull request to have them incorporated into the package. 

To see your changes live before opening a pull request, you can locally build the documentation. For this, you will need two terminal windows. Start by making sure you are in the `Subzero` folder in both terminals. Then, in both terminals, you will want to launch Julia.

In the first terminal:
```julia
pkg> activate docs
julia> include("docs/make.jl")
```

This will build the documents. It might take a little while.

```julia
pkg> activate docs
julia> using LiveServer
julia> serve(;dir = "docs/build")
```

This uses [`LiveServer.jl`](https://github.com/tlienart/LiveServer.jl) to launch a local webserver which you can visit at
[http://localhost:8000](http://localhost:8000). 

For more information on documentation see:

**Useful resources**
 - General information about documenting Julia code in the [Julia manual](https://docs.julialang.org/en/v1/manual/documentation/).
 - [Documentation for `Documenter.jl`](https://juliadocs.github.io/Documenter.jl/) which is used to render the HTML pages.
 - [Documentation for `Literate.jl`](https://fredrikekre.github.io/Literate.jl/v2/) which is used for tutorials/examples.

## Reporting issues

If you have found a bug or a problem with Subzero you can open an [issue](https://github.com/Subzero-Sea-Ice/Subzero.jl/issues/new). Try
to include as much information about the problem as possible and some code that
can be copy-pasted to reproduce it (see [How to create a Minimal, Reproducible
Example](https://stackoverflow.com/help/minimal-reproducible-example)).

If you can identify a fix for the bug you can submit a pull request without first opening an
issue, see [Code changes](#code-changes).

## Code changes

Bug fixes and improvements to the code, or to the unit tests are always welcome. If you have
ideas about new features or functionality it might be good to first open an
[issue](https://github.com/Subzero-Sea-Ice/Subzero.jl/issues/new) to get feedback before spending too much time implementing something.

When you are ready to make changes, check out the developer docs section of the documentation for insight into how the code is written and organized.

Remember to always include (when applicable): unit tests which exercises the new code, and updated documentation.

## Code Modularity

The most important thing to know about Subzero.jl is that the code has been written to be very modular, meaning that you are able to swap in and out different components of the `model`/`simulation` objects to gain new functionality. A user of the code should not need to modify the source code to add new science functionality. They should simply be able to add a new "method" to an existing "function" and swap that in seamlessly.

To do this, Subzero relies heavily on the Julia paradigm of multiple dispatch. If you don't know much about Julia, the first step would be to gain some intuition about [how to write good Julia code](https://modernjuliaworkflows.org/), how [multiple dispatch works](https://www.youtube.com/live/kc9HwsxE1OY), and [examples of multiple dispatch](https://www.matecdev.com/posts/julia-multiple-dispatch.html). 

You can see examples of this in the existing code, both in the model/simulation setup, as well as the scientific aspects of the code. For example, there are multiple types of boundary walls the user can choose from. A new user could easily add a new boundary type within their own code by declaring their new boundary a subtype of `AbstractBoundary` and implementing the two needed functions: `_update_boundary!(boundary::AbstractBoundary, Δt::Int)` and `_periodic_compat(boundary1::AbstractBoundary, boundary2::AbstractBoundary)`. To understand how multiple dispatch is set up within this code, please read through the [abstract_domains.jl](https://github.com/Subzero-Sea-Ice/Subzero.jl/blob/main/src/simulation_components/domain_components/abstract_domains.jl) file. 

Another example is the ability of the user to create and select new fracture criteria to determine when a floe within the simulation will break. This can be seen in the [fractures.jl](https://github.com/Subzero-Sea-Ice/Subzero.jl/blob/main/src/physical_processes/fractures.jl) file. Within the file, there are three subtypes of `AbstractFractureCriteria`, which each have a different method of the `update_criteria!` function.

With these examples in mind, if you are a developer are working with the code, and realize there is a functionality that the code doesn’t currently have, try to make it as modular as possible! If you want to add a new science functionality, this could mean adding a subtype to an existing abstract type and adding the needed methods for required functions (e.g. new fracture criteria) However, this could also mean adding in a new modular component to the code. The most important thing to remember is that if you start adding functionality to the code and start adding an if/else statement or some type of boolean flag, you should probably be using multiple dispatch instead! Make things modular where you can. 

## Acknowledgments

Thank you to the folks at [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) for having such amazing contributor docs. I took lots of inspiration and links from their wonderful, comprehensive page! Check them out for examples of documentation and guides done right!
