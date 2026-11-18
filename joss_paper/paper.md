---
title: 'Subzero.jl: A Coupled, Julia Version of a Discrete-Element Sea Ice Model'
tags:
  - Julia
  - oceanography
  - sea ice
  - dynamics
authors:
  - name: Skylar A. Gering
    orcid: 0000-0003-1974-3966
    corresponding: true
    affiliation: "1, 2"
  - name: Mukund Gupta
    orcid: 0000-0003-0181-9504
    affiliation: "1, 3"
  - name: Samuel Brenner
    orcid: 0000-0002-0826-1294
    affiliation: 1
  - name: Andrew F. Thompson
    orcid: 0000-0003-0322-4811
    affiliation: 1
affiliations:
 - name: California Institute of Technology, USA
   index: 1
 - name: Massachusetts Institute of Technology, USA
   index: 2
 - name: Delft University of Technology, NL
   index: 3
date: 1 January 2026
bibliography: paper.bib
---

# Summary

Subzero.jl is a discrete-element model (DEM) for simulating the dynamics of two-dimensional polygonal sea-ice floes. The model is a port of the MATLAB-based SubZero model [@montemuro2025] to the Julia programming language [@bezanson2017], with substantial re-engineering to increase performance, improve extensibility and usability, and provide native Julia infrastructure for two-way coupling with the external ocean dynamics model Oceananigans.jl [@ramadhan2020; @wagner2025]. Subzero.jl is available as a registered Julia package.

# Statement of need

Arctic sea ice extent and concentration have declined at rates that are commonly underestimated by climate projection models [@notz2020; @shu2020]. Potential sources of uncertainty arise from an inaccurate representation of interactions between the ocean and sea ice within these climate models, as well as the assumption of continuum sea-ice dynamics for the sake of reducing computational complexity.

Discrete-element models (DEMs), where each piece of sea ice is represented as an individual simulation element, all of which can dynamically interact, are used to explore these uncertainties and study fine-scale sea ice dynamics. While a range of sea ice DEMs are available [@hopkins2004a; @herman2013numerical; @rabatel2015a; @herman2016discrete; @damsgaard2018], many of these models make a range of simplifications, particularly in the geometric representation of the ice floes.

The MATLAB-based SubZero model [@montemuro2025] allows for complex and evolving floe shapes, but is computationally expensive. Moreover, due to the highly connected nature of sea ice floes and ocean processes [@horvat2018; @gupta2022; @gupta2024; @brenner2023c], there is a need for DEM simulations to be coupled to a dynamic ocean model to explore two-way feedbacks, which is a feature not readily available for most extant models. 

We present Subzero.jl, a native Julia [@bezanson2017] version of the MATLAB discrete-element model SubZero [@montemuro2025; @manucharyan2022b] that addresses both of these problems. 

# Functionality

Subzero.jl represents sea ice floes as polygonal elements that move in response to forcing from the atmosphere and the ocean [@manucharyan2022b]. Over time, a given floe may change its horizontal shape, adjust vertical thickness, and fracture into multiple pieces as a result of interactions with other floes and topographical elements. Subzero.jl improves upon the MATLAB version in three major ways: (i) a modular interface that reduces the need for users to modify source code, (ii) the ability to couple to a performant native-Julia ocean model, and (iii) significant enhancement in computational speed.

## Modular interface

This new Julia implementation provides a more extensible, user-friendly interface while preserving the scientific behavior of the previous MATLAB model. In the previous MATLAB workflow, simulation-specific settings were managed at the repository level, requiring copying and adapting directories to create new experiments. However, Subzero.jl is organized around composable object types (`structs`) in a hierarchical structure, with a top-level `Simulation` object containing `Floes`, `Domain`, `Atmospheric Forcing`, `Ocean Forcing`, `Physical Settings` and `Output Writer` components. Using dedicated run scripts, the user may modify object parameters to easily turn on/off physical processes such as fracturing, ridging, and welding. Moreover, the use of Julia's multiple dispatch allows flexibility to add more functionality without changing the underlying model structure. This flexibility is highlighted in the tutorial and examples provided with the documentation. Additionally, Subzero.jl includes a suite of unit and integration tests that verify consistency with the original MATLAB model and ensure conservation of energy and momentum. These changes allow for more dynamic and flexible usage and development than in the MATLAB version, while ensuring parity with the original code.

## Ice-Ocean Coupling

The code now also includes a framework that enables two-way mechanically coupled simulations with the Julia ocean model Oceananigans.jl [@ramadhan2020; @wagner2025]. In this configuration, sea-ice floes in Subzero.jl are forced by gridded ocean velocity fields supplied by Oceananigans.jl and interpolated onto a mesh carried by each floe. In turn, floe-resolved ice-ocean stresses computed by Subzero.jl are interpolated onto the ocean grid and returned to the ocean model as a spatially heterogeneous and temporally evolving surface boundary condition. The exchange of fields is implemented using a callback functionality in Oceananigans.jl, allowing serial coupling at user-defined frequencies. This coupled DEM-LES system resolves the mechanical interaction between evolving ocean fields and individual sea ice floes, enabling investigation of floe-scale modulation of upper-ocean dynamics [@brenner2025]. An example coupled simulation is demonstrated in \autoref{fig:oceananigans_coupled_example}.


\begin{figure}
  \includegraphics[width=\linewidth]{oceananigans_coupled_example.png}
  \caption{Example Subzero-Oceananigans coupled simulation. Left: surface fields of ocean vorticity, overlain by sea ice floes, each coloured by their vorticity. Right: the same ocean vorticity field with floes not plotted; the impacts of ice-ocean coupling are still evident in the "patchiness" of the field. Vorticity units are given in terms of the Coriolis frequency, $f$, which has a value of $1.26 \times 10^{-4} s^{-1}$ at $60^\circ$N.}
  \label{fig:oceananigans_coupled_example}
\end{figure}

## Performance improvements
Beyond this new functionality, Subzero.jl provides significant performance improvements relative to the original MATLAB implementation of the model.
We used the \texttt{shear\_flow} example from the [Subzero.jl repository](https://subzero-sea-ice.github.io/Subzero.jl/dev/examples/shear_flow/) as a basis for testing speed enhancements.
The simulation includes 50 floes at 75% sea ice concentration, forced by sheared ocean currents in a doubly-periodic domain; we also created a version with 1000 floes, and then implemented both versions in the MATLAB SubZero codebase.
We tested the models on the California Institute of Technology's "Resnick" High Performance Computing system, on a single compute node using 1-16 CPU cores, with each configuration repeated five times.
Results in \autoref{fig:speed_comparison} show comparisons of both end-to-end runtime (including initialization and I/O) and simulation time (time spent advancing the model state).  


\begin{figure}
  \includegraphics[width=\linewidth]{speed_comparison.png}
  \caption{Runtime of \texttt{shear\_flow} simulation in Julia Subzero.jl (blue) and MATLAB SubZero (red) for comparison. Left: 50 floe simulation. Right: 1000 floe simulation. Solid lines show end-to-end time while dashed lines show strictly simulation time.}
  \label{fig:speed_comparison}
\end{figure}


Subzero.jl achieves substantial speedups relative to the MATLAB implementation: 7.2-13.0 times faster (end-to-end) and 17.4-31.1 times faster (simulation) for 50 floes, increasing to 32.2-48.5 times (end-to-end) and 43.9-61.1 times (simulation) for 1000 floes. Parallel scaling is modest in both implementations, with simulation time speedups of 1.1 times (Julia) and 1.8 times (MATLAB) for 50 floes, and 1.9 times (Julia) and 2.3 times (MATLAB) for 1000 floes when increasing from 1 to 16 cores. The speed increase of Subzero.jl relative to the original MATLAB model will enable longer and more complex simulations, particularly at higher floe counts.


# Acknowledgements

Our work is supported by the Office of Naval Research (ONR) grant N00014-19-1-2421. The authors thank the creators of MATLAB code SubZero, Georgy Manucharyan and Brandon Montemuro, for their guidance and advice during the porting process. 


S.A.G. would also like to acknowledge funding support from the Schmidt Academy for Software Engineering from 2022-2024 and the MIT Presidential Fellowship for 2024-2025. The manuscript was completed in 2026 while S.A.G. was supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Department of Energy Computational Science Graduate Fellowship under Award Number DE-SC0026073. 



# References

