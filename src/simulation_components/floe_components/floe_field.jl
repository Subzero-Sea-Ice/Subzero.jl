export initialize_floe_field, AbstractFloeFieldGenerator, CoordinateListFieldGenerator, VoronoiTesselationFieldGenerator

"""
    YourFloeFieldGenerator{FT} <: AbstractFloeFieldGenerator{FT}

Each simulation run with requires a field of ice floes, which is simply a list of [`Floe`](@ref) structs.
However, there are various ways to create a floe field, and various configurations of floes that
user might want to create. 

To make it easier for the user to define a new floe creation functionality, this process uses multiple dispatch
on subtypes of the `AbstractFloeFieldGenerator` type. The generators hold the information needed to generate a floe
field - and the method is determined by the generator type. Right now, there are two subtypes implemented.

The first is the `CoordinateListFieldGenerator`. This generator takes in a list of floe coordinates (along with other arguments)
to create a list of floes. The second is the `VoronoiTesselationFieldGenerator`. This generator tesselates the domain with floes and then randomly keeps
various floes to match the user-provided floe concentration. 

Users and developers can create new generators to easily implement new floe generators. For example, it would be good to have a
generator that creates floes with a user-defined FSD.

## _API_
The following methods must be implemented for all subtypes:
- `_initialize_floe_field!(::Type{FT}, floe_arr, generator::CoordinateListFieldGenerator, domain; kwargs...)`: 
"""
abstract type AbstractFloeFieldGenerator{FT <: AbstractFloat}  end

"""
    initialize_floe_field(FT; generator, domain, supress_warnings, floe_settings, rng, kwargs...)

Create a field of floes using the provided generator. Regardless of generator, the floe array will be a vector of
[`Floe`](@ref) structs, but the exact shape/location/number of the floes will be determined by the generator, which must
be a subtype of the [`AbstractFloeFieldGenerator`](@ref). After floe creation, the unique, numerical ids of the floes within
the floe field will be set.

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `generator::AbstractFloeFieldGenerator`: generator type that determines how the floe field is generated
- `domain::Domain`: simulation domain
- `floe_settings::FloeSettings`: floe settings that determine individual floe characteristics - also needed to create a `Simulation` object
- `supress_warnings::Bool`: boolean flag on if warnings regarding floe area and centroid position should checked
- `rng::RandomNumberGenerator`:: random number generator needed to randomly assign floe characteristics (e.g height), depending on the generator
- `floe_bounds::Polys`: bounding box for floes to be generated within - only used for `VoronoiTesselationFieldGenerator`
- `kwargs...`: other keywords night be needed if new generator types are implemented. 
"""
function initialize_floe_field(
    ::Type{FT} = Float64;
    generator::AbstractFloeFieldGenerator{FT}, domain, floe_settings,
    supress_warnings = false, rng = Xoshiro(),
    kwargs...,
) where FT
    # initialize empty floe field array
    floe_arr = StructArray{Floe{FT}}(undef, 0)
    # add floes using whichever process
    _initialize_floe_field!(FT, floe_arr, generator, domain; floe_settings, rng, kwargs...)
    # warn about floes that don't match floe settings
    !supress_warnings && _warn_floe_validity(FT, floe_arr, domain, floe_settings)
    # add numerical ids to all floes
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end

#=
    _warn_floe_validity(FT, floe_arr, domain, floe_settings)

After floe field is created, check if the area of each of the floes is greater than the minimum floe area define by
the `floe_settings` and if the floe centroids are all within the domain. 

Additional checks could be added, or this function could dispatch off of the generator in the future if this needs to
be mroe specific to the floe generation method.
=#
function _warn_floe_validity(::Type{FT}, floe_arr, domain, floe_settings) where FT
    # find a reasonable minimum floe area given user input OR domain size
    min_floe_area = floe_settings.min_floe_area
    if floe_settings.min_floe_area ≤ 0
        min_floe_area = FT(4 * (domain.east.val - domain.west.val) * (domain.north.val - domain.south.val) / 1e4)
    end
    # warn about floes that are smaller than suggested minimum area
    if any(floe_arr.area .< min_floe_area)
        @warn "Some user input floe areas are less than the suggested minimum floe area."
    end
        # Warn about floes with centroids outside of domain
    floes_not_in_left_right = !all(domain.west.val .< first.(floe_arr.centroid) .< domain.east.val)
    floes_not_in_top_bottom = !all(domain.south.val .< last.(floe_arr.centroid) .< domain.north.val)
    if floes_not_in_left_right || floes_not_in_top_bottom
        @warn "Some floe centroids are out of the domain."
    end
    return
end

# Concrete subtype of AbstractFloeFieldGenerator - see below for documentation
struct CoordinateListFieldGenerator{FT} <: AbstractFloeFieldGenerator{FT}
    polys::Vector{Polys{FT}}
    hmean::FT
    Δh::FT
end

"""
    CoordinateListFieldGenerator{FT} <: AbstractFloeFieldGenerator{FT}

A concrete implementatin of [`AbstractFloeFieldGenerator`](@ref) that generates a floe field from a list of floe coordinates.
The floe coordinates are turned into polygons. The heights are a random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`.

The fields are:

- `polys::Vector{Polys{FT}}`: list of polygons in the form of [`PolyVec`](@ref)s. 
- `hmean::FT`: mean height of the generated floe field.
- `Δh::FT`: range of floe heights around the mean (i.e. random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`)

Here is how to construct a `CoordinateListFieldGenerator`:

    CoordinateListFieldGenerator([FT = Float64]; coords, hmean, Δh)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `coords::Vector{PolyVec{FT}}`: list of polygon coordinates in the form of [`PolyVec`](@ref)s
- `hmean::FT`: mean height of the generated floe field.
- `Δh::FT`: range of floe heights around the mean (i.e. random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`)
- `kwargs...`

!!! note
    It is now strongly reccomended to generate floe fields by first making a generator, but the functionality to create floes from a 
    list of coordiantes existed prior to the `AbstractFloeFieldGenerator` multiple dispatch setup. Therefore, the user can also call
    `initialize_floe_field(FT, coords, domain, hmean, Δh; kwargs...)` and that will automatically create a `CoordinateListFieldGenerator`
    to preserve existing functionality. However, this is NOT reccomended and only exists for backwards compatability. See source code if needed.

## _Examples_

```jldoctest
julia> coords = [[[(0.0, 0.0), (0.0, 100.0), (100.0, 100.0), (0.0, 0.0)]], [[(150.0, 0.0), (150.0, 100.0), (200.0, 100.0), (150.0, 0.0)]]];

julia> generator = CoordinateListFieldGenerator(Float64; coords, hmean = 1, Δh = 0.25)
CoordinateListFieldGenerator:
  ⊢Number of polygons: 2
  ⊢Mean height: 1.0
  ∟Height range: 0.25
```
"""
function CoordinateListFieldGenerator(::Type{FT} = Float64; coords, hmean, Δh) where FT
    polys = [make_polygon(valid_polyvec!(c), FT) for c in coords]
    return CoordinateListFieldGenerator{FT}(polys, hmean, Δh)
end

# Pretty printing for CoordinateListFieldGenerator
function Base.show(io::IO, generator::CoordinateListFieldGenerator{FT}; digits = 5) where FT
    overall_summary = "CoordinateListFieldGenerator:"
    total_number_summary = "Number of polygons: $(length(generator.polys))"
    mean_height_summary = "Mean height: $(round(mean(generator.hmean); digits))"
    variance_height_summary = "Height range: $(round(mean(generator.Δh); digits))"
    print(io, overall_summary, "\n",
        "  ⊢", total_number_summary, "\n",
        "  ⊢", mean_height_summary, "\n",
        "  ∟", variance_height_summary)
end

#=
    initialize_floe_field(FT, coords, domain, hmean, Δh; kwargs...)

Function to preserve backwards compatability that allows the creation of a floe field from coordiantes without first
making a CoordinateListFieldGenerator.
=#
function initialize_floe_field(::Type{FT}, coords::Vector{<:PolyVec}, domain, hmean, Δh; kwargs...) where FT
    generator = CoordinateListFieldGenerator(FT; coords, hmean, Δh)
    return initialize_floe_field(FT; generator, domain, kwargs...)
end

#=
    initialize_floe_field(FT, coords, domain, hmean, Δh; kwargs...)

Function to preserve backwards compatability that allows the creation of a floe field from coordiantes without first
making a CoordinateListFieldGenerator. Default value of `FT = Float64`.
=#
initialize_floe_field(coords::Vector{<:PolyVec}, args...; kwargs...) = initialize_floe_field(Float64, coords, args...; kwargs...)

#=
Create a field of floes from a list of polygons. This function dispatches off of CoordinateListFieldGenerator.
It also removes overlaps with topography and then just passes each polygon to _poly_to_floes!.
=#
function _initialize_floe_field!(
    ::Type{FT}, floe_arr, generator::CoordinateListFieldGenerator, domain;
    floe_settings, rng = Xoshiro(), kwargs...,
) where FT
    floe_polys = generator.polys
    # Remove overlaps with topography
    if !isempty(domain.topography)
        floe_polys = diff_polys(make_multipolygon(generator.polys), make_multipolygon(domain.topography.poly), FT)
    end
    # Turn polygons into floes
    for p in floe_polys
        _poly_to_floes!(FT, floe_arr, p, generator.hmean, generator.Δh, domain.east.val - domain.west.val;
            floe_settings, rng,
        )
    end
    floe_arr.id .= range(1, length(floe_arr))
    calculate_bonds!(floe_arr)
end

const FLOE_CONC_STR = "Concentrations per cell in a VoronoiTesselationFieldGenerator must be between 0 and 1. Values less than zero will be set to 0 and those greater than 1 will be set to 1."

# Concrete subtype of AbstractFloeFieldGenerator - see below for documentation
struct VoronoiTesselationFieldGenerator{FT} <: AbstractFloeFieldGenerator{FT}
    nfloes::Int
    concentrations::Matrix{FT}
    hmean::FT
    Δh::FT
    lloyd_iters::Int

    function VoronoiTesselationFieldGenerator{FT}(nfloes, concentrations, hmean, Δh, lloyd_iters = 0) where FT
        if !all(0 .<= concentrations .<= 1)
            @warn FLOE_CONC_STR
        end
        if lloyd_iters < 0
            @warn "lloyd_iters must be non-negative. Setting to 0."
            lloyd_iters = 0
        end
        clamped_concentrations = clamp.(concentrations, 0, 1)
        new{FT}(nfloes, clamped_concentrations, hmean, Δh, lloyd_iters)
    end
end

"""
    VoronoiTesselationFieldGenerator{FT} <: AbstractFloeFieldGenerator{FT}

A concrete implementatin of [`AbstractFloeFieldGenerator`](@ref) that generates a floe field using voronoi tesselation according to user
input number of floes and concentrations. The heights are a random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`.

The fields are:

- `nfloes::Int`:number of floes to try to create - note you might not end up with this number of floes.
        Topography in domain and multiple concentrations can decrease number of floes created
- `concentrations::Matrix{FT}`: matrix of concentrations to fill domain.
        If size(concentrations) = N, M then split the domain into NxM cells, each to be filled with the corresponding concentration. If concentration is
        below 0, it will default to 0. If it is above 1, itwill default to 1
- `hmean::FT`: mean height of the generated floe field.
- `Δh::FT`: range of floe heights around the mean (i.e. random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`)

Here is how to construct a `VoronoiTesselationFieldGenerator`:

    VoronoiTesselationFieldGenerator([FT = Float64]; nfloes, concentrations, hmean, Δh)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `nfloes::Int`:number of floes to try to create (see above for more info)
- `concentrations::Matrix{FT}`: matrix of concentrations to fill domain (see above for more info)
- `hmean::FT`: mean height of the generated floe field.
- `Δh::FT`: range of floe heights around the mean (i.e. random distribution around `hmean` between `hmean - Δh` and `hmean + Δh`)
- floe_bounds
- `kwargs...`

!!! note
    It is now strongly reccomended to generate floe fields by first making a generator, but the functionality to create floes from a 
    voronoi tesselation existed prior to the `AbstractFloeFieldGenerator` multiple dispatch setup. Therefore, the user can also call
    `initialize_floe_field(FT, nfloes, concentrations, hmean, Δh; kwargs...)` and that will automatically create a `VoronoiTesselationFieldGenerator`
    to preserve existing functionality. However, this is NOT reccomended and only exists for backwards compatability. See source code if needed.

## _Examples_

```jldoctest
julia> nfloes = 20;

julia> concentrations = [0.5];

julia> generator = VoronoiTesselationFieldGenerator(Float64; nfloes, concentrations, hmean = 1, Δh = 0.25)
VoronoiTesselationFieldGenerator:
  ⊢Requested number of floes: 20
  ⊢Requested concentrations [0.5;;]
  ⊢Mean height: 1.0
  ∟Height range: 0.25
```
"""
VoronoiTesselationFieldGenerator(::Type{FT} = Float64; nfloes, concentrations, hmean, Δh, lloyd_iters = 0, kwargs...) where FT =
    VoronoiTesselationFieldGenerator{FT}(nfloes, concentrations[:, :], hmean, Δh, lloyd_iters)

# Pretty printing for VoronoiTesselationFieldGenerator
function Base.show(io::IO, generator::VoronoiTesselationFieldGenerator{FT}; digits = 5) where FT
    overall_summary = "VoronoiTesselationFieldGenerator:"
    total_number_summary = "Requested number of floes: $(generator.nfloes)"
    concentration_summary = "Requested concentrations $(generator.concentrations)"
    mean_height_summary = "Mean height: $(round(mean(generator.hmean); digits))"
    variance_height_summary = "Height range: $(round(mean(generator.Δh); digits))"
    lloyd_summary = "Lloyd iterations: $(generator.lloyd_iters)"
    print(io, overall_summary, "\n",
        "  ⊢", total_number_summary, "\n",
        "  ⊢", concentration_summary, "\n",
        "  ⊢", mean_height_summary, "\n",
        "  ∟", variance_height_summary, "\n",
        "  ∟", lloyd_summary)
end

#=
    initialize_floe_field(FT, nfloes, concentrations, domain, hmean, Δh; kwargs...)

Function to preserve backwards compatability that allows the creation of a floe field from voronoi tesselation without first
making a VoronoiTesselationFieldGenerator.
=#
function initialize_floe_field(
    ::Type{FT}, nfloes::Int, concentrations, domain, hmean, Δh;
    kwargs...
) where FT
    generator = VoronoiTesselationFieldGenerator(FT; nfloes, concentrations, hmean, Δh, lloyd_iters)
    return initialize_floe_field(FT; generator, domain, kwargs...)
end

#=
    initialize_floe_field(FT, nfloes, concentrations, domain, hmean, Δh; kwargs...)

Function to preserve backwards compatability that allows the creation of a floe field from voronoi tesselation without first
making a VoronoiTesselationFieldGenerator. Default value of `FT = Float64`.
=#
initialize_floe_field(nfloes::Int, args...; kwargs...) = initialize_floe_field(Float64, nfloes, args...; kwargs...)

#=
Create a field of floes using voronoi tesselation. This function dispatches off of VoronoiTesselationFieldGenerator.
It also removes overlaps with topography and then tries to satisy each quardants concentration and also match the requested
total number of floes as closely as possible.
=#
function _initialize_floe_field!(
    ::Type{FT}, floe_arr, generator::VoronoiTesselationFieldGenerator, domain;
    floe_settings, rng,
    floe_bounds = _make_bounding_box_polygon(FT, domain.west.val, domain.east.val, domain.south.val, domain.north.val),
    kwargs...,
) where FT
    nfloes_added = 0 
    domain_poly = _make_bounding_box_polygon(FT, domain.west.val, domain.east.val, domain.south.val, domain.north.val)
    (bounds_xmin, bounds_xmax), (bounds_ymin, bounds_ymax) = GI.extent(domain_poly)
    FT_area_poly(p) = area_poly(p, FT)
    # Split domain into cells with given concentrations
    nrows, ncols = size(generator.concentrations)
    Lx = bounds_xmax - bounds_xmin
    Ly = bounds_ymax - bounds_ymin
    rowlen = Ly / nrows
    collen = Lx / ncols
    # Availible space in whole domain
    open_water = intersect_polys(floe_bounds, domain_poly, FT)
    if !isempty(domain.topography)
        open_water = diff_polys(make_multipolygon(open_water), make_multipolygon(domain.topography.poly), FT)
    end
    open_water_mp = make_multipolygon(open_water)
    open_water_area = area_poly(open_water_mp, FT)
    total_covered_water_area = (sum(generator.concentrations) / (nrows * ncols)) * open_water_area
    # Loop over cells
    for j in range(1, ncols)
        for i in range(1, nrows)
            c = generator.concentrations[i, j]
            if c > 0
                c = c > 1 ? 1 : c
                # Grid cell bounds
                xmin = bounds_xmin + collen * (j - 1)
                ymin = bounds_ymin + rowlen * (i - 1)
                # Open water in cell -> could make this into a helper function perhapes...
                cell_init = _make_bounding_box_polygon(FT, xmin, xmin + collen, ymin, ymin + rowlen)
                open_cell = intersect_polys(cell_init, open_water_mp, FT)
                open_cell_mpoly = make_multipolygon(open_cell)
                open_area = sum(FT_area_poly, open_cell; init = 0.0)
                # Generate coords with voronoi tesselation that fill the whole open space
                ncells = ceil(Int, generator.nfloes * ((open_area) / total_covered_water_area))
                floe_poly_list = _generate_voronoi_coords(FT, ncells, collen, rowlen, xmin, ymin,
                    open_cell_mpoly, rng, ncells; lloyd_iters = generator.lloyd_iters,
                )
                # determine which polygons to keep to meet concentration c
                if !isempty(floe_poly_list)
                    nfloes = length(floe_poly_list)
                    floe_idx = shuffle(rng, range(1, nfloes))
                    floes_area = FT(0.0)
                    # keep adding polygons as floes until concentration is met
                    while !isempty(floe_idx) && floes_area/open_area <= c
                        idx = pop!(floe_idx)
                        poly_pieces_list = intersect_polys(floe_poly_list[idx], open_cell_mpoly, FT)
                        for piece in poly_pieces_list
                            n_new_floes = _poly_to_floes!(
                                FT,
                                floe_arr,
                                piece,
                                generator.hmean,
                                generator.Δh,
                                domain.east.val - domain.west.val;
                                floe_settings,
                                rng,
                            )
                            floes_area += sum(Iterators.drop(floe_arr.area, nfloes_added))
                            nfloes_added += n_new_floes
                        end
                    end
                end
            end
        end
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    calculate_bonds!(floe_arr)
end

function Base.show(io::IO, floes::StructArray{Floe{FT}}; digits = 5) where FT
    overall_summary = "Floe List:"
    total_number_summary = "Number of floes: $(length(floes))"
    total_area_summary = "Total floe area: $(round(sum(floes.area); digits))"
    avg_height_summary = "Average floe height: $(round(mean(floes.height); digits))"
    print(io, overall_summary, "\n",
        "  ⊢", total_number_summary, "\n",
        "  ⊢", total_area_summary, "\n",
        "  ∟", avg_height_summary)
end