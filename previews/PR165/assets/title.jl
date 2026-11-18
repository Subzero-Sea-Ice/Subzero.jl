using Subzero, JLD2, CairoMakie, GeoInterfaceMakie, Random, Statistics

const FT = Float64
const Lx = 8e4
const Ly = 8e4
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20
const nΔt = 7500  # 150 frames, a 7.5 s loop at 20 fps

# Model instantiation
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

uvels = repeat(
    [range(0, 0.2, length = 21); range(0.2, 0, length = 20)],
    outer = (1, 41),
)
ocean = Ocean(;
    u = uvels',
    grid,
    v = 0.0,
    temp = 0.0,
)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)

# Domain creation
nboundary = PeriodicBoundary(North; grid)
sboundary = PeriodicBoundary(South; grid)
eboundary = PeriodicBoundary(East; grid)
wboundary = PeriodicBoundary(West; grid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe creation
floe_settings = FloeSettings()
floe_arr = initialize_floe_field(FT, 400, [0.85], domain, hmean, Δh; rng = Xoshiro(1), floe_settings)

# Model creation
model = Model(; grid, ocean, atmos, domain, floes = floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
dir = "output/title/"

# Output setup
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(floewriter)

simulation = Simulation(; model, consts, Δt, nΔt, writers, floe_settings, verbose = true, rng = Xoshiro(1))
run!(simulation)

function plot_logo(floe_fn, Lx, Ly, dir)
    # Floe Information
    file = jldopen(floe_fn)
    sim_polys = file["poly"]
    timesteps = keys(sim_polys)

    # Set up observables for recording (updated whenever `time[]` is set to a new value)
    time = Observable(timesteps[1])  # note these are strings as we index into a JLD2 file
    time_polys = @lift(sim_polys[$time])

    # Set up figure
    fig = Figure(; size = (800, 200))

    floe_rgb = (255, 255, 255)
    floe_color = RGBf((floe_rgb ./ 255)...)
    ocean_rgb = (0, 157, 196)
    ocean_color = RGBf((ocean_rgb ./ 255)...)

    ax = Axis(fig[1, 1]; limits = (0.0, Lx, 0.0, Ly), backgroundcolor = ocean_color)
    hidedecorations!(ax)

    # Plot starting state (floes + topography)
    poly!(time_polys, color = floe_color, strokecolor = :black, strokewidth = 0.5)  # floes
    text!(
        Lx/2, Ly/2;
        align = (:center, :center),
        text = "Subzero.jl",
        font = :bold,
        fontsize = 130,
        strokecolor = :white,
        strokewidth = 5,
    )

    # Keep the GIF small: render at 800 x 200 and snap every pixel to a 9-colour palette.
    # The palette holds white, black, the ocean blue and two blends between each pair,
    # so the anti-aliased edges stay smooth.
    base = [floe_rgb, (0, 0, 0), ocean_rgb]
    blend(a, b, f) = round.(Int, a .+ f .* (b .- a))
    palette = copy(base)
    for (i, a) in enumerate(base), b in base[(i + 1):end], f in (1/3, 2/3)
        push!(palette, blend(a, b, f))
    end
    nearest(c) = argmin(p -> sum(abs2, p .- c), palette)
    to_rgb(c) = (c = convert(RGBf, c); round.(Int, 255 .* (c.r, c.g, c.b)))

    # ffmpeg (bundled with Makie) needs the palette as an image of exactly 256 pixels
    palette_fn = dir*"palette.ppm"
    open(palette_fn, "w") do io
        padded = [palette; fill(palette[end], 256 - length(palette))]
        write(io, "P6\n16 16\n255\n", UInt8.(collect(Iterators.flatten(padded))))
    end

    # Create movie: pipe the quantized frames as raw RGB into ffmpeg
    CairoMakie.activate!(; px_per_unit = 1)
    w, h = 800, 200
    ffmpeg = CairoMakie.Makie.FFMPEG_jll.ffmpeg()
    cmd = `$ffmpeg -loglevel error -y -f rawvideo -pix_fmt rgb24 -s $(w)x$(h) -framerate 20
        -i pipe:0 -i $palette_fn -lavfi "paletteuse=dither=none:diff_mode=rectangle"
        -loop 0 $(dir*"title.gif")`
    open(cmd, "w") do io
        for t in timesteps
            time[] = t  # updating the time auto updates the related time_polys
            img = colorbuffer(fig)
            @assert size(img) == (h, w)
            frame = permutedims(img)  # ffmpeg expects the pixels row by row
            write(io, UInt8.(collect(Iterators.flatten(nearest.(to_rgb.(frame))))))
        end
    end

    # Close file
    close(file)
    rm(palette_fn)
end

plot_logo(dir*"floes.jld2", Lx, Ly, dir)