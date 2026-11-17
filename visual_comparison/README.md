# Visual comparison with the original CPU implementation

Side-by-side videos of four documentation examples, run with the original CPU implementation
(`main` at ae96ff2, before the GPU work) and with the current implementation (this checkout,
on a GPU if there is one). They're a quick check by eye that nothing looks broken.

| example | file in `docs/literate/examples/` |
|---|---|
| Shear Flow with Periodic Boundaries | `shear_flow.jl` |
| Simple Strait Simulation | `simple_strait.jl` |
| Floes Bounded by Ocean Currents | `forcing_contained_floes.jl` |
| Moving Boundaries | `moving_bounds.jl` |

## How to run

```sh
REPO=/path/to/Subzero.jl            # checked out at the commit to compare
WORK=$(mktemp -d)                   # work directory for environments, runs, and videos
$REPO/visual_comparison/setup.sh "$WORK"   # worktree of the original + env-original and env-current
cd "$WORK" && $REPO/visual_comparison/run.sh
git -C $REPO worktree remove "$WORK/subzero-original"   # clean up afterwards
```

`setup.sh` takes about 10 minutes (precompiling CairoMakie and the GPU package). It adds
CUDA.jl to the current environment if `nvidia-smi` is found and AMDGPU.jl if `rocminfo` or
`rocm-smi` is; choose with `VIS_GPU` (e.g. `VIS_GPU=AMDGPU`, or `""` for none). Its second
argument is the original commit (default ae96ff2).

`run.sh` takes a few minutes for 1000 timesteps. The videos are in `videos/<example>.mp4`.

| Variable | Meaning | Default |
|---|---|---|
| `VIS_STEPS` | number of timesteps | `1000` |
| `VIS_EVERY` | timesteps per frame | `10` |
| `VIS_EXAMPLES` | examples to run | all four |
| `VIS_BACKEND` | `CPU` runs the current version on the CPU backend, even with a GPU | GPU if there is one |

## How it works

- `examples.jl` has the setup of each example, as in the documentation, written so that it
  works with both versions of Subzero.
- The examples make their floes with a Voronoi tessellation, which isn't reproducible with a
  fixed rng. So `make_floes.jl` makes them once with the original version and saves them
  (`floes/<example>.jld2`), and both versions start from these floes. Later runs reuse them;
  delete the file to make new ones.
- `run_example.jl` runs an example with one version, on 1 thread (the original's
  multi-threaded collisions aren't reproducible), and saves every `VIS_EVERY` timesteps the
  floes' outlines (with their periodic ghosts, as the documentation's animations show them),
  ids, and centroids, as plain arrays that both versions can read (`runs/`).
- `render.jl` draws three panels: the original, the current version, and both on top of each
  other (original outlines orange and wide, current outlines cyan and narrow). Moving
  boundaries are drawn as grey bands where the walls have moved in. Below them, two
  graphs of the differences in position against time, on a log and on a linear scale, draw
  themselves as the video plays.

## How to read the videos

The line above the panels gives the number of floes in each run and the largest distance
between the centroids of the same floe in both runs. That distance is only over the initial
floes that are still in both runs: fracture and ridging make new floes with new ids, and
once these events differ between the runs, the same new id isn't the same floe.

The graphs below the panels show the largest and the median of these distances against
time. The log scale shows the early growth from tiny differences; the linear scale shows
how many metres the later growth adds, which looks flat on the log scale. A dashed line
marks when the runs start to have different floes.

The runs don't stay identical. The current implementation rounds slightly differently (the
whole-simulation reference test allows a relative difference of 1e-5 after 100 timesteps), and
collisions make such differences grow. With fracture and ridging, a small difference can tip a
discrete decision: in the Simple Strait, both runs fracture the same floe the same way at
timestep 75, but by timestep 100 ridging has gone differently for the new pieces, and from
timestep 160 the runs have different floes. On the CPU backend, the current version gives
exactly the same results as on the GPU, so this comes from the implementation, not the GPU.

So late differences are expected. Look for things that only happen in one run: floes that
jump, vanish, overlap, or stop moving, or differences that appear early, before collisions
could have made them grow.
