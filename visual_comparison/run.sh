#!/bin/bash
# Run the examples with the original and the current implementation and render side-by-side
# videos. Run from the work directory made by setup.sh. Output: videos/<example>.mp4, the
# saved runs in runs/, and the output of each run in logs/.
#
# VIS_STEPS: number of timesteps (default 1000). VIS_EVERY: timesteps per frame (default 10).
# VIS_EXAMPLES: examples to run (default: shear_flow simple_strait forcing_contained_floes
# moving_bounds).
# VIS_BACKEND=CPU: run the current version on the CPU backend, even if there is a GPU.
# The floes of each example are made once (floes/<example>.jld2) and reused by later runs;
# delete the file to make new ones. Both versions run on 1 thread, so that the original's
# multi-threaded collisions don't make its runs differ from each other.
set -eo pipefail
DIR=$(dirname "$(realpath "$0")")
REPO=$(git -C "$DIR" rev-parse --show-toplevel)
[ -f setup-done ] || { echo "Run setup.sh first and run this from its work directory"; exit 1; }
STEPS=${VIS_STEPS:-1000}
EVERY=${VIS_EVERY:-10}
EXAMPLES=${VIS_EXAMPLES:-"shear_flow simple_strait forcing_contained_floes moving_bounds"}
ORIGINAL_COMMIT=$(git -C subzero-original rev-parse --short HEAD)
CURRENT_COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
[ -n "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] && CURRENT_COMMIT="$CURRENT_COMMIT-dirty"
mkdir -p floes runs/original runs/current logs
for ex in $EXAMPLES; do
    [ -f "floes/$ex.jld2" ] || julia -t 1 --project=env-original "$DIR/make_floes.jl" "$ex" "floes/$ex.jld2"
    VIS_COMMIT=$ORIGINAL_COMMIT julia -t 1 --project=env-original "$DIR/run_example.jl" \
        "$ex" original "$STEPS" "$EVERY" "floes/$ex.jld2" "runs/original/$ex.jld2" 2>&1 | tee "logs/${ex}_original.txt"
    VIS_COMMIT=$CURRENT_COMMIT julia -t 1 --project=env-current "$DIR/run_example.jl" \
        "$ex" current "$STEPS" "$EVERY" "floes/$ex.jld2" "runs/current/$ex.jld2" 2>&1 | tee "logs/${ex}_current.txt"
done
# shellcheck disable=SC2086
julia --project=env-current "$DIR/render.jl" runs videos $EXAMPLES 2>&1 | tee logs/render.txt
