#!/bin/bash
# Create the work directory: a worktree of main and one Julia environment per Subzero
# version. Usage: setup.sh <workdir> [main-ref]
set -e
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
WORK=$1; REF=${2:-main}
mkdir -p "$WORK" && cd "$WORK"
git -C "$REPO" worktree add --detach "$PWD/subzero-main" "$REF"
julia --project=env-branch -e "using Pkg; Pkg.develop(path=\"$REPO\"); Pkg.add([\"CUDA\", \"BenchmarkTools\", \"JLD2\", \"StructArrays\", \"KernelAbstractions\", \"Adapt\", \"GeoInterface\"]); Pkg.precompile()"
julia --project=env-main -e "using Pkg; Pkg.develop(path=\"$PWD/subzero-main\"); Pkg.add([\"BenchmarkTools\", \"JLD2\", \"StructArrays\", \"GeoInterface\"]); Pkg.precompile()"
