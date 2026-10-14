#!/bin/bash
# Create the work directory: a worktree of main and one Julia environment per Subzero
# version. Usage: setup.sh <workdir> [main-ref]
# The GPU package for the branch environment is CUDA if nvidia-smi is found, AMDGPU if
# rocminfo or rocm-smi is found. Set BENCH_GPU (e.g. CUDA, AMDGPU, or "CUDA AMDGPU") to choose.
set -e
if [ -z "${BENCH_GPU+x}" ]; then
    BENCH_GPU=""
    command -v nvidia-smi > /dev/null && BENCH_GPU="CUDA"
    { command -v rocminfo || command -v rocm-smi; } > /dev/null && BENCH_GPU="$BENCH_GPU AMDGPU"
fi
GPU_PKGS=$(for p in $BENCH_GPU; do printf '"%s", ' "$p"; done)
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
WORK=$1; REF=${2:-main}
mkdir -p "$WORK" && cd "$WORK"
git -C "$REPO" worktree add --detach "$PWD/subzero-main" "$REF"
julia --project=env-branch -e "using Pkg; Pkg.develop(path=\"$REPO\"); Pkg.add([$GPU_PKGS\"BenchmarkTools\", \"JLD2\", \"StructArrays\", \"KernelAbstractions\", \"Adapt\", \"GeoInterface\"]); Pkg.precompile()"
julia --project=env-main -e "using Pkg; Pkg.develop(path=\"$PWD/subzero-main\"); Pkg.add([\"BenchmarkTools\", \"JLD2\", \"StructArrays\", \"GeoInterface\"]); Pkg.precompile()"
