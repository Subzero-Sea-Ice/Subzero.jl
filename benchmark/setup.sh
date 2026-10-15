#!/bin/bash
# Create the work directory: a worktree of the baseline and one Julia environment per
# Subzero version. Usage: setup.sh <workdir> [baseline-ref]
# The baseline is main before the GPU work (ae96ff2) by default.
# The GPU package for the current environment is CUDA if nvidia-smi is found, AMDGPU if
# rocminfo or rocm-smi is found. Set BENCH_GPU (e.g. CUDA, AMDGPU, or "CUDA AMDGPU") to choose.
set -e
# Precompile once at the end instead of after each Pkg.add. This also avoids a bug in
# Julia 1.13, where the precompile after adding CUDA fails with "Unknown package: cuFFT".
export JULIA_PKG_PRECOMPILE_AUTO=0
if [ -z "${BENCH_GPU+x}" ]; then
    BENCH_GPU=""
    command -v nvidia-smi > /dev/null && BENCH_GPU="CUDA"
    { command -v rocminfo || command -v rocm-smi; } > /dev/null && BENCH_GPU="$BENCH_GPU AMDGPU"
fi
GPU_PKGS=$(for p in $BENCH_GPU; do printf '"%s", ' "$p"; done)
ADD_GPU=""
[ -n "$GPU_PKGS" ] && ADD_GPU="Pkg.add([$GPU_PKGS]);"
DIR=$(dirname "$(realpath "$0")")
REPO=$(git -C "$DIR" rev-parse --show-toplevel)
WORK=$1; REF=${2:-ae96ff2}
mkdir -p "$WORK/env-baseline" "$WORK/env-current" && cd "$WORK"
git -C "$REPO" worktree add --detach "$PWD/subzero-baseline" "$REF"
cp "$DIR/Project.toml" env-baseline/
cp "$DIR/Project.toml" env-current/
julia --project=env-current -e "using Pkg; Pkg.develop(path=\"$REPO\"); $ADD_GPU Pkg.instantiate(); Pkg.precompile()"
julia --project=env-baseline -e "using Pkg; Pkg.develop(path=\"$PWD/subzero-baseline\"); Pkg.instantiate(); Pkg.precompile()"
# Written last, so a work directory without it is incomplete (snellius.sh checks for it)
touch setup-done
