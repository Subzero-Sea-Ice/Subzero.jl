#!/bin/bash
# Create the work directory for the visual comparison: a worktree of the original CPU
# implementation and one Julia environment per version. Usage: setup.sh <workdir> [original-ref]
# The original is main before the GPU work (ae96ff2) by default. The current version is this
# checkout of the repository. The current environment gets CUDA if nvidia-smi is found and
# AMDGPU if rocminfo or rocm-smi is found; set VIS_GPU (e.g. CUDA, AMDGPU, or "") to choose.
set -e
# Precompile once at the end instead of after each Pkg.add
export JULIA_PKG_PRECOMPILE_AUTO=0
if [ -z "${VIS_GPU+x}" ]; then
    VIS_GPU=""
    command -v nvidia-smi > /dev/null && VIS_GPU="CUDA"
    { command -v rocminfo || command -v rocm-smi; } > /dev/null && VIS_GPU="$VIS_GPU AMDGPU"
fi
PKGS=$(for p in CairoMakie $VIS_GPU; do printf '"%s", ' "$p"; done)
DIR=$(dirname "$(realpath "$0")")
REPO=$(git -C "$DIR" rev-parse --show-toplevel)
WORK=$1; REF=${2:-ae96ff2}
mkdir -p "$WORK/env-original" "$WORK/env-current" && cd "$WORK"
git -C "$REPO" worktree add --detach "$PWD/subzero-original" "$REF"
cp "$DIR/Project.toml" env-original/
cp "$DIR/Project.toml" env-current/
julia --project=env-current -e "using Pkg; Pkg.develop(path=\"$REPO\"); Pkg.add([$PKGS]); Pkg.instantiate(); Pkg.precompile()"
julia --project=env-original -e "using Pkg; Pkg.develop(path=\"$PWD/subzero-original\"); Pkg.instantiate(); Pkg.precompile()"
# Written last, so a work directory without it is incomplete
touch setup-done
