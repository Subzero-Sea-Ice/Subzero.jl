#!/bin/bash
# Run the benchmarks on Snellius (SURF) on a quarter of a gpu_a100 node: 18 cores, one
# A100 and 120 GiB of memory, the smallest unit there (128 SBU per hour). Submit from the
# repository, checked out at the commit to benchmark:
#
#     sbatch benchmark/snellius.sh
#
# On a quarter of an H100 node (16 cores and one H100):
#
#     sbatch --partition=gpu_h100 --cpus-per-task=16 benchmark/snellius.sh
#
# The first job makes the work directory with setup.sh: the baseline worktree and the Julia
# environments, with CUDA.jl. That needs internet access. If the compute nodes don't have
# it, run this on a login node first (it doesn't need a GPU):
#
#     BENCH_GPU=CUDA benchmark/setup.sh ~/subzero-benchmark
#
# To compare with the results of other machines, copy the input they used (floes_10k.jld2)
# into the work directory first. Otherwise the first job makes a new one.
#
# The results go to benchmark/results/snellius-<partition>-<cores>c/ in the repository.
#
#SBATCH --job-name=subzero-benchmark
#SBATCH --partition=gpu_a100
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=04:00:00
#SBATCH --output=subzero-benchmark-%j.out
set -eo pipefail

# Julia from juliaup: curl -fsSL https://install.julialang.org | sh
export PATH="$HOME/.juliaup/bin:$PATH"
command -v julia > /dev/null || { echo "julia not found: install juliaup or load a Julia module"; exit 1; }

REPO=$(git -C "${SLURM_SUBMIT_DIR:-.}" rev-parse --show-toplevel)
WORK=${WORK:-$HOME/subzero-benchmark}
export BENCH_MACHINE=${BENCH_MACHINE:-snellius-${SLURM_JOB_PARTITION}-${SLURM_CPUS_PER_TASK}c}

echo "Job $SLURM_JOB_ID on $(hostname), partition $SLURM_JOB_PARTITION, $SLURM_CPUS_PER_TASK CPUs"
echo "Subzero: $REPO at $(git -C "$REPO" rev-parse --short HEAD), work directory: $WORK"
lscpu | grep -E "^(Model name|Socket|Core|Thread|NUMA node)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv

if [ ! -f "$WORK/setup-done" ]; then
    [ -d "$WORK/env-current" ] && { echo "$WORK exists but its setup didn't finish: remove it, run git worktree prune and resubmit"; exit 1; }
    BENCH_GPU=CUDA "$REPO/benchmark/setup.sh" "$WORK"
fi
cd "$WORK"
"$REPO/benchmark/run.sh"
