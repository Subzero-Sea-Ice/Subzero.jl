#!/bin/bash
# Generate the input and run all benchmarks. Run from the work directory made by setup.sh.
# Results go to results/. BENCH_THREADS sets the number of threads (default: all CPUs).
# The main benchmarks are a baseline that doesn't change when the branch changes, so they
# only run if their results don't exist yet. Delete results/main_*.txt to rerun them.
set -e
DIR=$(dirname "$(realpath "$0")")
THREADS=${BENCH_THREADS:-$(nproc)}
[ -f floes_10k.jld2 ] || julia -t "$THREADS" --project=env-branch "$DIR/make_floes.jl" 10000 100 floes_10k.jld2
mkdir -p results
for t in 1 "$THREADS"; do
    [ -f "results/main_t$t.txt" ] || julia -t "$t" --project=env-main "$DIR/bench_main.jl" > "results/main_t$t.txt" 2>&1
done
BENCH_BACKENDS=CPU julia -t 1 --project=env-branch "$DIR/bench_branch.jl" > results/branch_t1.txt 2>&1
julia -t "$THREADS" --project=env-branch "$DIR/bench_branch.jl" > "results/branch_t$THREADS.txt" 2>&1
