#!/bin/bash
# Generate the input and run all benchmarks. Run from the work directory made by setup.sh.
# Results go to results/.
set -e
DIR=$(dirname "$(realpath "$0")")
[ -f floes_10k.jld2 ] || julia -t 8 --project=env-branch "$DIR/make_floes.jl" 10000 100 floes_10k.jld2
export BENCH_SIZES=1000,10000,100000 BENCH_SECONDS=15
mkdir -p results
julia -t 1  --project=env-main "$DIR/bench_main.jl" > results/main_t1.txt 2>&1
julia -t 20 --project=env-main "$DIR/bench_main.jl" > results/main_t20.txt 2>&1
BENCH_BACKENDS=CPU julia -t 1 --project=env-branch "$DIR/bench_branch.jl" > results/branch_t1.txt 2>&1
julia -t 20 --project=env-branch "$DIR/bench_branch.jl" > results/branch_t20.txt 2>&1
