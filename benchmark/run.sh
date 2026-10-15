#!/bin/bash
# Generate the input and run the benchmarks. Run from the work directory made by setup.sh.
# Results go to results/ next to this script, the output of each run to logs/.
# BENCH_STEPS sets the steps to run (default: all). The baseline doesn't change, so its
# results are only made if they don't exist yet (BENCH_FORCE=1 reruns them).
#
# BENCH_CPUS sets the CPUs that Julia runs on (a taskset list, e.g. 0-11). By default, on
# CPUs with fast and slow cores (Intel P-cores and E-cores), only the fast cores, and
# otherwise all CPUs. BENCH_THREADS sets the number of threads (default: the number of
# CPUs in BENCH_CPUS). The 1-thread runs use the same CPUs, so they run on a fast core.
set -eo pipefail
DIR=$(dirname "$(realpath "$0")")
if [ -z "${BENCH_CPUS+x}" ]; then
    BENCH_CPUS=""
    [ -f /sys/devices/cpu_core/cpus ] && BENCH_CPUS=$(cat /sys/devices/cpu_core/cpus)
fi
export BENCH_CPUS
PIN=""
[ -n "$BENCH_CPUS" ] && PIN="taskset -c $BENCH_CPUS"
THREADS=${BENCH_THREADS:-$($PIN nproc)}
echo "CPUs: ${BENCH_CPUS:-all}, threads: 1 and $THREADS"
[ -f "${BENCH_FLOES:-floes_10k.jld2}" ] || $PIN julia -t "$THREADS" --project=env-current "$DIR/make_input.jl" 10000 100 floes_10k.jld2
mkdir -p logs
for t in 1 "$THREADS"; do
    $PIN julia -t "$t" --project=env-baseline "$DIR/bench.jl" baseline 2>&1 | tee "logs/baseline_t$t.txt"
done
BENCH_BACKENDS=CPU $PIN julia -t 1 --project=env-current "$DIR/bench.jl" current 2>&1 | tee logs/current_t1.txt
$PIN julia -t "$THREADS" --project=env-current "$DIR/bench.jl" current 2>&1 | tee "logs/current_t$THREADS.txt"
julia "$DIR/summarize.jl"
