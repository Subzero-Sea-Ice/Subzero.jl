#!/bin/bash
# Generate the input and run the benchmarks. Run from the work directory made by setup.sh.
# Results go to results/$BENCH_MACHINE/ next to this script, the output of each run to logs/.
# BENCH_STEPS sets the steps to run (default: all). The baseline doesn't change, so its
# results are only made if they don't exist yet (BENCH_FORCE=1 reruns them).
#
# BENCH_MACHINE names the machine (default: the host name). Results of different machines
# can't be compared, so each machine has its own results folder.
#
# BENCH_CPUS sets the CPUs that Julia runs on (a taskset list, e.g. 0-11; Linux only). By
# default, on Intel CPUs with fast and slow cores (P-cores and E-cores), only the fast
# cores, and otherwise all CPUs. BENCH_THREADS sets the number of threads (default: the
# number of CPUs in BENCH_CPUS). The 1-thread runs use the same CPUs.
set -eo pipefail
DIR=$(dirname "$(realpath "$0")")
if [ -z "${BENCH_MACHINE:-}" ]; then
    BENCH_MACHINE=$(hostname -s)
    echo "BENCH_MACHINE isn't set, using the host name: $BENCH_MACHINE"
fi
export BENCH_MACHINE
if [ -z "${BENCH_CPUS+x}" ]; then
    BENCH_CPUS=""
    [ -f /sys/devices/cpu_core/cpus ] && BENCH_CPUS=$(cat /sys/devices/cpu_core/cpus)
    # Other CPUs with cores of different speeds (e.g. AMD Zen 5 + Zen 5c, ARM big.LITTLE)
    # aren't detected: warn if the maximum frequencies of the CPUs differ by more than 20%
    if [ -z "$BENCH_CPUS" ] && ls /sys/devices/system/cpu/cpu*/cpufreq/cpuinfo_max_freq > /dev/null 2>&1; then
        sort -n /sys/devices/system/cpu/cpu*/cpufreq/cpuinfo_max_freq | awk 'NR == 1 {min = $1} {max = $1}
            END {if (min < 0.8 * max) print "Warning: the CPUs have different maximum frequencies (" min " to " max " kHz). Set BENCH_CPUS to run on the fast ones only."}'
    fi
    # Threads that aren't pinned move between sockets, and so does their memory
    if [ -d /sys/devices/system/cpu/cpu0/topology ]; then
        SOCKETS=$(cat /sys/devices/system/cpu/cpu*/topology/physical_package_id | sort -u | wc -l)
        [ "$SOCKETS" -gt 1 ] && echo "Warning: $SOCKETS CPU sockets. Set BENCH_CPUS to the CPUs of one socket (see lscpu) for less noise."
    fi
fi
export BENCH_CPUS
PIN=""
if [ -n "$BENCH_CPUS" ]; then
    command -v taskset > /dev/null || { echo "BENCH_CPUS needs taskset (Linux)"; exit 1; }
    PIN="taskset -c $BENCH_CPUS"
fi
if command -v nproc > /dev/null; then
    NCPUS=$($PIN nproc)
else  # macOS
    NCPUS=$(sysctl -n hw.ncpu)
fi
THREADS=${BENCH_THREADS:-$NCPUS}
echo "Machine: $BENCH_MACHINE, CPUs: ${BENCH_CPUS:-all}, threads: 1 and $THREADS"
[ -f "${BENCH_FLOES:-floes_10k.jld2}" ] || $PIN julia -t "$THREADS" --project=env-current "$DIR/make_input.jl" 10000 100 floes_10k.jld2
mkdir -p logs
for t in 1 "$THREADS"; do
    $PIN julia -t "$t" --project=env-baseline "$DIR/bench.jl" baseline 2>&1 | tee "logs/baseline_t$t.txt"
done
BENCH_BACKENDS=CPU $PIN julia -t 1 --project=env-current "$DIR/bench.jl" current 2>&1 | tee logs/current_t1.txt
$PIN julia -t "$THREADS" --project=env-current "$DIR/bench.jl" current 2>&1 | tee "logs/current_t$THREADS.txt"
julia "$DIR/summarize.jl"
