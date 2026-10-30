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
    # The CPUs this process may use (e.g. those of a Slurm job), one per line
    if [ -f /proc/self/status ]; then
        ALLOWED=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status |
            tr , '\n' | awk -F- '{for (c = $1; c <= ($2 == "" ? $1 : $2); c++) print c}')
        SYS=/sys/devices/system/cpu
        # Other CPUs with cores of different speeds (e.g. AMD Zen 5 + Zen 5c, ARM
        # big.LITTLE) aren't detected: warn if their maximum frequencies differ by more than 20%
        if [ -z "$BENCH_CPUS" ] && [ -f $SYS/cpu0/cpufreq/cpuinfo_max_freq ]; then
            for c in $ALLOWED; do cat $SYS/cpu$c/cpufreq/cpuinfo_max_freq; done | sort -n |
                awk 'NR == 1 {min = $1} {max = $1}
                END {if (min < 0.8 * max) print "Warning: the CPUs have different maximum frequencies (" min " to " max " kHz). Set BENCH_CPUS to run on the fast ones only."}'
        fi
        # Threads that aren't pinned move between sockets, and so does their memory
        if [ -f $SYS/cpu0/topology/physical_package_id ]; then
            SOCKETS=$(for c in $ALLOWED; do cat $SYS/cpu$c/topology/physical_package_id; done | sort -u | wc -l)
            if [ "$SOCKETS" -gt 1 ]; then
                echo "Warning: the CPUs are on $SOCKETS sockets. Set BENCH_CPUS to the CPUs of one socket (see lscpu) for less noise."
            fi
        fi
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
START=$SECONDS
minutes() { echo "$(( ($1 + 30) / 60 )) min"; }
# Run a command with its output, start time and duration in the log file $1, e.g. to choose
# the time limit of a Slurm job
logged() {
    local log=$1 start=$SECONDS
    shift
    {
        echo "Started $(date '+%F %T')"
        "$@" 2>&1
        echo "Finished $(date '+%F %T') after $(minutes $((SECONDS - start)))"
    } | tee "$log"
}
mkdir -p logs
# The input is made with 1 thread and the baseline, so it doesn't depend on the number of
# threads or on the commit being benchmarked (see make_input.jl)
[ -f "${BENCH_FLOES:-floes_10k.jld2}" ] || logged logs/make_input.txt $PIN julia -t 1 --project=env-baseline "$DIR/make_input.jl" 10000 100 floes_10k.jld2
for t in 1 "$THREADS"; do
    logged "logs/baseline_t$t.txt" $PIN julia -t "$t" --project=env-baseline "$DIR/bench.jl" baseline
done
logged logs/current_t1.txt env BENCH_BACKENDS=CPU $PIN julia -t 1 --project=env-current "$DIR/bench.jl" current
logged "logs/current_t$THREADS.txt" $PIN julia -t "$THREADS" --project=env-current "$DIR/bench.jl" current
julia "$DIR/summarize.jl"
echo "All benchmarks took $(minutes $((SECONDS - START)))"
