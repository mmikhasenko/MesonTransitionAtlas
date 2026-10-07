#!/bin/sh
# Full precompute, from the root of this repository:
#
#   sh run.sh [N]        N worker processes for the strong pass (default 3)
#
# 1. Solve all spectra once (cache/states.jls).
# 2. N workers fill cache/strong/, one file per decaying parent.
# 3. A final run computes the radiative and annihilation widths, reads the
#    strong cache, and writes data/transitions.json and site/data.js.
#
# Each Julia process needs 1-2 GB; on an 8 GB machine keep N ≤ 3.
set -e
N=${1:-3}
LOG=cache/logs
mkdir -p "$LOG"
julia --project=. -e 'using Pkg; Pkg.instantiate()'
[ -f cache/states.jls ] || SOLVE_ONLY=1 julia --project=. compute.jl > "$LOG/solve.log" 2>&1
i=0
while [ "$i" -lt "$N" ]; do
  STRONG_WORKER="$i/$N" julia --project=. compute.jl > "$LOG/worker-$i.log" 2>&1 &
  i=$((i + 1))
done
wait
julia --project=. compute.jl > "$LOG/assemble.log" 2>&1
tail -1 "$LOG/assemble.log"
