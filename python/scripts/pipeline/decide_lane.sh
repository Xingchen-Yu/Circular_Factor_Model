#!/usr/bin/env bash
# After the main 2-chain runs: per House, extend the better chain (+40k, once), or run 4 extra seeds, or nothing.
#   ROOT=<python dir> decide_lane.sh 100 101 ...
cd "${ROOT:-/workspace/Circular_Factor_Model/python}" || exit 1; export PATH=$HOME/.local/bin:$PATH
D=$(dirname "$(readlink -f "$0")")
for c in "$@"; do
  out=$(uv run python "$D/decide.py" "$c" 2>&1); act=$(echo "$out" | head -1)
  echo "$(date '+%F %T') $(echo "$out" | tail -1) -> $act"
  case "$act" in
    extend*) k=${act#extend }
      [ -f "results/full_ext/c$k/H$c/summary.json" ] && { echo "  already extended"; continue; }
      uv run python scripts/run_congress.py "$c" --resume "results/full/c$k" --samples 40000 \
        --out "results/full_ext/c$k" > "logs/ext_H${c}_c$k.log" 2>&1
      echo "$(date '+%F %T') H$c c$k extended: $(tail -n1 logs/ext_H${c}_c$k.log)";;
    brute)
      for s in 5 6 7 8; do
        [ -f "results/brute/c$s/H$c/summary.json" ] && continue
        uv run python scripts/run_congress.py "$c" --bill-flip --bill-shift --burnin 20000 --samples 10000 \
          --seed $((2020 + s)) --out "results/brute/c$s" > "logs/brute_H${c}_c$s.log" 2>&1 &
      done; wait
      echo "$(date '+%F %T') H$c brute done: $(tail -qn1 logs/brute_H${c}_c*.log | grep -o 'chi_0 = [0-9.]*' | tr '\n' ' ')";;
  esac
done
echo "$(date '+%F %T') decide lane finished"
