#!/usr/bin/env bash
# Extend the best of all chains (main + brute) by 40k draws, once.  pick_lane.sh 102 103 ...
cd "${ROOT:-/workspace/Circular_Factor_Model/python}" || exit 1; export PATH=$HOME/.local/bin:$PATH
D=$(dirname "$(readlink -f "$0")")
for c in "$@"; do
  out=$(uv run python "$D/pick.py" "$c" 2>&1); best=$(echo "$out" | head -1); k=${best##*/}
  echo "$(date '+%F %T') H$c pick:"; echo "$out" | tail -n +2
  [ -f "results/full_ext/$k/H$c/summary.json" ] && { echo "  already extended"; continue; }
  uv run python scripts/run_congress.py "$c" --resume "$best" --samples 40000 --out "results/full_ext/$k" > "logs/ext_H${c}_$k.log" 2>&1
  echo "$(date '+%F %T') H$c $best extended: $(tail -n1 logs/ext_H${c}_$k.log)"
done
echo "$(date '+%F %T') pick lane finished"
