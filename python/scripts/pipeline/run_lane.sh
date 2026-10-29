#!/usr/bin/env bash
# Run whole Houses one after another, 2 chains each (flip + shift), then diagnostics.
#   run_lane.sh 119 117 115 ...
cd /home/josephyu90/circular_factor_model_remake/python || exit 1
for c in "$@"; do
  if [ -f "results/full/c2/H$c/summary.json" ] && [ -f "results/full/c1/H$c/summary.json" ]; then
    echo "$(date '+%F %T') H$c already done, skipping"; continue
  fi
  echo "$(date '+%F %T') H$c start"
  for s in 1 2; do
    uv run python scripts/run_congress.py "$c" --bill-flip --bill-shift \
      --burnin 20000 --samples 40000 --seed $((2020 + s)) \
      --out "results/full/c$s" > "logs/full_H${c}_c$s.log" 2>&1 &
  done
  wait
  for s in 1 2; do
    grep -qE "Traceback|Error" "logs/full_H${c}_c$s.log" && echo "$(date '+%F %T') H$c chain $s FAILED (see logs/full_H${c}_c$s.log)"
  done
  uv run python scripts/chain_diagnostics.py results/full/c*/H"$c" > "results/full/diag_H$c.txt" 2>&1
  echo "$(date '+%F %T') H$c done: $(grep -E 'chi_0 split' results/full/diag_H$c.txt)"
done
echo "$(date '+%F %T') lane finished"
