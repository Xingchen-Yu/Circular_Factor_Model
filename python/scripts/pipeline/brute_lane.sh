#!/usr/bin/env bash
# Brute-force fix: 4 extra chains per House (seeds 2025-2028, 20k burn-in + 10k draws), one House at a time.
#   brute_lane.sh 100 101 ...
export PATH=$HOME/.local/bin:$PATH; cd /workspace/Circular_Factor_Model/python || exit 1
for c in "$@"; do
  echo "$(date '+%F %T') brute H$c start"
  uv run python -c "from circfactor.data import load_congress; load_congress($c, \"House\")" > /dev/null || echo "H$c download FAILED"
  for s in 5 6 7 8; do
    [ -f "results/brute/c$s/H$c/summary.json" ] && continue
    uv run python scripts/run_congress.py "$c" --bill-flip --bill-shift \
      --burnin 20000 --samples 10000 --seed $((2020 + s)) \
      --out "results/brute/c$s" > "logs/brute_H${c}_c$s.log" 2>&1 &
  done
  wait
  echo "$(date '+%F %T') brute H$c done: $(tail -qn1 logs/brute_H${c}_c*.log | grep -o 'chi_0 = [0-9.]*' | tr '\n' ' ')"
done
echo "$(date '+%F %T') brute lane finished"
