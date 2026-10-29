#!/usr/bin/env bash
# One-time setup on a RunPod pod: uv, repo, environment, GPU check.
set -euo pipefail
cd /workspace
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
[ -d Circular_Factor_Model ] || git clone https://github.com/Xingchen-Yu/Circular_Factor_Model.git
cd Circular_Factor_Model/python
git pull --ff-only
uv sync
uv run python -c "import jax; print(jax.devices())"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
mkdir -p logs results
