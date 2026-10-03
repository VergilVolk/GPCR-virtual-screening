#!/usr/bin/env bash
# Stage-3 runtime compatibility repair.
#   Frozen spec   : project/environment_drugclip_cpu.yml   (NOT modified)
#   Override      : MKL = 2024.0.x   (single authorized compatibility override)
#   Unchanged     : pytorch = 2.0.1, all other frozen pins
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
MM=/root/miniforge3/bin/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
YML=$SRC/project/environment_drugclip_cpu.yml
ENVPY=$MAMBA_ROOT_PREFIX/envs/drugclip-cpu/bin/python

echo "=== frozen environment spec (unchanged) ==="
sha256sum "$YML"
echo "=== state before override ==="
"$MM" list -n drugclip-cpu 2>/dev/null | grep -Ei '^(pytorch|mkl|intel-openmp|_openmp_mutex|llvm-openmp|python) ' || true
"$ENVPY" -c "import torch" 2>&1 | tail -1

echo "=== applying the single authorized override: mkl=2024.0.* ==="
"$MM" install -y -n drugclip-cpu -c conda-forge --channel-priority=flexible 'mkl=2024.0.*' 2>&1 | tail -8

echo "=== state after override ==="
"$MM" list -n drugclip-cpu 2>/dev/null | grep -Ei '^(pytorch|mkl|intel-openmp|_openmp_mutex|llvm-openmp|python) ' || true

echo "=== torch import test ==="
"$ENVPY" - <<'PY'
import sys
try:
    import torch
    print("TORCH_IMPORT_OK", torch.__version__, "cuda", torch.cuda.is_available())
except Exception as e:
    print("TORCH_IMPORT_FAIL", type(e).__name__, e)
    sys.exit(9)
PY
echo "MKL_OVERRIDE_DONE"