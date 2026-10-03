#!/usr/bin/env bash
# Stage-3A environment completion: provide the Intel OpenMP runtime that the
# pytorch 2.0.1 cpu build links against (undefined symbol iJIT_NotifyEvent),
# then retry the frozen Uni-Core install step verbatim.
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
MM=/root/miniforge3/bin/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
UNICORE=$SRC/project/tools/Uni-Core
ENVPY=$MAMBA_ROOT_PREFIX/envs/drugclip-cpu/bin/python

echo "=== [1/3] add intel-openmp to complete the frozen env's runtime deps ==="
"$MM" install -y -n drugclip-cpu -c conda-forge --channel-priority=flexible intel-openmp || {
  echo "INTEL_OPENMP_INSTALL_FAILED"; exit 8; }

echo "=== [2/3] retry Uni-Core install (frozen step verbatim) ==="
pushd "$UNICORE" >/dev/null
"$ENVPY" setup.py install --disable-cuda-ext || { popd >/dev/null; echo "UNICORE_INSTALL_FAILED"; exit 7; }
popd >/dev/null

echo "=== [3/3] environment identity ==="
"$ENVPY" - <<'PY'
import json, sys
info = {"python": sys.version.split()[0]}
for m in ("rdkit", "torch", "unicore", "lmdb", "numpy", "scipy", "pandas", "sklearn", "biopandas"):
    try:
        mod = __import__(m); info[m] = getattr(mod, "__version__", "ok")
    except Exception as e:
        info[m] = f"FAIL {type(e).__name__}: {e}"
try:
    import torch; info["torch_cuda"] = torch.cuda.is_available()
except Exception: pass
print(json.dumps(info, indent=2))
PY
echo "STAGE3A_ENV_READY"
