#!/usr/bin/env bash
# Stage-3A environment build (final).
# The frozen YAML is used UNCHANGED. Only the solver's channel priority differs from
# the frozen setup script's plain invocation, because under the default priority the
# frozen spec is unsolvable on today's channels ("pytorch=2.0.1 ... not installable").
# Steps 2-3 below are verbatim the frozen script's remaining steps.
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
MM=/root/miniforge3/bin/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
YML=$SRC/project/environment_drugclip_cpu.yml
DRUGCLIP=$SRC/project/tools/DrugCLIP
UNICORE=$SRC/project/tools/Uni-Core
ENVPY=$MAMBA_ROOT_PREFIX/envs/drugclip-cpu/bin/python

echo "=== frozen YAML sha256 (unchanged) ==="
sha256sum "$YML"

if [ ! -x "$ENVPY" ]; then
  echo "=== [1/3] create drugclip-cpu from the frozen YAML ==="
  "$MM" create -y -n drugclip-cpu -f "$YML" --channel-priority=flexible || {
    echo "ENV_CREATE_FAILED"; exit 5; }
fi
"$ENVPY" -V || { echo "ENV_PYTHON_MISSING"; exit 5; }

echo "=== [2/3] audited CPU device-placement patch (frozen script verbatim) ==="
"$ENVPY" "$SRC/project/scripts/prepare_drugclip_cpu_checkout.py" \
  --checkout "$DRUGCLIP" --audit "$DRUGCLIP/artifacts/cpu_patch_audit.json" || {
    echo "CPU_PATCH_FAILED"; exit 6; }

echo "=== [3/3] install Uni-Core --disable-cuda-ext (frozen script verbatim) ==="
pushd "$UNICORE" >/dev/null
"$ENVPY" setup.py install --disable-cuda-ext || { popd >/dev/null; echo "UNICORE_INSTALL_FAILED"; exit 7; }
popd >/dev/null

echo "=== environment identity ==="
"$ENVPY" - <<'PY'
import json, sys, rdkit, torch
info = {"python": sys.version.split()[0], "rdkit": rdkit.__version__,
        "torch": torch.__version__, "cuda": torch.cuda.is_available()}
try:
    import unicore; info["unicore"] = "ok"
except Exception as e:
    info["unicore"] = f"FAIL {type(e).__name__}: {e}"
for m in ("lmdb", "numpy", "scipy", "pandas", "sklearn"):
    try:
        mod = __import__(m); info[m] = getattr(mod, "__version__", "ok")
    except Exception as e:
        info[m] = f"FAIL {type(e).__name__}"
print(json.dumps(info, indent=2))
PY
echo "STAGE3A_ENV_READY"
