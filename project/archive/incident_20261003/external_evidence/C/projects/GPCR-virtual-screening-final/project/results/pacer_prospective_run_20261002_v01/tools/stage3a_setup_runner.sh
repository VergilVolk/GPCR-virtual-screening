#!/usr/bin/env bash
# Stage-3A environment build (retry).
#   bzip2 is not installed in this WSL distro, so micromamba is provided by conda.
#   The frozen setup script is a CRLF checkout, which breaks 'set -o pipefail' in
#   WSL; a CRLF-stripped copy is used. Both hashes are recorded.
set -u
export MAMBA_ROOT_PREFIX=/root/micromamba
SRC=/mnt/c/projects/GPCR-virtual-screening-final
FROZEN=$SRC/project/scripts/setup_drugclip_cpu_wsl.sh
LFCOPY=/tmp/setup_drugclip_cpu_wsl.lf.sh

echo "=== frozen script hashes ==="
sha256sum "$FROZEN"

if [ ! -x /root/miniforge3/bin/micromamba ]; then
  echo "=== installing micromamba via conda ==="
  /root/miniforge3/bin/conda install -y -n base -c conda-forge micromamba || {
    echo "MICROMAMBA_INSTALL_FAILED"; exit 4; }
fi
export MAMBA_EXE=/root/miniforge3/bin/micromamba
"$MAMBA_EXE" --version || { echo "MICROMAMBA_NOT_RUNNABLE"; exit 4; }

echo "=== git safe.directory (repos owned by the Windows user, WSL runs as root) ==="
echo "HOME=$HOME  uid=$(id -u)"
printf '[safe]\n\tdirectory = *\n' > "${HOME:-/root}/.gitconfig"
cat "${HOME:-/root}/.gitconfig"
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0=safe.directory
export GIT_CONFIG_VALUE_0='*'
git config --get safe.directory || true
git -C "$SRC/project/tools/DrugCLIP" rev-parse HEAD
git -C "$SRC/project/tools/Uni-Core" rev-parse HEAD

echo "=== solver configuration (frozen YAML is realisable only with flexible channel priority) ==="
printf 'channel_priority: flexible\n' > "${MAMBA_ROOT_PREFIX}/.mambarc"
cat "${MAMBA_ROOT_PREFIX}/.mambarc"

echo "=== CRLF-stripped copy ==="
tr -d '\r' < "$FROZEN" > "$LFCOPY"
sha256sum "$LFCOPY"

echo "=== running the frozen setup script (LF copy) ==="
bash "$LFCOPY" "$SRC"
rc=$?
echo "SETUP_EXIT=$rc"
if [ "$rc" -ne 0 ]; then echo "STAGE3A_ENV_BUILD_FAILED"; exit "$rc"; fi

echo "=== environment identity ==="
/root/micromamba/envs/drugclip-cpu/bin/python - <<'PY'
import json, sys, rdkit, torch
info = {"python": sys.version.split()[0], "rdkit": rdkit.__version__,
        "torch": torch.__version__, "cuda": torch.cuda.is_available()}
try:
    import unicore; info["unicore"] = "ok"
except Exception as e:
    info["unicore"] = f"FAIL {type(e).__name__}: {e}"
print(json.dumps(info))
PY
echo "STAGE3A_ENV_READY"
