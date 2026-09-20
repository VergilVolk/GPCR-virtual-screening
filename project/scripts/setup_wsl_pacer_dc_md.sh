#!/usr/bin/env bash
#
# Bootstrap the PACER-DC molecular-dynamics environment on a Linux host (WSL2).
# =============================================================================
#
# Why this script exists
# ----------------------
# `project/environment_pacer_dc_md.yml` is necessary but not sufficient. Bringing
# the PACER-DC MD pipeline up from scratch also needs:
#
#   1. a *Linux* conda — AmberTools has no win-64 build, and `openmmforcefields`
#      hard-depends on it, so the frozen environment file cannot solve on Windows;
#   2. the repository on a native Linux filesystem — MD I/O over /mnt/c is far too
#      slow for trajectories;
#   3. `importlib_resources`, which openmmforcefields 0.14.1 imports unconditionally
#      but does not declare as a dependency;
#   4. six experimental structures (five from RCSB, one OPM-oriented) plus three
#      ligand template SDFs from the RCSB Chemical Component Dictionary;
#   5. the OPM-oriented 7TRS, for which no static download URL is documented.
#
# Every one of those was found the hard way; this script encodes the result.
#
# Usage
# -----
# Inside WSL (or any Linux host) with the repository reachable:
#
#   bash /mnt/c/projects/GPCR-virtual-screening/project/scripts/setup_wsl_pacer_dc_md.sh
#
# Overridable via environment variables:
#   SRC_REPO    source checkout        (default /mnt/c/projects/GPCR-virtual-screening)
#   WORK_ROOT   native FS parent       (default /root/work)
#   CONDA_ROOT  conda installation     (default /root/miniforge3)
#   ENV_NAME    conda environment      (default pacer-dc-md)
#
# Idempotent: each stage is skipped when its output already exists, so rerunning
# after a partial failure resumes rather than restarting.
#
set -euo pipefail

SRC_REPO="${SRC_REPO:-/mnt/c/projects/GPCR-virtual-screening}"
WORK_ROOT="${WORK_ROOT:-/root/work}"
CONDA_ROOT="${CONDA_ROOT:-/root/miniforge3}"
ENV_NAME="${ENV_NAME:-pacer-dc-md}"
REPO="$WORK_ROOT/GPCR-virtual-screening"
ENV_PREFIX="$CONDA_ROOT/envs/$ENV_NAME"

say()  { printf '\n=== %s ===\n' "$*"; }
ok()   { printf '  [ok]   %s\n' "$*"; }
skip() { printf '  [skip] %s\n' "$*"; }
die()  { printf '  [FAIL] %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 0. Preconditions
# ---------------------------------------------------------------------------
say "0. preconditions"
[[ "$(uname -s)" == "Linux" ]] || die "this script targets Linux; run it inside WSL2"
[[ -d "$SRC_REPO" ]] || die "source repository not found at $SRC_REPO (set SRC_REPO=...)"
ok "running on Linux, source repository visible"

# ---------------------------------------------------------------------------
# 1. conda (miniforge)
# ---------------------------------------------------------------------------
say "1. miniforge"
if [[ -x "$CONDA_ROOT/bin/conda" ]]; then
  skip "conda already present: $("$CONDA_ROOT/bin/conda" --version)"
else
  curl -fsSL --retry 3 --connect-timeout 20 -o /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash /tmp/miniforge.sh -b -p "$CONDA_ROOT"
  "$CONDA_ROOT/bin/conda" config --set channel_priority strict
  ok "installed $("$CONDA_ROOT/bin/conda" --version) at $CONDA_ROOT"
fi

# ---------------------------------------------------------------------------
# 2. Repository on a native Linux filesystem
# ---------------------------------------------------------------------------
say "2. repository on native filesystem"
fs_type="$(df -T "$WORK_ROOT" 2>/dev/null | tail -1 | awk '{print $2}' || echo unknown)"
if [[ -d "$REPO/.git" || -d "$REPO/project" ]]; then
  skip "repository already at $REPO"
else
  mkdir -p "$WORK_ROOT"
  cp -a "$SRC_REPO" "$REPO"
  ok "copied to $REPO ($(du -sh "$REPO" | cut -f1))"
fi
case "$(df -T "$REPO" | tail -1 | awk '{print $2}')" in
  9p|drvfs) die "repository sits on $fs_type (Windows mount); copy it into the WSL filesystem" ;;
  *) ok "filesystem is $(df -T "$REPO" | tail -1 | awk '{print $2}')" ;;
esac

# ---------------------------------------------------------------------------
# 3. conda environment
# ---------------------------------------------------------------------------
say "3. conda environment '$ENV_NAME'"
if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  skip "environment already exists at $ENV_PREFIX"
else
  "$CONDA_ROOT/bin/conda" env create -f "$REPO/project/environment_pacer_dc_md.yml" -y
  ok "environment created"
fi
# Repair the undeclared runtime dependency (see header note 3).
if ! "$ENV_PREFIX/bin/python" -c "import importlib_resources" 2>/dev/null; then
  "$CONDA_ROOT/bin/conda" install -n "$ENV_NAME" -c conda-forge importlib_resources -y
  ok "installed importlib_resources"
else
  skip "importlib_resources already present"
fi

# ---------------------------------------------------------------------------
# 4. Verification
# ---------------------------------------------------------------------------
say "4. verification"
# Delegate to the repository's own audit rather than duplicating its checks.
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

for b in sqm antechamber tleap parmchk2; do
  command -v "$b" >/dev/null || die "$b missing from PATH — is the environment activated?"
done
ok "AmberTools binaries on PATH (required for AM1-BCC charges)"

"$ENV_PREFIX/bin/python" - <<'PYEOF' || die "import check failed"
import importlib, openmm
for m in ["openmm", "openmm.app", "openff.toolkit", "openff.units",
          "openmmforcefields.generators", "pdbfixer", "parmed",
          "MDAnalysis", "rdkit", "numpy", "pandas", "scipy"]:
    importlib.import_module(m)
print("  [ok]   all pipeline imports resolve; openmm", openmm.version.version)
print("  [ok]   platforms:", [openmm.Platform.getPlatform(i).getName()
                              for i in range(openmm.Platform.getNumPlatforms())])
PYEOF

# ---------------------------------------------------------------------------
# 5. Experimental structures
# ---------------------------------------------------------------------------
say "5. structures from RCSB"
cd "$REPO"
mkdir -p project/data/pdb/m4_ligands
for id in 7TRS 7TRQ 7TRP 7V68 7V6A; do
  dest="project/data/pdb/${id}.pdb"
  if [[ -s "$dest" ]]; then skip "$id.pdb"; continue; fi
  curl -fsSL --retry 3 -o "$dest" "https://files.rcsb.org/download/${id}.pdb"
  ok "$id.pdb ($(wc -c < "$dest") bytes)"
done

say "5b. ligand templates from the RCSB Chemical Component Dictionary"
for ccd in ACH 2CU 5XI; do
  dest="project/data/pdb/m4_ligands/${ccd}_ideal.sdf"
  if [[ -s "$dest" ]]; then skip "${ccd}_ideal.sdf"; continue; fi
  curl -fsSL --retry 3 -o "$dest" "https://files.rcsb.org/ligands/download/${ccd}_ideal.sdf"
  ok "${ccd}_ideal.sdf ($(wc -c < "$dest") bytes)"
done

# ---------------------------------------------------------------------------
# 6. OPM-oriented 7TRS
# ---------------------------------------------------------------------------
# The membrane builder aligns the receptor to the OPM frame, so it needs the
# OPM-oriented coordinates. OPM's web front end is JS-rendered and has no
# documented static URL, so use moleculekit's own OPM fetcher.
# Installed into the *base* environment on purpose: it is a one-off fetch tool
# and must not perturb the frozen MD environment.
say "6. OPM-oriented 7TRS"
OPM_DEST="project/data/pdb/7TRS_OPM.pdb"
if [[ -s "$OPM_DEST" ]]; then
  skip "7TRS_OPM.pdb already present"
else
  "$CONDA_ROOT/bin/python" -m pip install --quiet moleculekit
  "$CONDA_ROOT/bin/python" - "$OPM_DEST" <<'PYEOF'
import sys
from moleculekit.opm import get_opm_pdb
mol, thickness = get_opm_pdb("7TRS")
mol.write(sys.argv[1])
print(f"  [ok]   OPM 7TRS: {mol.numAtoms} atoms, thickness {thickness:.2f} A")
PYEOF
  ok "wrote $OPM_DEST"
fi
ca_r="$(awk '/^ATOM/ && substr($0,22,1)=="R" && substr($0,13,4) ~ /CA/' "$OPM_DEST" | wc -l)"
[[ "$ca_r" -gt 0 ]] || die "no chain-R CA atoms in $OPM_DEST; the membrane alignment would fail"
ok "chain R CA atoms: $ca_r"

# ---------------------------------------------------------------------------
# 7. Readiness audits
# ---------------------------------------------------------------------------
say "7. readiness audits"
"$ENV_PREFIX/bin/python" project/scripts/audit_pacer_dc_ligand_parameterization.py >/dev/null
"$ENV_PREFIX/bin/python" project/scripts/audit_pacer_dc_md_environment.py >/dev/null
"$ENV_PREFIX/bin/python" - <<'PYEOF'
import json, pathlib
p = json.loads(pathlib.Path("project/results/pacer_dc_ligand_parameterization_audit.json").read_text())
e = json.loads(pathlib.Path("project/results/pacer_dc_md_environment_audit.json").read_text())
print(f"  ligand parameterization : {p['success_count']}/{p['ligand_count']} passed, all_passed={p['all_passed']}")
print(f"  production_ready        : {e['production_ready']}  ({e['blocking_reason']})")
if not p["all_passed"]:
    raise SystemExit("ligand parameterization smoke test did not pass")
PYEOF

say "done"
cat <<'EOF'
  The MD environment is ready. Remaining phase-1 blockers are scientific, not
  environmental:
    - no candidate has complete matched-context evidence
    - fewer than three independent PAM positive-control chemotypes
    - no confirmed binding-but-non-PAM control
  See project/docs/PACER_DC_GPU_HANDOFF.md before launching production runs.
EOF
