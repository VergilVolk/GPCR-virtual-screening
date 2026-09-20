# PACER-DC GPU Phase-1 handoff

## Purpose

Run the frozen 18-job pilot: six matched contexts × three independent replicas × 5 ns. This stage tests trajectory integrity and whether a reproducible dynamic signal exists. It is not PAM validation.

## Required environment

- Linux;
- one CUDA-capable GPU per task;
- OpenMM with the `CUDA` platform;
- the packages frozen in `environment_pacer_dc_md.yml`.

### The environment MUST be activated, not just installed

AmberTools supplies the `sqm` / `antechamber` binaries that OpenFF uses for
AM1-BCC partial charges. Those binaries live in `$CONDA_PREFIX/bin`, so calling
the interpreter by absolute path without activating makes OpenFF fall back to
RDKit-only charge methods and every ligand fails with:

```
ValueError: No registered toolkits can provide the capability "assign_partial_charges"
  partial charge method 'am1bcc' is not available from RDKitToolkitWrapper
```

This looks like a broken environment but is only a PATH problem. Always run:

```bash
source <prefix>/etc/profile.d/conda.sh
conda activate pacer-dc-md
```

Confirm before doing anything else:

```bash
for b in sqm antechamber tleap parmchk2; do printf '%-10s %s\n' "$b" "$(command -v $b || echo MISSING)"; done
```

`scripts/audit_pacer_dc_md_environment.py` reports `tleap: null` and
`production_ready: false` when the environment was not activated. That is a
tooling artefact, not a real readiness failure.

### Preflight

```bash
python - <<'PY'
from openmm import Platform
print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
assert "CUDA" in [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
PY
```

Note that listing `CUDA` only proves the platform was compiled in. It does not
prove the platform can run. Builds whose CUDA toolkit is newer than the
installed driver fail later with
`CUDA_ERROR_UNSUPPORTED_PTX_VERSION`, so also create a real context:

```bash
python - <<'PY'
from openmm import Platform, System, LangevinMiddleIntegrator, Context, unit
p = Platform.getPlatformByName("CUDA")
s = System(); s.addParticle(1.0)
i = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
Context(s, i, p)
print("CUDA context OK")
PY
```

## Slurm launch

From the repository root:

```bash
mkdir -p project/results/pacer_dc_production_v01
export PACER_PYTHON=/absolute/path/to/pacer-dc-md/bin/python
sbatch project/scripts/run_pacer_dc_gpu_phase1.sbatch
```

The array uses `project/config/pacer_dc_gpu_phase1_manifest.csv`. Do not edit seeds, contexts, or target lengths after launch.

## Resume behavior

Rerunning the same array is safe. Each task first loads `checkpoint.chk`; if a different GPU/node cannot read it, the runner falls back to the portable `latest_state.xml` plus the progress ledger. Completed jobs exit without recomputation.

## Frozen Phase-1 gate

After all jobs finish:

```bash
python project/scripts/audit_pacer_dc_production_phase1.py
```

All 18 jobs must pass:

- completed target length;
- finite trajectory log;
- median temperature 285–315 K;
- at least 95% of reports within 270–330 K;
- median density 0.90–1.10 g/mL;
- absolute box-volume drift no greater than 5%;
- expected frame count.

If any task fails, do not compute paired biological differences until the cause is resolved.

## Endpoint extraction

For the 5 ns pilot:

```bash
python project/scripts/extract_pacer_dc_production_endpoints.py \
  --evidence-level production_pilot
```

The resulting records cannot satisfy `pacer_dc_score.py`, by design. Formal `trajectory` evidence requires the longer registered sampling, convergence audit, and explicit override.

## Claim boundary

Passing Phase-1 means the systems are numerically suitable for extension. It does not demonstrate conformational convergence, distinguish PAM from ago-PAM, establish predictive performance, or confirm a candidate molecule.
