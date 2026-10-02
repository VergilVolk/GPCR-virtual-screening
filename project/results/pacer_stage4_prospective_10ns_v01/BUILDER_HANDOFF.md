# Stage4 prospective builder handoff

Implementation complete; full membrane build and CUDA smoke have not run locally.
No equilibration, release, production runner or server connection is invoked by
either new script. Production remains NOT_STARTED and launch readiness is false.

## Commands on the server

```bash
cd /data/GPCR-virtual-screening
conda activate pacer-dc-md
python project/scripts/build_pacer_stage4_prospective_membrane.py --all --dry-run
python project/scripts/build_pacer_stage4_prospective_membrane.py --candidate PACER0073
python project/scripts/run_pacer_stage4_membrane_smoke.py --candidate PACER0073 --device-index 0
python project/scripts/build_pacer_stage4_prospective_membrane.py --all
python project/scripts/run_pacer_stage4_membrane_smoke.py --all --device-index 0
```

Builder dry-run validates all three frozen candidates even when a single
candidate is selected. Missing OPM/probe assets are reported, never substituted.
`--prepare-only` writes rigidly oriented receptor and ligand inputs, with no
membrane build or dynamics. Both real build and prepare-only require the formal
OPM and ACh assets. Complete cluster systems can be reused after validating all
source/output hashes; partial cluster builds fail instead of mixing bases.

Required server inputs (in addition to frozen Stage3D tables/poses/receptors):

- `project/data/pdb/7TRS_OPM.pdb`, receptor chain R.
- `project/data/pdb/7TRS.pdb` and `project/data/pdb/m4_ligands/ACH_ideal.sdf`.
- `project/results/pacer_dc_openmm_reference_qc_v01/probe_only/minimized.pdb`:
  receptor chain E, ACh chain F. Explicit CLI overrides for reference receptor
  path/chain are available; no automatic fallback to a different reference.
- Activated `pacer-dc-md`, with the existing AmberTools binaries on PATH for
  AM1-BCC. No packages are installed by these scripts.

## Input and topology policy

Candidate coordinates are recovered directly from frozen PDBQT MODEL 1 using
the complete SMILES IDX mapping. Heavy atom count/elements, full state graph,
parent graph, formal charge and hashes must agree. PACER0027 retains its frozen
+2 microstate; it is not parameterized as the neutral parent. Undefined stereo
in SMILES is not promoted to chemical identity from 3D coordinates.

The source receptor includes native ACH despite its receptor filename. Remove
it from the protein/shared base, then add the mapped historical ACh only to
probe contexts. Two existing ACE/NME-capped fragments are retained. The source
TER records precede NME caps; the preparation copy moves separators after NME
so each cap is bonded to its own fragment. It does not bridge the missing loop,
delete the four additional receptor residues, or change protein coordinates.
Three-cluster ff14SB template tests validate this repair without dynamics.

Frozen receptor -> OPM uses 270 normalized-sequence CA correspondences, computes
Kabsch rotation/translation and RMSD, and preserves all 274 receptor residues.
Candidate uses this identical transform. Historical QC ACh is extracted using
the old 7TRS bond-assignment/OpenFF helpers, its QC receptor is aligned to the
corresponding frozen cluster CA atoms, then the same cluster -> OPM transform is
applied. Probe-only and candidate-probe receive the same mapped ACh template.

## Inherited scientific protocol

Existing forcefield, hydrogen-QC, positional restraint and solvent-clash helpers
are imported from `build_pacer_dc_membrane_reference.py`. Each cluster builds
one shared POPC/water/0.15 M ion base and forks four contexts. Protein ff14SB,
lipid17, TIP3P, OpenFF unconstrained 2.2.1 / AM1-BCC, PME/1 nm/5e-4, HBonds,
300 K, 1 bar XYIsotropic/ZFree membrane barostat and initial restraint k=1000
remain inherited. Only CPU restrained minimization (100 iterations) runs during
build. Parameter charges, input/minimized PDB, System/State XML, log and hashed
build receipts are written to the new namespace.

CUDA smoke stays bounded to 500 x 0.5 fs = 0.25 ps, with the historical
50/100/200/300 K schedule and disabled barostat. It rejects nonfinite energy,
forces, coordinates and constraint errors before stepping. Historical deep
minimization threshold/tolerance/iteration count are retained. It writes DCD,
CSV, refined/final states and a checkpoint, checking checkpoint roundtrip
without extending the 500-step length. These are smoke receipts, not equilibrated
production starting states. Replica equilibration/release and production
launch adaptation are deliberately outside this task.

## Local validation and limits

Syntax/import/CLI/dry-run, frozen identity and graph checks pass. Contract tests
cover corruption rejection, proper rigid transforms, capped fragment topology
and native probe exclusion for all clusters, and nonfinite smoke rejection.
Local full membrane build is blocked by absent OPM PDB, ACh template and QC
probe PDB; local CUDA smoke is NOT_RUN. No claim of a built 12-system campaign
or CUDA validation is made. Server synchronization is required.
