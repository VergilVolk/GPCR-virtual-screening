# FINAL_CODE_FREEZE_CANDIDATE_REPORT

**Phase:** FINAL_PIPELINE_INTEGRATION_PREFLIGHT - minimal engineering integration
**Date:** 2026-10-02
**Status:** **READY_TO_TAG_CODE_FREEZE**
**Tag created:** NO (as instructed) - no tag, no push, no prospective generation

---

## 0. Scope and prohibitions honoured

Authorised: write **only** inside the new integration worktree.
Not executed: molecule generation, candidate DrugCLIP batch scoring, docking campaign,
MD production, new-candidate PACER-FKG inference, final freeze tag, push.

Integration worktree: `C:\projects\GPCR-virtual-screening-final`
Branch: `integration/final-freeze-v01`
Base commit: `18100f9b6b42f27e43760daa5b3dfdd028b330d3`

The other three worktrees were **not touched**:

| Worktree | HEAD before | HEAD after | dirty |
|---|---|---|---|
| `GPCR-virtual-screening` | `5e86149b` | `5e86149b` | 0 |
| `GPCR-virtual-screening-fkg-v02` | `853aaa12` | `853aaa12` | 0 |
| `GPCR-virtual-screening-encoder` | `0a4b4f23` | `0a4b4f23` | 0 |

---

## A. WSL historical assets protected and inventoried

Per correction 1, `/root/work/GPCR-virtual-screening` is treated as a **historical
execution workspace**. Nothing in it was reverted, overwritten, moved or deleted.
A read-only manifest / backup receipt was produced instead.

| Item | Value |
|---|---|
| Branch / HEAD | `main` / `5ef9904cc3ca688fac7e4521473d54f3cb09c345` (2026-09-20 17:51) |
| Relationship to the Windows repo | **same repository** - `5ef9904` is an ancestor of every current branch |
| Working-tree size | 24,266,102,850 bytes (~22.6 GiB) |
| `project/results` | **514 files, 24,250,250,806 bytes** |
| Git state | 333 porcelain entries (327 ` M`, 6 `??`), empty stash |
| Real content differences | **5 files** (the other 322 differ only by line endings / file mode) |

Receipts written to `project/integration_freeze/receipts/`:
`WSL_GIT_HEAD_v01.txt`, `WSL_GIT_STATUS_PORCELAIN_v01.txt` (333 lines),
`WSL_GIT_DIFF_NUMSTAT_v01.txt`, `WSL_GIT_DIFF_STAT_IGNORE_CR_v01.txt`,
`WSL_DIRTY_CODE_DIFF_v01.patch` (full patch for the 5 files),
`WSL_DIRTY_FILES_BLOB_HASHES_v01.tsv`, `WSL_RESULTS_SHA256_v01.txt` (**all 514 files**),
`WSL_RESULTS_FILES_v01.tsv`, `WSL_TOPFILES_v01.tsv`, `WSL_WORK_EXTRAS_SHA256_v01.txt`,
`WSL_GIT_{BRANCH,BRANCHES,REMOTES,LOG1,STASH}_v01.txt`.
Human-readable summary: `WSL_HISTORICAL_WORKSPACE_RECEIPT_v01.md`.

**Unique assets found only in WSL** (not in any Windows worktree):
`pacer_dc_production_v01` (70 files, 13.34 GB, 5 ns pilot set, 6 systems x 3 replicas),
`pacer_dc_restraint_release_v01` (84 files, 1.79 GB, **the only complete 6-system set**),
`pacer_dc_membrane_reference_v01` (22 files, 0.99 GB, **the only complete 6-system set**),
plus `pacer_dc_recovery*`, `pacer_dc_four_context_v01`, `pacer_dc_short_equilibration_v01`,
`pacer_dc_membrane_smoke_v01`, `pacer_dc_unified_qc_v01`, `pacer_dc_common_protein_v01`,
`pacer_dc_md_environment_audit.json`.

**Declared hazard preserved as found:** the WSL copies of
`build_pacer_dc_membrane_reference.py` and `build_pacer_dc_openmm_reference_systems.py`
carry abandoned in-progress edits (the ligand `definitions` block removed) that exist
in **no commit**. They are recorded, not repaired, and are **not** the code that produced
the historical membrane references.

**Forward rule:** production MD will use a clean copy of the code-freeze tree
re-synchronised onto a native Linux filesystem, never this tree.

---

## B / C. Integration worktree and branch

| Step | Result |
|---|---|
| `git worktree add --detach C:\projects\GPCR-virtual-screening-final 18100f9` | exit 0 |
| `git switch -c integration/final-freeze-v01` | exit 0 |
| `git worktree list` | 4 worktrees, FINAL at `18100f9 [integration/final-freeze-v01]` |

`18100f9` was chosen because the preflight proved it is a **strict descendant of every
relevant head** (`5e86149`, `853aaa1`, `0a4b4f2`, `c19b744`, `a2e01fe`, `07ce8ae`, `6d12bce`),
so **no merge was needed and none was performed**. Delta vs the frozen Stage-B head:
374 added, 3 modified, 0 deleted.

### Asset links (NTFS junctions, declared)

The worktree carries the integrated **tracked** tree only. Ignored assets are reached
through junctions; none of the targets is modified:

| Junction in FINAL | Target |
|---|---|
| `project/cache` | `-fkg-v02\project\cache` |
| `project/results/pacer_fkg_v02_longmd_v01` | `-fkg-v02` counterpart |
| `project/results/pacer_dc_membrane_reference_v01` | `-fkg-v02` counterpart |
| `project/results/pacer_dc_cm00734_stage_b_20ns_v01/production` | `GPCR-virtual-screening` counterpart |
| `project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/track_b_fkg_v02` | `-fkg-v02` counterpart |
| `project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/topology_sanitized` | `-fkg-v02` counterpart |
| `project/data/pdb` | `GPCR-virtual-screening\project\data\pdb` |

---

## D / E. Manifest-driven `frozen_runner`

New package `project/pacer_fkg_v02/frozen_runner/` - **no frozen numeric was edited**.

| File | Role |
|---|---|
| `engines.py` | import-only re-export layer; `assert_frozen_numerics()`, `frozen_numerics_report()`, `frozen_module_ledger()`, `cached_anchor()` |
| `spec.py` | fail-closed run-spec loader/validator + `authenticate_inputs()` |
| `paths.py` | output-path policy; protected-root list **imported** from the frozen guard |
| `blocks.py` | parameterised block constructor + `prove_equivalence()` |
| `io_utils.py` | safe `.npy` loading (`load_npy_copy`, `check_npy`) |
| `fkg_phase1.py` | generic C1-BS256 extraction driver |
| `fkg_phase2a.py` | generic frozen normalization + RFF apply |
| `fkg_phase2b.py` | generic graph diffusion + region pooling + contrasts |
| `cli.py` | `python -m project.pacer_fkg_v02.frozen_runner.cli --spec ... --phase {1,2a,2b,all} --run|--verify` |
| `regression.py` | CM00734 payload-identity harness |

### Frozen values are imported, never redefined

```
REQUIRED_FREEZE_SHA256 b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd
REPLICAS (1,2,3) | ALPHA 0.65 | DIFFUSION_STEPS 20 | N_RESIDUES 270 | RFF_FEATURES 512
BRANCHES STATE_MOTION {width 512, bandwidth 32.30188361260893, seed 272340}
         SIGNED_DRIFT {width 256, bandwidth 29.29934899925964, seed 272084}
CONTRASTS Delta_PAM {CP:+1,P:-1} | Delta_AGO {C:+1,A:-1} | Delta_INT {CP:+1,P:-1,C:-1,A:+1}
P2: N_FRAMES 1000 | BLOCK_FRAMES 20 | N_BLOCKS 50 | BASE_SEED 271828
```

### Fail-closed validators

Run specs must declare `replicas == [1,2,3]` and the frozen paired seed group
`{1:27101, 2:38201, 3:49301}` (required for context matching with the historical
shared controls); `block_frames == 20`; `n_frames \u2208 {400,1000}`;
`blocks_per_trajectory * block_frames == n_frames`; `calibration.usage \u2208 {forbidden, verify_only}`
with `r2_fitting_allowed == false`; nine forbidden scope claims must all be `false`;
every declared input hash must be a 64-hex digest and must match on disk.
Output paths may never fall inside a protected historical root and may never overwrite
an existing artifact.

### Correction 2 (R2 semantics) implemented exactly

R2 is **never used for fitting**. It is still mapped through the frozen
normalization/RFF state and reported as a descriptive replica result, preserving
Stage-A/Stage-B semantics. The primary cross-replica metric remains the **R1/R3
direction cosine**; `calibration.usage` for the regression spec is `forbidden`.

---

## F. DrugCLIP apply-only production path

New package `project/drugclip_freeze/`.

### Model A - official frozen apply-only

| Must-restore asset | Identity | State |
|---|---|---|
| DrugCLIP source | `bowen-gao/DrugClip @ 7a3a3fa33673f8668c811790f2e4681c98af44ef` | absent |
| Uni-Core source | `dptech-corp/Uni-Core @ 44f6386f4dcd7137fc1e5d5e768117d635d64a26` | absent |
| Official checkpoint | `checkpoint_best.pt`, 1,183,713,459 B, sha256 `dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e` | absent |
| `drugclip-cpu` WSL env + micromamba | per `setup_drugclip_cpu_wsl.sh` | absent |
| Molecule / pocket LMDBs | **rebuilt deterministically** (1 conformer, seed 20260922) | n/a |
| Source PDBs 7TRQ / 7TRP / 7TRS | hashes pinned | **present, hash-verified** |
| Model B projection weights | 9 files, 2.63 MB each | **present (tracked)** |

`model_a.validate_archive()` enforces: required archive keys, unique molecule ids,
`(n,128)`/`(n,512)` shapes, `scores == (n_pockets, n_molecules)`, finiteness, the
non-empty M4 pocket set, and **empty `missing_checkpoint_keys` /
`unexpected_checkpoint_keys`**. The frozen CPU patch in
`prepare_drugclip_cpu_checkout.py` applies exactly four device-placement
substitutions, each required to occur exactly once.

### Model B - frozen projection-only

`projection.py` re-expresses the deployed forward pass functionally
(`linear1 512x512 -> ReLU -> linear2 128x512 -> L2-normalise`), validating every
state dict key-by-key and shape-by-shape. `model_b.score()` consumes Model A's
frozen 512-D `molecule_representations` / `pocket_representations`, takes the
**official max over the declared M4 pockets**, and averages the three deployed seeds.

- the backbone is never re-run;
- no gradient is computed;
- **no training module is imported** into the production path.

### Dual-model decision

`old_top50 INTERSECT new_top25`, Model A = M4 decision line, Model B = second opinion.
**Score fusion is forbidden** and `assert_no_fusion` guards the decision inputs.
`verify_committed_shortlist()` re-derives the rule against the committed
`candidate_dual_model_shortlist_v02.csv`; the test suite proves the rule reproduces
(`_merge == "both"` exactly when `drugclip_rank <= 50 and new_rank <= 25`).

---

## G. Run-spec / manifest schema

- `project/integration_freeze/schemas/final_pipeline_run_spec_v1.schema.json` -
  JSON Schema draft-07 with `const` constraints on every frozen field
  (seeds, replicas, block geometry, freeze SHA256, encoder checkpoint SHA256).
- `project/integration_freeze/run_specs/CM00734_REGRESSION_v01.json` - the
  verification/replay spec; `scope.regression_replay_only = true`,
  `scope.new_md_production_run = false`, `calibration.usage = "forbidden"`,
  and 7 declared inputs each with an independently re-verified SHA256.

---

## H. CM00734 regression verification (correction 3)

**Verification / replay only. No new MD was run.**

### H.1 Frozen anchor

`engines.verify_frozen_anchor()` passes **strictly and unmodified**:
`{"status": "VERIFIED_FROZEN", "phase1_pairs": 12, "trajectory_blocks": 600,
"calibration_blocks": 200, "fitting_replicas": [2], "excluded_replicas": [1, 3],
"v01_protected_files": 41}`.

### H.2 Block-constructor equivalence

`prove_equivalence()` compared the parameterised constructor against **both** frozen
implementations on random data:

| Check | Result |
|---|---|
| `phase2_1000x50_STATE_MOTION_bit_identical` | true |
| `phase2_1000x50_SIGNED_DRIFT_bit_identical` | true |
| `stage_b_400x20_STATE_MOTION_bit_identical` | true |
| `stage_b_400x20_SIGNED_DRIFT_bit_identical` | true |
| endpoint identity (1000/50 and 400/20) | true |

### H.3 Payload identity

Frozen `STAGE_B_MATCHED_20NS_RESULTS_v01.json` SHA256 =
`be663c7d45cc4ade1e053029e6aa09e65ca24c1e2fc22dae5398a26bf49d7d99` - **matches the
recorded anchor**.

| Metric | Value |
|---|---|
| **Exact fields compared** | **1,632** |
| **Exact fields mismatched** | **0** |
| Cosine fields compared | 108 |
| Cosine max abs diff | **6.44e-08** |
| Cosine sign mismatches | **0** |
| Cosine 6-dp rendering differences | 2 (1 ULP, last digit) |
| **`scientific_payload_identical`** | **true** |

**Primary Stage-B signal** (`STATE_MOTION / CM00734 / Delta_INT / compound110_extension`):

| | value |
|---|---|
| frozen | `-0.3268392457993834` |
| recomputed | `-0.32683923840522766` |
| abs diff | `7.39e-09` |
| sign preserved | **true** |
| rendered at 6 dp | `-0.326839` (both) |

### H.4 Declared comparison policy

Per correction 3, the acceptance criterion is **scientific payload identity**, not byte
identity of receipt files. The policy is explicit in the receipt:
all descriptor / block-summary quantities must match **exactly**; the single tolerated
quantity is the scalar `R1_R3_direction_cosine` (one float32 inner product, hence
summation-order noise of order 1e-8) with an absolute tolerance of **1e-6**, i.e. four
orders of magnitude below the 6-decimal precision at which the frozen reports publish
these cosines. Sign preservation is a hard gate; 6-dp rendering equality is
**reported as evidence, not gated**, because a 1-ULP difference can flip the last digit.

Non-scientific provenance fields deliberately excluded from the comparison are listed
in the receipt: `created_at`, `upstream.phase2a_receipt_sha256`,
`vector_artifacts[].{path,bytes,sha256}`, `graph.path`, run-level `matched_design`
labels, `version`, `schema`, `status`, `run_id`.

### H.5 Frozen encoder replay (Phase 1)

`verify_phase1_replay.py` re-ran the **frozen C1-BS256 extractor on CUDA** over the
first four frames of two existing CM00734 trajectories and compared with the frozen
Phase-1 cache:

| System | Replica | relative RMSE | limit | passed |
|---|---|---:|---:|---|
| `CM00734__candidate_no_probe` | 1 | 1.770e-07 | 1e-4 | **true** |
| `CM00734__candidate_no_probe` | 2 | 1.764e-07 | 1e-4 | **true** |

Checkpoint SHA256 verified, `strict_load_ok = true`,
`geom2vec_source_equivalent = true`, `bitwise_identical = false` with max abs
difference 2.1e-06 - i.e. GPU float non-determinism far inside the frozen
relative-RMSE acceptance criterion.

---

## I. Tests and freeze-consistency checks

| Suite | Result |
|---|---|
| `project/tests/test_frozen_runner.py` | **16 / 16 passed** |
| `project/tests/test_drugclip_freeze.py` | **9 / 9 passed** |

Highlights: `test_frozen_anchor_verifies` (strict anchor),
`test_frozen_numerics_match_freeze_manifest`,
`test_frozen_module_ledger_matches_recorded`, both bit-identity block tests,
six spec-rejection tests, three path-policy tests,
`test_pinned_pdbs_match`, `test_committed_shortlist_reproduces_rule`.

Freeze-consistency: `git diff` for the four frozen v02 modules is **empty**;
their LF-normalised SHA256s are recorded in `FROZEN_MODULE_LEDGER_v01.json`:

```
run_phase1_bs256.py       d3ed9764ff813b9f164e139668cab2f55d7112d7b04283812d4a5fe3beee2271
run_phase2_calibration.py ccf8af5247ff2c7325f1c451d36943f2e293dee0dff48625031d8bb10a2d6e11
run_phase3_fkg.py         147cb31363fe733c72719494dd282c80c86a61e1b733ad2568b94eca3b1d224c
path_guard.py             144e3bf7bb45523246c45844c8a3a33b20819b4289c169314c12655bc8341a33
```

---

## NEW DEFECT FOUND AND REPAIRED (must be carried into the tag message)

### The frozen v02 anchor is line-ending sensitive

`V02_FREEZE_MANIFEST.json` records a raw SHA256 **and** a byte count for every artifact
it protects. Those values were produced from the `-fkg-v02` working tree. A fresh
Windows checkout under `core.autocrlf=true` materialises CRLF, so the raw hash and the
byte count differ and `run_phase2_calibration.verify()` fails closed:

```
RuntimeError: frozen artifact mismatch:
  .../project/pacer_fkg_v02/run_phase2_calibration.py
```

**The freeze pins a MIXED line-ending state:**

| File | Pinned as | Pinned bytes | LF bytes |
|---|---|---:|---:|
| `run_phase2_calibration.py` | **LF** | 21,545 | 21,545 |
| `analyze_pacer_fkg_g2b.py` | **CRLF** | 15,688 | 15,387 |
| `analyze_pacer_fkg_g2c.py` | **CRLF** | 17,260 | 16,955 |

Consequently **no blanket `.gitattributes` rule can fix this** - forcing `eol=lf`
would repair the first file and break the other two. A human decision is required.

**Repair applied in this worktree** (an environment normalisation, not a content change):
`project/integration_freeze/tools/normalize_frozen_eol.py` rewrites a file to LF **only
when** `lf_sha256(file) == recorded sha256 AND lf_bytes(file) == recorded bytes`, so the
operation is self-proving. It examined 68 protected records: **67 already satisfied,
1 normalised** (`run_phase2_calibration.py`, 21,993 -> 21,545 bytes), **0 unsatisfied**,
**0 missing**. Evidence: `FROZEN_EOL_NORMALISATION_{DRYRUN,APPLIED}_v01.json`.

Proof that no frozen content changed:
`git hash-object` (working tree) = `6baa42d6a2cdda8e17c95c43e0c97f282406ef72` =
`git rev-parse HEAD:<file>`; `git diff --stat` is empty.

---

## Obligations that the freeze tag must carry

1. **Line endings.** Before running the frozen chain from any fresh checkout, verify the
   three `source_code` artifacts against the freeze manifest and run
   `normalize_frozen_eol.py` if needed. The underlying mixed-EOL defect needs a human
   ruling (see F-1 in the risks).
2. **DrugCLIP restoration.** Model A cannot be executed until the four assets in
   section F are restored and the archive validator has been run once.
   `drugclip_freeze/RUNBOOK.md` is the frozen procedure.
3. **IFP** is declared **OPTIONAL STRUCTURAL EVIDENCE**; **multi-conformation GaMD
   gating** is declared **BLOCKED / HISTORICAL for v01** (no provenance pin exists for
   the third-party ensemble repository). No substitute method may be introduced.
4. **Generation** is a new prospective run. It must not claim byte-identical recovery of
   the historical generated library; its input provenance and every generation parameter
   must be frozen into the run spec before the first execution.
5. **Claim boundary.** PACER-FKG output for a new candidate is prospective
   representation-space mechanistic / prioritization evidence only. No magnitude
   threshold, no PAM probability, no outcome-driven rule.

---

## Deliverable inventory (all new, all untracked, nothing committed)

```
project/pacer_fkg_v02/frozen_runner/   __init__.py engines.py spec.py paths.py blocks.py
                                       io_utils.py fkg_phase1.py fkg_phase2a.py fkg_phase2b.py
                                       cli.py regression.py
project/drugclip_freeze/               __init__.py projection.py model_a.py model_b.py
                                       dual_model.py apply.py RUNBOOK.md
project/integration_freeze/            README.md FROZEN_MODULE_LEDGER_v01.json
                                       RECEIPT_INDEX_v01.json
                                       WSL_HISTORICAL_WORKSPACE_RECEIPT_v01.md
                                       schemas/final_pipeline_run_spec_v1.schema.json
                                       run_specs/CM00734_REGRESSION_v01.json
                                       tools/{normalize_frozen_eol,run_regression,
                                              verify_phase1_replay,diagnose_regression}.py
                                       tools/_shim_helper.py
                                       receipts/ (23 files)
project/tests/                         test_frozen_runner.py test_drugclip_freeze.py
```

`project/integration_freeze/RECEIPT_INDEX_v01.json` records every receipt with its
size and SHA256.
