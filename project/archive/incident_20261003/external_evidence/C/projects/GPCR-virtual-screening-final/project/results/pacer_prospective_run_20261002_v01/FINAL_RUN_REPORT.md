# FINAL RUN REPORT - PACER prospective run 2026-10-02 (Stages 1-3)

## Final status

| Stage | Status |
|---|---|
| STAGE 1 - Prospective molecule generation | **STAGE1_GENERATION_COMPLETE** |
| STAGE 2 - Physicochemical / drug-likeness filtering | **STAGE2_FILTERING_COMPLETE** |
| STAGE 3 - Frozen DrugCLIP dual-model screening | **STAGE3_DRUGCLIP_COMPLETE** |

Execution was strictly sequential and fail-closed. All three stages completed. No docking, IFP,
diversity selection, MD or PACER-FKG was executed.

Production baseline: branch `integration/final-freeze-v01`, commit
`ddc11f970f508d89a75dac8ca30c41cfcd361731`, tag `pacer-code-freeze-v01`.
The frozen scientific assets were not modified.

---

## Stage 1 - generation

| Field | Value |
|---|---|
| Implementation | `project/scripts/fragment_generation.py` (unchanged) |
| Seeds | 12, derived from the tracked frozen dataset `modeling_potency_exact_calcium.csv` (top-12 by max pEC50, deduplicated by molecule) |
| Fallback seeds used | **no** |
| Effective RNG pin | `PYTHONHASHSEED=0` |
| Raw generated | **5000** (generation cap reached; 17,146 / 200,000 join attempts) |
| Valid | **2605** |
| Invalid | **2395** (pharmacophore / drug-likeness rejections inside the generator) |
| Duplicates | **0** (raw and canonical) |
| Unique valid | **2605** |
| Output | `project/results/generated/generated_pam_analogs.csv` |
| Output SHA256 | `a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d` |
| QC | **PASS** - 12/12 seeds reconstructed, 0 invalid SMILES, 0 duplicates, 0 rule violations, determinism proven by two byte-identical runs |

**Defect found:** the frozen generator is non-deterministic under default hash randomisation
(two unpinned runs: 3212 vs 930 passing molecules). `PYTHONHASHSEED` is the effective random seed
because the fragment pool is built by iterating a `set` of strings. Pinned per the instruction to
freeze random seeds before execution; no scientific parameter changed.

---

## Stage 2 - filtering

| Field | Value |
|---|---|
| Implementation | `project/scripts/build_candidate_portfolio.py` (unchanged) |
| Input | Stage-1 output of this run only, SHA256 `a832c9cb...` |
| Input count | **2605** |
| Pass count | **200** (frozen quota: 160 local + 40 exploratory) |
| Rejection count | **2405** |
| Rejections by primary reason | exact known molecule 13; outside applicability domain 931; drug-likeness 0; strict-inactive risk 15; portfolio quota cap 1446 (sums exactly to 2405) |
| Rejections by rule (conjunctive) | exact known 13; PAINS 0; druglike composite 0; outside AD 931; risk >= 0.45: 55 |
| Eligible before quota | 1646 |
| Unique final molecules | **200** |
| Candidate IDs | `PACER0001` .. `PACER0200`, deterministic |
| Output | `project/results/pacer_candidates_v01/predock_portfolio.csv` |
| Output SHA256 | `0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b` |
| QC | **PASS** - 200 unique molecules, 0 duplicate IDs, 0 independent rule violations, exact rejection accounting |

**Environment provisioning:** `lightgbm 4.6.0` installed into conda env `chrm4_vs` with `--no-deps`
(absent from every environment; the frozen script imports it at module level). No existing package
changed. The affected model is reference-only and is not read by the eligibility filter.

---

## Stage 3 - frozen DrugCLIP

### 3A asset restoration / validation

| Asset | Required | Result |
|---|---|---|
| DrugCLIP source | `7a3a3fa33673f8668c811790f2e4681c98af44ef` | **VERIFIED** |
| Uni-Core source | `44f6386f4dcd7137fc1e5d5e768117d635d64a26` | **VERIFIED** |
| Official checkpoint | 1,183,713,459 bytes, sha256 `dc2c76d0...f667e` | **VERIFIED - exact byte and SHA256 match** |
| Audited CPU patch | 4 device-placement substitutions, each unique | **APPLIED** (architecture and tensors unchanged) |
| Environment `drugclip-cpu` | from `environment_drugclip_cpu.yml` unchanged | **FAILED - identity not realisable** |

The frozen environment spec cannot be realised on this host:

- with the default channel priority micromamba reports
  `pytorch =2.0.1 * is not installable because it conflicts with any installable versions previously reported`;
- with flexible channel priority the environment resolves and installs, but the pytorch 2.0.1 CPU
  binary cannot load: `libtorch_cpu.so: undefined symbol: iJIT_NotifyEvent` (the Intel OpenMP runtime
  it links against is no longer distributed through the channels used).

Substituting a different pytorch build would change environment identity and was **not** performed.

### 3B / 3C / 3D

**NOT STARTED.** The apply-only code path is present and unit-tested
(`project/drugclip_freeze/`: 9/9 tests), but Model A inference requires the `drugclip-cpu`
environment.

### Remediation (human decision required)

- **A** authorize a change to `environment_drugclip_cpu.yml` (e.g. a pytorch version whose CPU runtime
  installs today) - changes environment identity;
- **B** supply a pre-built `drugclip-cpu` environment (conda-pack / micromamba export) that realises the
  frozen spec;
- **C** run Stage 3 on a host where the pinned pytorch 2.0.1 CPU runtime is available.

In all cases the checkpoint, repository commits and the frozen apply-only code path must remain as
verified here.

---

## Preservation

The frozen PACER-FKG modules and `V02_FREEZE_MANIFEST.json` are untouched. Nothing was committed,
tagged or pushed. All run outputs live under `project/results/` (gitignored). The two third-party
checkouts `project/tools/DrugCLIP` and `project/tools/Uni-Core` are untracked and must never be committed.


---

# STAGE 3 RESULT - STAGE3_DRUGCLIP_COMPLETE (2026-10-02, resumed run)

## Runtime compatibility repair

| Item | Value |
|---|---|
| Frozen environment spec | `project/environment_drugclip_cpu.yml` (sha256 `24781f5e90494f2d8e48239defa6ffa658f3e90c79bec9a526a6f2943ad3ecd9`, **not modified**) |
| Original failure | `libtorch_cpu.so: undefined symbol: iJIT_NotifyEvent` |
| Authorized override | **mkl = 2024.0.x** (single override) |
| PyTorch | **2.0.1** (unchanged) |
| Result | `import torch` PASS, torch 2.0.1, cuda False |
| Extra solver flag required | `--channel-priority=flexible` (the frozen spec is unsolvable under default priority) |

Environment record: `stage3/ENV_EXPLICIT_lock.txt`, `stage3/ENV_EXPORT_final.yml`, python 3.10.21,
unicore 0.0.1. Frozen identities verified at execution: DrugCLIP `7a3a3fa3…`, Uni-Core
`44f6386f…`, checkpoint 1,183,713,459 bytes / `dc2c76d0…`.

## Compatibility smoke - PASS

torch import · torch==2.0.1 · DrugCLIP imports · checkpoint loaded ·
`missing_checkpoint_keys == []` · `unexpected_checkpoint_keys == []` ·
1 molecule x 3 frozen M4 pockets finite (score range -0.1895 .. -0.1294) · no optimizer/training path.
Receipt: `stage3/STAGE3_RUNTIME_COMPATIBILITY_RECEIPT_v01.json`.

## Model A - frozen official apply-only

200/200 Stage-2 candidates scored, device cpu, `missing_checkpoint_keys == []`,
`unexpected_checkpoint_keys == []`, all three frozen M4 pockets (7TRQ/7TRP/7TRS),
finite scores spanning -0.2816 .. 0.2567. Archive `stage3/embeddings.npz`.
No candidate selected or modified during Model-A execution.

## Model B - frozen projection-only

Three deployed famaug seeds (20260925/20260926/20260927) applied to Model A's frozen 512-D
representations; official max pooling over the declared M4 pockets. No backbone fitting.

## Dual-model decision - `old_rank <= 50 AND new_rank <= 25`

| | |
|---|---|
| Stage-2 input molecules | 200 |
| old-top50 | 50 |
| new-top25 | 25 |
| **intersection** | **5** |
| fusion used | **no** |
| deterministic candidate ids | **true** (PACER0001..PACER0200) |

**Dual-model shortlist** (ordered by the Model-A decision line):

| # | candidate_id | Model-A rank | Model-B rank |
|---|---|---|---|
| 1 | PACER0108 | 21 | 14 |
| 2 | PACER0157 | 24 | 6 |
| 3 | PACER0080 | 25 | 18 |
| 4 | PACER0149 | 35 | 4 |
| 5 | PACER0089 | 38 | 15 |

All forbidden operations recorded false: PocketRealign fitting, adapter training, fine-tuning,
optimizer creation, backward passes, weight updates, threshold tuning, score fusion.

## Output SHA256

```
stage3/embeddings.npz                       f0a71918b707c13053444a7de9f41d47487201b92f8b79f5ccc5940f1fdcef3a
STAGE3_OUTPUTS/model_a_scores.csv           a263a866db1e4ccdf6ac05d83f81cae3d29704d2ba00cb0f13ee67d538ee05ba
STAGE3_OUTPUTS/model_b_scores.csv           4323130ee12c2774c8ba217bcd5235360e9e7f0ad288116b15db0350867c3364
STAGE3_OUTPUTS/model_b_scores.json          faac31c6a87b9979000590bc274bedbf64cfd648a0ba47eb755bf4975a4b1ca6
STAGE3_OUTPUTS/STAGE3_DUAL_MODEL_DECISION_v01.json   d8a68aeedd8a50038f1acacea821468d9019573f40a128a2ffaad6e5b1c631c2
STAGE3_OUTPUTS/STAGE3_DUAL_MODEL_SHORTLIST_v01.csv   10fc0c8ed474a309271b0817bff38e2080847d3ca0a7272f6f9313a8fe7265ac
STAGE3_DRUGCLIP_RECEIPT_v01.json            502cd42d13820317c324d3bd602e645e3501a8ebdaed6e11563195a5d6c78474
```

Full manifest (33 entries): `STAGE3_SHA256_MANIFEST_v01.txt`.

**Claim boundary:** frozen DrugCLIP apply-only binding-compatibility evidence plus the frozen
dual-model prioritization rule only. Not a PAM classifier, not potency, not experimental confirmation.

## HARD STOP

Docking, IFP, diversity selection, MD and PACER-FKG were **not** executed.
