# PACER prospective run 2026-10-02 (Stages 1-3)

Authorized scope: **Stage 1 generation -> Stage 2 physchem/drug-likeness filtering ->
Stage 3 frozen DrugCLIP dual-model screening**. Strictly sequential and fail-closed.
No docking, IFP, multi-conformation analysis, diversity selection beyond the frozen
DrugCLIP decision, MD preparation, MD or PACER-FKG.

## Production baseline

| | |
|---|---|
| Branch | `integration/final-freeze-v01` |
| Commit | `ddc11f970f508d89a75dac8ca30c41cfcd361731` |
| Tag | `pacer-code-freeze-v01` |
| PACER-FKG freeze anchor | `b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e346e11bd` |

## Directory layout

| Path | Contents |
|---|---|
| `STAGE1_RUN_SPEC_v01.json` | run spec, frozen before execution, incl. the PYTHONHASHSEED pin and the defect log |
| `STAGE1_SEED_DERIVATION_RECEIPT_v01.json` | seed derivation from tracked frozen data |
| `STAGE1_GENERATION_RECEIPT_v01.json` | Stage-1 deliverable |
| `STAGE1_QC_v01.json` | Stage-1 QC gates |
| `STAGE2_QC_v01.json` | Stage-2 QC + exact rejection accounting |
| `STAGE2_FILTERING_RECEIPT_v01.json` | Stage-2 deliverable |
| `STAGE3_DRUGCLIP_RECEIPT_v01.json` | Stage-3 deliverable |
| `stage3/` | Model A run inputs/outputs (LMDBs, embeddings.npz) |
| `STAGE3_OUTPUTS/` | model_a_scores.csv, model_b_scores.csv, dual-model decision, shortlist |
| `tools/` | run scaffolding and QC scripts |
| `logs_*.txt` | raw command output |

## Declared deviations and environment provisioning

1. **PYTHONHASHSEED pin.** The frozen generator is non-deterministic under default hash
   randomisation (two unpinned runs gave 3212 vs 930 passing molecules). `PYTHONHASHSEED=0`
   is declared in the run spec as the effective RNG pin. No scientific parameter changed.
2. **lightgbm installed** into conda env `chrm4_vs` with `--no-deps` (absent from every env;
   the frozen Stage-2 script imports it at module level). No existing package changed; the
   affected model is reference-only and not read by the eligibility filter.
3. **bzip2 absent in WSL** - micromamba was installed through the existing conda instead of
   the `.tar.bz2` route used by the frozen setup script.
4. **Frozen shell script is a CRLF checkout** - a CRLF-stripped copy is executed; both hashes
   are recorded. The tracked file is not modified.
5. **git safe.directory** is set because the cloned third-party repositories are owned by the
   Windows user while WSL runs as root.
