# FINAL FREEZE RECORD — CHRM4 PAM Final Modeling Dataset v1.0

**Freeze date: 2026-08-16**

## 1. Freeze decision

The dataset **CHRM4 PAM Final Modeling Dataset v1.0** is declared the official
frozen modeling dataset of the CHRM4 PAM activity-evidence pipeline.

**This version is scientifically identical to Freeze Candidate 2.**
**No scientific rows or modeling labels changed during final freeze.**

The freeze is an administrative/provenance freeze of
`CHRM4_PAM_FINAL_DATASET_FC2/` (FC2). No scientific interpretation changed; no
modeling CSV row, label, activity, potency, canonical identity, measurement
grouping, exclusion, or eligibility changed; no FC3 was created.

Source freeze candidate: `CHRM4_PAM_FINAL_DATASET_FC2/`.

Freeze procedure:

1. Byte-identical copy of the 10 scientific CSVs from FC2 into
   `CHRM4_PAM_FINAL_DATASET_v1.0/`.
2. Frozen documentation (`README.md`, `DATA_DICTIONARY.md`,
   `FREEZE_REPORT.md`, this record).
3. Final manifest (`manifest_sha256.json`) recording SHA-256 of every frozen
   output, the authoritative scientific inputs, the FC2 builder/validator/test
   file, and the FC2 source manifest/hash.
4. Final validation (byte-identity + invariant verification + FC2 validator +
   FC2 tests).

## 2. Modeling safety rules (mandatory, apply to every downstream use)

- **ROW-LEVEL RANDOM TRAIN/VALIDATION/TEST SPLIT IS PROHIBITED.**
- **Minimum isolation unit: canonical_molecule_id.**
- **The same canonical molecule must never occur in multiple partitions.**
- For prospective structure-based evaluation, **scaffold-based splitting is recommended**.
- For potency modeling:
  - `modeling_potency_exact_calcium.csv` = **preferred assay-harmonized QSAR
    starting view**;
  - `modeling_potency_exact.csv` = broader combined exact-potency view;
    **assay-family/context covariates must be retained**.
- **Do not describe representation duplicates as biological replicates.**

## 3. Headline counts (identical to FC2)

| metric | count |
|---|---|
| Exact potency source representations | **712** |
| Exact potency modeling measurements | **472** |
| Exact potency canonical molecules | **460** |
| Calcium exact modeling measurements | **437** |
| Calcium exact canonical molecules | **430** |
| Strict PAM positive canonical molecules | **463** |
| Strong inactive canonical molecules | **66** |
| 10uM threshold-positive canonical molecules | **463** |
| 10uM threshold-negative canonical molecules | **70** |
| Excluded-only canonical molecules | **118** |
| Extended-only canonical molecules | **155** |
| Surviving canonical molecules with excluded rows | **53** |
| Unresolved-structure parents | **51** |

## 4. Scientific freeze basis — final five divergent-group resolutions

| group | classification | final resolution |
|---|---|---|
| MG00039 | `B_distinct_measurements_should_split` | split into two distinct measurements: **410.0 nM** (ACT_16608408) and **1202.26 nM** (ACT_16608432, pEC50 5.92) |
| MG00067 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression** |
| MG00260 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression** |
| MG00290 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression** |
| MG00461 | `A_correct_value_identified` | **103.0 nM retained as the single supported measurement** (ACT_18137897, pEC50 6.9872) |

No divergent pair is averaged; no pair is treated as a biological replicate.
The six excluded rows and the one rejected representation remain traceable in
`activities.csv`, `modeling_extended.csv`, and `excluded_or_unresolved.csv`.

**The 3 excluded exact-value molecules (MG00067 / MG00260 / MG00290) remain
valid PAM positives / 10 µM positives** where direction/threshold class is
scientifically stable.

## 5. Known exclusions (frozen)

- **RU001**: rat 400 nM must not be used as human exact potency.
- **RU002**: species-conflated 629 nM must not be used as human exact potency.
- **RU011**: unresolved primary-review boundary; excluded from strict views.
- **RU018 / PARENT0872**: identity-conflated narrative entity; excluded.
- **PARENT0805**: 13l / VU6008810 identity confirmed.
- **PARENT0820**: 14o / VU6008677 identity confirmed.

Strict-leak gates verified = 0 for RU001, RU002, RU011, PARENT0872.

## 6. Validation summary

- scientific CSVs byte-identical to FC2: **10/10 verified**
- exact measurements = 472 / exact molecules = 460: verified
- calcium measurements = 437 / calcium molecules = 430: verified
- strong inactive molecules = 66 / threshold-negative molecules = 70: verified
- RU001 / RU002 / RU011 / PARENT0872 strict leaks = 0: verified
- no divergent group incorrectly averaged: verified
- canonical structure deduplication preserved: verified
- representation terminology preserved: verified
- all hashes resolve: verified
- deterministic freeze procedure: verified (fresh FC2 builder rebuild is
  byte-identical for all scientific CSVs)
- FC2 validator (`scripts/validate_chrm4_final_dataset_fc2.py`): **PASS**
- FC2 tests (`tests/test_chrm4_final_dataset_fc2.py`): **PASS**
- final v1.0 validation (`scripts/validate_chrm4_final_dataset_v1_0.py`): **PASS**

## 7. Provenance

- Source freeze candidate: `CHRM4_PAM_FINAL_DATASET_FC2/`
- FC2 builder: `scripts/build_chrm4_final_dataset_fc2.py`
- FC1 builder (reused logic): `scripts/build_chrm4_final_dataset.py`
- FC2 validator: `scripts/validate_chrm4_final_dataset_fc2.py`
- FC2 test file: `tests/test_chrm4_final_dataset_fc2.py`
- Normalized v2.1 activity input:
  `targets/CHRM4-PAM/normalized/v2_1/activities_v2_1.jsonl`
- Human overlays: `human_audit/label_blind_primary_review_001/`
  (`human_decisions.jsonl`, `human_row_decisions.jsonl`,
  `human_identity_decisions.jsonl`, `human_decision_reassessments.jsonl`,
  `negative_set_purity_audit.csv`, `fc1_divergent_duplicate_reconciliation.csv`)
- Final manifest: `manifest_sha256.json` (records all hashes, including the FC2
  source manifest/hash).

---

STATUS: FROZEN
DATASET: CHRM4 PAM Final Modeling Dataset v1.0
SOURCE CANDIDATE: FC2
SCIENTIFIC CHANGES FROM FC2: NONE
DATE: 2026-08-16
