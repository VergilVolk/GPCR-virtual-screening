# CHRM4 PAM Final Modeling Dataset v1.0 — FREEZE REPORT

**Status: FROZEN** — 2026-08-16

## 0. What this freeze is

This is the **administrative / provenance freeze** of Freeze Candidate 2
(FC2) as the official **CHRM4 PAM Final Modeling Dataset v1.0**.

**This version is scientifically identical to Freeze Candidate 2. No scientific
rows or modeling labels changed during final freeze.**

The freeze procedure was:

1. Byte-identical copy of every scientific CSV from
   `CHRM4_PAM_FINAL_DATASET_FC2/` into `CHRM4_PAM_FINAL_DATASET_v1.0/`
   (`activities.csv`, `compounds.csv`, `assays.csv`, `sources.csv`,
   `modeling_potency_exact.csv`, `modeling_potency_exact_calcium.csv`,
   `modeling_pam_vs_inactive.csv`, `modeling_threshold_10uM.csv`,
   `modeling_extended.csv`, `excluded_or_unresolved.csv`).
2. Frozen documentation (`README.md`, `DATA_DICTIONARY.md`,
   `FREEZE_REPORT.md`, `FINAL_FREEZE_RECORD.md`) reflecting the frozen status.
3. A new final manifest (`manifest_sha256.json`) recording SHA-256 of every
   frozen output, the authoritative scientific inputs, the FC2
   builder/validator/test file, and the FC2 source manifest/hash.
4. Final validation: byte-identity vs FC2, full invariant re-verification,
   FC2 validator and FC2 test suite re-run.

No FC3 was created. No scientific interpretation, CSV row, label, activity,
potency, canonical identity, measurement grouping, exclusion, or eligibility
was changed.

## 1. Headline counts (v1.0, identical to FC2)

| view | count |
|---|---|
| exact source representations (activity rows) | **712** |
| exact modeling measurements (groups) | **472** |
| exact canonical molecules | **460** |
| calcium exact modeling measurements | **437** |
| calcium exact canonical molecules | **430** |
| strict PAM positive canonical molecules | **463** |
| strong inactive canonical molecules | **66** |
| 10uM threshold-positive canonical molecules | **463** |
| 10uM threshold-negative canonical molecules | **70** |
| excluded-only canonical molecules | **118** |
| extended-only canonical molecules | **155** |
| surviving canonical molecules with excluded rows | **53** |
| unresolved-structure parents | **51** |

Additional row-level context (unchanged from FC2): activities.csv = 2115 rows;
CORE activity rows = 868; binary positive rows = 718; strong inactive rows =
67; threshold-positive rows = 718; threshold-negative rows = 78; extended rows
= 1080; excluded/unresolved activity rows = 265; excluded canonical molecules
(≥1 excluded row) = 190; parent-ID collisions resolved = 2;
context-conflicted molecules = 0.

## 2. Exact potency view

- exact source representations (activity rows): **712**
- exact modeling measurements (groups): **472**
  - multi-representation groups: **240** (240
    `rounding_or_pEC50_conversion_duplicate`)
  - single-representation groups: **232** (229 `single`, 2
    `distinct_measurements_split` (MG00039a/MG00039b), 1
    `divergent_duplicate_resolved_single_value` (MG00461))
- exact canonical molecules: **460**
- assay families (rows): GIRK/thallium 2, calcium 665, ERK 24, unspecified 2,
  cAMP 19
- the combined exact view provides maximum coverage; **assay_family must be
  considered** in downstream QSAR.

### Calcium-restricted exact view

`modeling_potency_exact_calcium.csv` = 665 calcium rows, 437 modeling
measurements, 430 canonical molecules. **Preferred assay-harmonized QSAR
starting view** (`assay_family == calcium` only, same FC2 scientific gates).
RU003 (GIRK/thallium) is absent by design (assay family), not by eligibility;
it remains in the combined exact view.

## 3. Scientific freeze basis — final five divergent-group resolutions

| group | classification | disposition |
|---|---|---|
| MG00039 | `B_distinct_measurements_should_split` | **split into two distinct measurements**: `MG00039a` ACT_16608408 = **410.0 nM** (HTS-hit characterization, compound 5 / VU0473619) and `MG00039b` ACT_16608432 = **1202.26 nM** / pEC50 5.92 (resynthesized-batch confirmation). No averaging. |
| MG00067 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression**; both source rows preserved in extended/excluded provenance (no choosing/averaging 159 nM vs pEC50 7.00 ⇔ 100 nM). |
| MG00260 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression**; both source rows preserved (no choosing/averaging 29 nM vs pEC50 7.74 ⇔ 18.2 nM). |
| MG00290 | `C_unresolved_conflicting_duplicate` | **excluded from exact regression**; both source rows preserved (no choosing/averaging 37 nM vs pEC50 7.12 ⇔ 75.86 nM). |
| MG00461 | `A_correct_value_identified` | **103.0 nM retained as the single supported measurement** (ACT_18137897, pEC50 6.9872), corroborated by primary narrative; rejected representation ACT_18137896 (pEC50 7.20 / 63.1 nM) preserved in provenance only (no second measurement, no average). |

None of the five divergent pairs is treated as a biological replicate; no pair
is averaged. The six excluded rows and the one rejected representation remain
fully traceable to their source `activity_id`s in `activities.csv`,
`modeling_extended.csv` and `excluded_or_unresolved.csv`.

**The 3 exact-excluded molecules (MG00067 / MG00260 / MG00290) remain valid PAM
positives / 10 µM positives** where direction/threshold class is scientifically
stable: PAM direction is primary-source confirmed, and all conflicting values
are below the 10 µM threshold.

## 4. Known exclusions (frozen)

- **RU001**: rat 400 nM must **not** be used as human exact potency.
- **RU002**: species-conflated 629 nM must **not** be used as human exact
  potency.
- **RU011**: unresolved primary-review boundary; **excluded from strict views**
  (non-blocking).
- **RU018 / PARENT0872**: identity-conflated narrative entity; **excluded**
  (contributes zero canonical modeling molecules).
- **PARENT0805**: `13l` / VU6008810 identity **confirmed**
  (thieno[3,2-d]pyrimidine), retained.
- **PARENT0820**: `14o` / VU6008677 identity **confirmed**
  (furo[3,2-d]pyrimidine), retained.

## 5. Strict-leak gates (all verified = 0)

- RU001 strict leaks: **0**
- RU002 strict leaks: **0**
- RU011 strict leaks: **0**
- PARENT0872 strict leaks: **0**
- No wrong-species numeric activity enters human-M4 strict potency.
- Canonical structure deduplication is mandatory and preserved.
- The strong inactive set excludes all threshold-only censored rows.
- The rejected MG00461 representation is absent from all strict views.

## 6. Modeling safety rules (frozen)

- **ROW-LEVEL RANDOM TRAIN/VALIDATION/TEST SPLIT IS PROHIBITED.**
- **Minimum isolation unit: `canonical_molecule_id`.**
- **The same canonical molecule must never occur in multiple partitions.**
- For prospective structure-based evaluation, **scaffold-based splitting is
  recommended**.
- `modeling_potency_exact_calcium.csv` = preferred assay-harmonized QSAR
  starting view.
- `modeling_potency_exact.csv` = broader combined exact-potency view;
  assay-family/context covariates must be retained.
- **Do not describe representation duplicates as biological replicates.**
- This frozen dataset ships no train/validation/test partitions.

## 7. Validation performed

1. All 10 scientific CSVs verified **byte-identical** to FC2.
2. Final manifest hashes verified against the frozen files and the
   authoritative inputs.
3. Full invariant re-verification on the frozen files (counts 712/472/460,
   437/430, 463/66, 463/70; divergent-group dispositions; no replicated
   averaging; strict leaks = 0; calcium-only restriction; canonical
   deduplication; representation terminology; excluded semantics).
4. FC2 validator (`scripts/validate_chrm4_final_dataset_fc2.py`) re-run:
   passes.
5. FC2 test suite (`tests/test_chrm4_final_dataset_fc2.py`) re-run: passes.
6. Deterministic freeze procedure: the frozen scientific CSVs are byte-identical
   to a fresh deterministic FC2 builder rebuild (no timestamps; every byte a
   pure function of the immutable inputs, the FC1 builder, and the FC2 builder).

## 8. Determinism and traceability

- builder (source of scientific content): `scripts/build_chrm4_final_dataset_fc2.py`
  (reuses `scripts/build_chrm4_final_dataset.py`)
- FC2 validator: `scripts/validate_chrm4_final_dataset_fc2.py`
- FC2 tests: `tests/test_chrm4_final_dataset_fc2.py`
- final validation: `scripts/validate_chrm4_final_dataset_v1_0.py`
- manifest: `manifest_sha256.json` (hashes of every frozen output, the
  authoritative scientific inputs, FC2 builder/validator/test file, and the
  FC2 source manifest/hash)
- every modeling row traces to a source `activity_id`; no undocumented
  synthetic rows are introduced.

## 9. Freeze conclusion

The dataset is scientifically identical to FC2, all frozen invariants hold,
all hashes resolve, and the source candidate remains fully validated.
**CHRM4 PAM Final Modeling Dataset v1.0 is frozen.**
