# DATA DICTIONARY — CHRM4 PAM Final Modeling Dataset v1.0

**FROZEN 2026-08-16** from Freeze Candidate 2 (FC2).

> This version is scientifically identical to Freeze Candidate 2. No
> scientific rows or modeling labels changed during final freeze.

This dictionary is scientifically identical to the FC2 dictionary; only the
dataset name and freeze status differ.

All CSV files are UTF-8, comma-delimited, header row, `\n` line endings.

## Train/test split safety (applies to every modeling view)

- **ROW-LEVEL RANDOM TRAIN/VALIDATION/TEST SPLIT IS PROHIBITED.**
- **Minimum isolation unit: canonical_molecule_id.**
- **The same canonical molecule must never appear in multiple train,
  validation, or test partitions.**
- For prospective structure-based evaluation, **scaffold-based splitting is recommended**
  (cluster by canonical scaffold before assigning partitions).
- This frozen dataset does **NOT** generate train/test splits.

## activities.csv (one row per CHRM4 source activity; 2115 rows)

Preserved source-fact / normalized v2.1 fields:

| column | meaning |
|---|---|
| `activity_id` | stable activity identifier |
| `parent_id` | frozen deduplicated parent entity id (NOT molecule identity) |
| `source_id` | source identifier (PMID_…, CHEMBLDOC…, SRC…) |
| `assay_id` | assay registry id (blank for narrative rows) |
| `receptor_species` | as-recorded species annotation |
| `assay_type` / `readout` / `cell_system` / `orthosteric_agonist` / `agonist_conc_or_ecxx` | assay context |
| `endpoint_type` / `endpoint_relation` / `potency_nM` / `pPotency` / `endpoint_value` / `endpoint_unit` | quantitative facts |
| `efficacy_pct` / `fold_shift` / `modulation_direction` / `intrinsic_activity` | efficacy facts |
| `legacy_pharmacological_class` | frozen legacy label; PRESERVED AS SOURCE FACT, NEVER modeling truth |
| `quantitative_validity` | exact / censored |
| `evidence_tier` | CORE / EXTENDED / EXPLORATORY |
| `v21_modeling_role`, `v21_exclusion_reason` | frozen v2.1 modeling role (renamed from `modeling_role`) |
| `v21_human_review_*` | frozen v2.1 embedded review fields (frozen baseline; NOT authoritative) |
| `raw_text_value` / `raw_standard_text_value` / `raw_activity_comment` | raw text semantics |
| `qualitative_result` / source / confidence | qualitative result |
| `evidence_supported_role` / `evidence_supported_role_basis` | v2.1 evidence role |
| `compound_id` / `chembl_molecule_id` | structure identifiers |
| `source_location`, `activity_comment`, `notes` | provenance text |

Derived (current-overlay) fields:

| column | meaning |
|---|---|
| `canonical_molecule_id` | canonical molecule (CM…) from free-base InChIKey |
| `canonical_inchikey` | canonical free-base InChIKey |
| `assay_family` | calcium / cAMP / G-protein / ERK / beta-arrestin / GIRK-thallium / … |
| `assay_context` | readout\|cell\|agonist\|ECxx |
| `final_target_valid` | yes (all rows in this CHRM4-scoped dataset) |
| `final_species_valid` | yes / no / unresolved |
| `final_identity_valid` | yes / no |
| `final_direction` | PAM / no_potentiation / unresolved / not_determined / … |
| `final_quantitative_type` | exact / censored / narrative_exact / not_quantitative / … |
| `final_modeling_role` | exact_potency / explicit_inactive / censored_weak_pam / threshold_only_bound / not_determined / unresolved_primary_review_nonblocking / excluded_nonhuman_misassigned / excluded_species_misassigned / identity_conflation_excluded / unresolved_contradictory_censored / identity_failure_no_canonical_structure / not_quantitative_unresolved / censored_pam_potency / efficacy_only / auxiliary / fold_shift / NAM_potency / quantitative_potency_other / extended_exact_secondary / **divergent_duplicate_excluded_exact** / **divergent_duplicate_rejected_representation** |
| `final_exclusion_reason` | why excluded (blank if not excluded) |
| `potency_regression_eligible` | yes / no |
| `pam_binary_eligible` | yes / no |
| `threshold_10uM_eligible` | yes / no |
| `unit_human_decision` | unit decision applied |
| `row_human_decision` | row decision applied |
| `identity_decision` | identity decision applied |
| `decision_provenance` | overlay files consumed for this row |
| `fc2_divergent_group_id` | FC1 divergent measurement group (MG00039/MG00067/MG00260/MG00290/MG00461) for the ten reconciled rows, else blank |
| `fc2_divergent_classification` | A_correct_value_identified / B_distinct_measurements_should_split / C_unresolved_conflicting_duplicate for the ten reconciled rows, else blank |
| `fc2_divergent_disposition` | split_into_two_measurements / exclude_group_from_exact_regression / keep_single_correct_measurement for the ten reconciled rows, else blank |

## compounds.csv (canonical molecule level)

| column | meaning |
|---|---|
| `canonical_molecule_id` | CM… id |
| `canonical_inchikey` / `canonical_smiles` / `canonical_inchi` | canonical (free-base) identity |
| `identity_status` | confirmed / confirmed_structure / conflated_member_excluded / unresolved_structure |
| `canonical_representative_parent_id` | deterministic representative parent |
| `all_parent_ids` / `all_chembl_ids` / `all_names` / `all_aliases` | aggregated member identifiers |
| `identity_decision` / `identity_decision_provenance` | identity overlay trace |
| `member_notes` | per-member chembl/class |
| `modeling_eligibility_note` | e.g. PARENT0872 exclusion note |
| `collision_report` | parent-ID collisions resolved to one canonical structure |

Unresolved-structure parent entities are appended with empty
`canonical_molecule_id` and `identity_status = unresolved_structure`.

## assays.csv / sources.csv

Registry context plus activity counts. Unknown assay metadata is left blank
(never silently filled).

## modeling_potency_exact.csv (combined exact view)

Exact human-M4 PAM EC50 regression view. Columns: `activity_id`,
`measurement_group_id`, `fc1_measurement_group_id`, `replicate_count`
(FC1-compatible alias), `representation_count`, `measurement_representation_type`
(single / rounding_or_pEC50_conversion_duplicate / distinct_measurements_split /
divergent_duplicate_resolved_single_value), `aggregation_method`
(log_mean_pEC50 / single), `aggregated_pEC50`, `aggregated_EC50_nM`,
molecule/assay/potency fields, `pEC50`, eligibility and decision provenance,
and divergent-duplicate classification/disposition for the reconciled rows.
Measurement groups require identical canonical molecule + source + assay +
endpoint + species + assay context, except the MG00039 split (two deliberately
separate groups). **Combined exact view: maximum coverage; assay_family must be
considered.** Row-level random splitting is prohibited; minimum isolation unit
is `canonical_molecule_id`.

## modeling_potency_exact_calcium.csv (calcium-restricted exact view)

Identical schema and scientific gates to `modeling_potency_exact.csv`,
restricted to `assay_family == calcium` (dominant harmonized calcium-
mobilization exact PAM subset). **Calcium-restricted exact view: preferred more
assay-harmonized QSAR view where appropriate.** RU003 is absent only because it
is GIRK/thallium; it remains in the combined view.

## modeling_pam_vs_inactive.csv

Strict binary classification. `class_label` = positive / negative.
Negatives carry `negative_evidence_type` (A_primary_confirmed_inactive /
B_database_text_confirmed_inactive), `negative_evidence_tier` and
`negative_evidence_provenance_group`. `molecule_strict_label` =
positive / negative / excluded_conflicted. The six divergent-duplicate rows are
positive (`positive_evidence_type =
human_M4_PAM_direction_confirmed_exact_value_unresolved`); the rejected MG00461
representation is absent.

## modeling_threshold_10uM.csv

`class_label` = threshold_positive / threshold_negative.
`endpoint_relation_recorded` vs `endpoint_relation_human_resolved` distinguish
frozen encoding from human-resolved right-censored interpretation.
Threshold-negatives are never called inactive. ACT_18852473 is threshold-negative
(stored value exactly 10000 nM; encoded threshold; one of the 17 EXTENDED
contribution; not an additional 78th source category).

## modeling_extended.csv

Non-strict scientifically useful evidence. `extended_category` describes the
evidence type; `strict_view_membership` lists which strict views the row also
appears in; `species_caveat` flags non-human narrative values (RU001/RU002).
The six excluded divergent rows and the one rejected MG00461 representation are
retained here for provenance.

## excluded_or_unresolved.csv

Every excluded/unresolved activity row:
`activity_id`, `parent_id`, `source_id`, `assay_id`, `canonical_molecule_id`,
`canonical_inchikey`, `exclusion_category`, `final_modeling_role`,
`exclusion_reason`, `decision_provenance`, `molecule_survives_elsewhere`,
`notes`. `molecule_survives_elsewhere` distinguishes strict/extended survival.

## manifest_sha256.json

SHA-256 manifest of all frozen outputs, the FC2 builder/validator/test file,
the authoritative scientific inputs (normalized v2.1 activity file, human
decision overlays, negative-set purity audit, divergent-duplicate
reconciliation), and the FC2 source manifest/hash. Deterministic (no
timestamps).
