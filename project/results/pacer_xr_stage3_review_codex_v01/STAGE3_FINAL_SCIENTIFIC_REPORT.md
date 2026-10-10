# PACER-XR Stage 3 final scientific report

## Conclusion

The frozen Cascade is supported as a reproducible, traceable ensemble-prioritization strategy, but its top end has clear single-condition dependence. Four Cascade Top10 candidates (PACERGEN00123, PACERGEN00462, PACERGEN02015, PACERGEN02355) fall to pure six-channel ranks 100, 104, 122 and 66. The pure-XR Top3 (PACERGEN01495, PACERGEN01580, PACERGEN01593) have balanced six-channel percentile profiles and complete representative pose contact coverage. Neither observation predicts PAM activity without external validation.

## Integrity and coverage

- Verified frozen ranking SHA256: `6fff7e34d94466ff288b52281b1c3a64fef4ea6cd9da84ca45bc703429a04559`.
- 20 Cascade candidates were analysed. The three pure-XR candidates are within that set at Cascade 10-12.
- Glide: 1 expected/stored requests, 39 CONTACTS_EXTRACTED statuses, 472 contact rows from 5 Poseviewer files. One of 40 intended Glide contexts is unavailable: PACERGEN01954 PDB/7TRS.
- Vina: 40/40 PDB/BEmin contexts and 482 contact rows. In total, 79 real pose-context contact sets were used; no missing pose was fabricated.

## Structural and scoring evidence

- Contact records are <=4.0 A heavy-atom geometries. They repeatedly locate PDB poses in the stated 7TRS pocket region, but do not establish hydrogen bonds, affinity, mechanism, or pose identity.
- PDB versus ensemble scoring is compared with frozen channel percentiles, not raw energy. Ensemble residue numbers are not mapped to 7TRS numbers.
- PACERGEN02015 additionally has Glide/Vina mean divergence +0.217. PACERGEN02432 is ensemble-enhanced (E-PDB +0.357) but still pure-XR #10.

## Recommendation for discussion

Retain the frozen Cascade unchanged for study continuity, and present its Top1% head with explicit gate/ensemble-dependence annotations. Use pure-XR Top3 as structurally better-balanced comparators. Require independent pose assessment and experimental functional or biophysical data before claims of prospective ranking performance or any revision of the frozen strategy.

## Files and provenance

- STAGE3_CANDIDATE_COMPARISON.md: Top20 candidate-level scores, raw values, receptor contexts, and all contact records.
- XR_TOP3_STRUCTURAL_REVIEW.md: detailed pure-XR Top3 review.
- CASCADE_STRATEGY_ASSESSMENT.md: all supplied ablation scenarios and interpretation.
- candidate_evidence_table.csv: machine-readable candidate evidence.
- Inputs: full200 ranking, raw scores, Glide/Vina channel tables, fusion audit, and the supplied ablation diagnostics, rank, pose, contact, status, request and summary files.
