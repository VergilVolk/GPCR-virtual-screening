# Data notes

## Packaged example

- `candidate_library.csv`: 12 generated M4 candidates selected from the frozen 2,605-molecule library.
- `m4_candidate_representations.npz`: DrugCLIP molecule representations for those 12 structures.
- `m4_pocket_representations.npz`: the matching M4 pocket representation.
- `six_channel_docking_summary.csv`: complete Glide/Vina PDB, BEmin and BEavg scores for the three candidates entering the functional stage.
- `six_channel_provenance.json`: source branch, channel contract and missing-score policy.
- `four_context_validation.csv`: frozen retrospective controls for the dynamic endpoint.
- `prospective_four_context_summary.csv`: honest summary of the three prospective candidates.

## Provenance

Candidate structures were generated and filtered within this project. Protein structures were obtained from the Protein Data Bank. Public ligand annotations and benchmark labels were derived from their cited primary sources and public databases. Upstream DrugCLIP representations follow the DrugCLIP release and license. AutoDock Vina was used under its open-source license.

Raw PDB structures, MD trajectories and public benchmark archives are not duplicated here. Their source identifiers, preparation scripts and checksums are retained in the complete repository. No hidden test labels were used to tune the packaged candidate ranking.

The six candidate channels were produced by the production PACER-XR workflow on branch `backup/pacer-xr-emergency-8mol-v01`. The packaged rows correspond to PACER0010, PACER0027 and PACER0073 by exact canonical SMILES. All three Glide and all three Vina channels are present; no channel was substituted or imputed.

## Leakage controls

The binding adapter uses an M4 leave-one-target-out protocol. Evaluation targets are excluded from adapter fitting, and chemical-series checks are reported separately. MD trajectory windows are repeated observations from the same molecule and are never counted as independent compounds.
