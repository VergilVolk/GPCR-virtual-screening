# Data notes

## Packaged example

- `candidate_library.csv`: 12 generated M4 candidates selected from the frozen 2,605-molecule library.
- `m4_candidate_representations.npz`: DrugCLIP molecule representations for those 12 structures.
- `m4_pocket_representations.npz`: the matching M4 pocket representation.
- `ensemble_docking_summary.csv`: ten-conformation Vina summary and structural quality gates.
- `four_context_validation.csv`: frozen retrospective controls for the dynamic endpoint.
- `prospective_four_context_summary.csv`: honest summary of the three prospective candidates.

## Provenance

Candidate structures were generated and filtered within this project. Protein structures were obtained from the Protein Data Bank. Public ligand annotations and benchmark labels were derived from their cited primary sources and public databases. Upstream DrugCLIP representations follow the DrugCLIP release and license. AutoDock Vina was used under its open-source license.

Raw PDB structures, MD trajectories and public benchmark archives are not duplicated here. Their source identifiers, preparation scripts and checksums are retained in the complete repository. No hidden test labels were used to tune the packaged candidate ranking.

## Leakage controls

The binding adapter uses an M4 leave-one-target-out protocol. Evaluation targets are excluded from adapter fitting, and chemical-series checks are reported separately. MD trajectory windows are repeated observations from the same molecule and are never counted as independent compounds.
