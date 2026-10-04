# Data and artifact provenance

## Frozen screening records

The files in `data/frozen/` were copied from the final project run on 2026-10-04. `manifest.json` records their byte sizes and SHA-256 hashes. The package verifies canonical molecular graphs when joining Module 2 and Module 3, rather than trusting reused display identifiers.

## DrugCLIP

The model used for the final Module 2 ranking is the three-seed 2023 M4-held-out projection ensemble in `models/used_module2/`. This identity is confirmed by forward inference over all 2,605 packaged molecule representations: its scores and top-200 order reproduce the frozen table.

The four-GPCR CE plus hardest-pocket triplet experiment and the 2026 family-augmentation experiment are retained under `models/reference_not_used/`. They were not used to produce the frozen top-200 ranking.

DrugCLIP source: `bowen-gao/DrugCLIP`, archived project commit `7a3a3fa33673f8668c811790f2e4681c98af44ef`. The upstream base checkpoints and representations are not redistributed because of size and upstream terms.

## Structures and docking

M4 experimental coordinates were obtained from RCSB PDB. The six structural channels are Glide-PDB, Glide-BEmin, Glide-BEavg, Vina-PDB, Vina-BEmin and Vina-BEavg. Glide is commercial and is not redistributed. Vina is external open-source software. The package retains the six raw score columns and recomputes their rank fusion.

## Four-context dynamics

The four systems are apo, ACh only, candidate only, and candidate plus ACh. The compact package includes the frozen analysis record but not the large prepared membrane systems or trajectories. A missing trajectory is never treated as a negative functional result.

The upstream Geom2Vec ViSNet checkpoint used for frozen feature extraction is bundled at `project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth` (SHA-256 `b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`, upstream commit `371d642ec1061664f16e49fcac702d07fc8d0b51`). It was not trained by this project.

## Limitations

The exact public-database release identifiers, raw console logs, original hardware and wall-clock time were not retained for every historical run. These gaps are disclosed and are not reconstructed retrospectively.
