# Data and artifact provenance

## Frozen screening records

The files in `data/frozen/` were copied from the final project run on 2026-10-04. `manifest.json` records their byte sizes and SHA-256 hashes. The package verifies canonical molecular graphs when joining Module 2 and Module 3, rather than trusting reused display identifiers.

## DrugCLIP

The model used for the final Module 2 ranking is the three-seed GPCR-adapted projection ensemble in `models/used_module2/`. Its training and validation record is `logs/module2_training_validation.json`.

The M4 leave-one-target-out 2023 experiment and the 2026 family-augmentation experiment are retained under `models/reference_not_used/`. They were not used to produce the frozen top-200 ranking. Their benchmark records are retained under `logs/`.

DrugCLIP source: `bowen-gao/DrugCLIP`, archived project commit `7a3a3fa33673f8668c811790f2e4681c98af44ef`. The upstream base checkpoints and representations are not redistributed because of size and upstream terms.

## Structures and docking

M4 experimental coordinates were obtained from RCSB PDB. The six structural channels are Glide-PDB, Glide-BEmin, Glide-BEavg, Vina-PDB, Vina-BEmin and Vina-BEavg. Glide is commercial and is not redistributed. Vina is external open-source software. The package retains the six raw score columns and recomputes their rank fusion.

## Four-context dynamics

The four systems are apo, ACh only, candidate only, and candidate plus ACh. The compact package includes the frozen analysis record but not the large prepared membrane systems or trajectories. A missing trajectory is never treated as a negative functional result.

## Limitations

The exact public-database release identifiers, raw console logs, original hardware and wall-clock time were not retained for every historical run. These gaps are disclosed and are not reconstructed retrospectively.
