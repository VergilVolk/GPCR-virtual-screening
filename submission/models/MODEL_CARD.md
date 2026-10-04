# DrugCLIP projection models

## Model used in the frozen Module 2 run

`used_module2/` contains the three projection checkpoints used to rank the 2,605 generated molecules and produce `results/02_drugclip_top200.csv`.

- backbone space: archived DrugCLIP 2023 representations
- task: four-GPCR pocket retrieval (B2AR, CCR2, M2R and M4R)
- adaptation: rank-4 projection updates, three seeds
- seeds: 20260925, 20260926 and 20260927
- ensemble: arithmetic mean of the three cosine scores

The checkpoints estimate molecule-pocket compatibility. They do not predict PAM efficacy.

## Reference models not used for the frozen ranking

`reference_not_used/drugclip2023_m4_loto/` contains an M4 leave-one-target-out experiment.
`reference_not_used/drugclip2026_family_aug/` contains the 2026 family-augmentation experiment.
Both are supplied for audit and comparison only. Neither directory contributed to `02_drugclip_top200.csv` or the final candidate order.

The large upstream DrugCLIP base checkpoints are not redistributed in this package. Their local copies were about 1.18 GB each and exceed GitHub's file-size limit. The compact projection checkpoints, hashes and validation reports are included.

## Known limits

The original hardware and wall-clock time were not retained for every historical training job. This missing metadata is not reconstructed. The package preserves the final checkpoint bytes, seeds, training report and frozen inference outputs.
