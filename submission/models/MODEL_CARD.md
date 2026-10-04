# DrugCLIP projection models

## Model used in the frozen Module 2 run

`used_module2/` contains the three M4-held-out projection checkpoints that reproduce `results/02_drugclip_top200.csv`.

- backbone space: archived DrugCLIP 2023 representations
- task: transfer to the held-out M4 pocket from B2AR, CCR2 and M2R
- adaptation: rank-8 residual projection updates, three seeds
- seeds: 20260925, 20260926 and 20260927
- ensemble: arithmetic mean of the three cosine scores

The checkpoints estimate molecule-pocket compatibility. They do not predict PAM efficacy.

## Reference models not used for the frozen ranking

`reference_not_used/gpcr_retrieval_hard_triplet/` contains a separate four-GPCR CE plus hardest-pocket triplet experiment.
`reference_not_used/drugclip2026_family_aug/` contains the 2026 family-augmentation experiment.
Both are supplied for audit and comparison only. Neither directory contributed to `02_drugclip_top200.csv` or the final candidate order.

The large upstream DrugCLIP base checkpoints are not redistributed in this package. Their local copies were about 1.18 GB each and exceed GitHub's file-size limit. The compact projection checkpoints, hashes and validation reports are included.

## Known limits

The package includes the 2,605 molecule representations and M4 pocket representation needed to recompute Module 2. The original hardware and wall-clock time were not retained and are not reconstructed.
