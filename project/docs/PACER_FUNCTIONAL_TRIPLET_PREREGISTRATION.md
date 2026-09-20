# PACER Functional Triplet v0.1 preregistration

Frozen before the first model run on 2026-08-30.

## Question

Does functional hard-negative metric learning improve new-series M4 PAM potency ranking beyond the frozen PACER-FS centered-LightGBM baseline and the same neural regressor without triplet supervision?

## Outer protocol

- Data: 430 potency molecules.
- Outer unit: complete `source_component`; evaluate the same 11 components with `n>=12`.
- The held-out component contributes no labels to training or triplet construction.
- Three anchors are chosen at the 25/50/75 percentiles of an outer-training-only absolute LightGBM prediction.
- Metrics use the remaining query molecules only.
- The three labels calibrate only a constant pEC50 offset; they cannot alter rank.

## Frozen models

1. Absolute-LightGBM + anchor offset.
2. PACER-FS centered-LightGBM + anchor offset.
3. MLP-Reg: centered potency regression, no metric loss.
4. MLP-RandomTriplet: the same MLP and positive pairs, with same-series non-positive random negatives.
5. PACER-FM-HardTriplet: the same MLP with same-series, chemically close potency cliffs as hard negatives.

Neural models use ECFP4 plus nine training-standardized descriptors, a 256-64 encoder, 32-dimensional normalized metric embedding, 160 full-batch epochs, AdamW 1e-3, weight decay 1e-4, cosine-triplet margin 0.2, metric-loss weight 0.5, and three fixed seeds 17/29/43. No hyperparameter is selected from outer results.

## Pair definitions

- Positive: same training `source_component`, Tanimoto >=0.50, absolute delta pEC50 <=0.30.
- Hard negative: same training `source_component`, Tanimoto >=0.55, absolute delta pEC50 >=1.0.
- Random negative control: same component and absolute delta pEC50 >0.30, without the hard-negative similarity/cliff requirement.

## Endpoints

- Primary: macro within-series Spearman, paired series-bootstrap delta versus PACER-FS.
- Secondary: worst-series Spearman, macro MAE after three-anchor offset, and held-out query cliff direction accuracy.
- Internal Go: macro Spearman above PACER-FS, paired 95% CI lower bound >0, worst-series no worse, and cliff direction accuracy above MLP-Reg.
- Promotion remains false until an untouched external campaign supports the frozen model.

Triplets are constructed after every outer split. A positive result is internal development evidence, not external SOTA or confirmation of generated PAMs.
