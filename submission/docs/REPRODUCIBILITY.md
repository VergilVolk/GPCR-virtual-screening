# Reproducibility

## Frozen result replay

```bash
python predict.py
python verify_submission.py
```

This path is CPU-only. It performs Module 2 forward inference on 2,605 molecules, verifies the top-200 order, checks canonical molecular identity, recomputes six-channel fusion, re-exports four-context evidence, and compares the final table with `results/results.csv`.

## Model retraining

```bash
python train.py --representations REPRESENTATIONS.npz --projection BASE_PROJECTION.pt --pairs PAIRS.csv --reference-predictions REFERENCE.csv
```

The source code and final compact projection checkpoints are included. The large upstream DrugCLIP checkpoint and archived training representations are not included.

## Full production rerun

A complete rerun additionally requires prepared M4 structures, AutoDock Vina, a licensed Glide installation, membrane MD systems and trajectories. The included production scripts do not substitute one docking engine for another.

## Historical limitation

Artifact hashes, seeds, validation reports and frozen outputs are preserved. Raw console logs, hardware and wall-clock time were not retained for every historical run.
