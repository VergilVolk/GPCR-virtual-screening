# M4 PAM virtual screening

This package contains the code, frozen inputs, model adapters and result tables for the four-stage Track 3 workflow.

## Workflow

1. M4 PAM candidate generation and medicinal-chemistry filtering.
2. M4 pocket retrieval with a GPCR-adapted DrugCLIP model.
3. Six-channel docking: Glide and Vina on the experimental structure, ensemble best-energy summary (BEmin), and ensemble mean-energy summary (BEavg).
4. Four-context receptor dynamics: apo, ACh, candidate, and candidate plus ACh.

The first three stages prioritize plausible binders. Stage 4 tests reproducible conditional structural responses. None of these calculations alone proves PAM efficacy.

## Quick verification

Tested with Python 3.11 on Windows 11. The frozen replay runs on CPU.

```bash
pip install -r requirements.txt
python predict.py
python verify_submission.py
```

`predict.py` verifies the frozen data hashes, checks molecular identity across stages, recomputes the six-channel rank fusion, exports the three four-context records, and writes `results/results.csv`.

## Results

- `results/01_generated_candidates_2605.csv`: generated valid molecules.
- `results/01_candidate_audit_2592.csv`: deduplicated candidate audit after removal of known reference molecules.
- `results/02_drugclip_top200.csv`: reranked top 200 from the 2,605-molecule library.
- `results/03_six_channel_8.csv`: five leading candidates plus three historical anchors, with all six docking channels.
- `results/04_four_context_3.csv`: four-context evidence for the three anchors.
- `results/results.csv`: standardized eight-candidate evidence table.

The module files are retained separately because the final table is a summary, not a replacement for the intermediate evidence.

## Models

`models/used_module2/` contains the three projection checkpoints used for the frozen Module 2 ranking.

Two comparison sets are also included and explicitly marked as unused:

- `models/reference_not_used/drugclip2023_m4_loto/`
- `models/reference_not_used/drugclip2026_family_aug/`

They did not contribute to the top-200 table or final ranking. See `models/MODEL_CARD.md` and `logs/training_manifest.json` for hashes and scope.

The upstream 1.18 GB DrugCLIP base checkpoints are not redistributed. The compact adapted projections used by this project are included.

## Entrypoints

```bash
# Recheck generated structures and medicinal-chemistry rules
python module1_generate_filter.py --input data/frozen/generated_candidates_2605.csv

# Score precomputed DrugCLIP representations with the packaged model
python screen.py --molecules MOLECULES.npz --pockets M4_POCKETS.npz

# Recompute the six-channel fusion
python module3_ensemble_dock.py --input data/frozen/six_channel_top5_plus_anchors.csv

# Export four-context evidence
python module4_four_context.py --summary data/frozen/four_context_results.json

# Retrain the GPCR-adapted projections when the original representations are available
python train.py --representations REPRESENTATIONS.pt --projection BASE_PROJECTION.pt --benchmark BENCHMARK.csv
```

Full docking requires AutoDock Vina and licensed Schrodinger Glide. Full membrane MD requires OpenMM and prepared systems. These expensive jobs are not launched by the CPU replay; their frozen outputs and the production scripts used to create them are included.

## Source and logs

`src/production/` contains the retained generation, DrugCLIP training, Glide/Vina, six-channel fusion and four-context analysis scripts. `logs/` contains the model manifest, the actual Module 2 validation report, comparison-model benchmarks and replay audit.

Raw console output, original hardware and wall-clock time were not preserved for every historical experiment. This is stated as a provenance limitation rather than reconstructed after the fact.

Glide is commercial software and is not redistributed. DrugCLIP and public molecular data remain subject to their original licenses and terms.
