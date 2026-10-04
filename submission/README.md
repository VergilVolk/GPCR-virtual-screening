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

`verify_submission.py` also performs a real Module 2 forward pass for all 2,605 molecules with the three packaged weights and checks the reproduced scores and top-200 order. `predict.py` verifies frozen data, recomputes six-channel fusion, exports four-context records, and writes `results/results.csv`.

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

- `models/reference_not_used/gpcr_retrieval_hard_triplet/`
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
python train.py --representations REPRESENTATIONS.npz --projection BASE_PROJECTION.pt --pairs PAIRS.csv --reference-predictions REFERENCE.csv
```

Full docking requires AutoDock Vina and licensed Schrodinger Glide. Full membrane MD requires OpenMM and prepared systems. These expensive jobs are not launched by the CPU replay; their frozen outputs and the production scripts used to create them are included.

The frozen ViSNet checkpoint used by Module 4 is included under `project/tools/geom2vec_checkpoints/` and verified by SHA-256. A fresh trajectory run still requires the prepared membrane systems, MD trajectories and the matching scientific Python environment maintained by the team.

A traceable four-context MD input example is included in
`data/examples/module4_pacer0073_r1/`: 20 frames per context from the actual
PACER0073 replica-1 production runs, with unchanged corresponding PDB topologies.
See its README for sample verification and regeneration commands. This sample
is also consumed by `run_module4_example.py`, which performs actual frozen
encoding, numerical-state application, four-context contrasts and graph/region
analysis. It is not a full three-replica Stage4 result reproduction. Install the
dedicated environment and run the example as described in
`docs/MODULE4_REPRODUCIBILITY.md`; do not mix its PyTorch 2.6 environment with the
PyTorch 2.8 result-replay environment.

## Source and logs

`src/production/` contains the retained generation, DrugCLIP training, Glide/Vina, six-channel fusion and four-context analysis scripts. `logs/` contains the model manifest, the actual Module 2 validation report, comparison-model benchmarks and replay audit.

Raw console output, original hardware and wall-clock time were not preserved for every historical experiment. This is stated as a provenance limitation rather than reconstructed after the fact.

Glide is commercial software and is not redistributed. DrugCLIP and public molecular data remain subject to their original licenses and terms.
