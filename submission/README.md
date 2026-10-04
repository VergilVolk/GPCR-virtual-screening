# M4 PAM virtual screening code

This directory contains the four modules used in our Track 3 workflow. Screening results,
candidate rankings, docking poses and MD trajectories are not included.

## Installation and check

Tested on Windows 11 with Python 3.11. The code-level smoke test runs on CPU.

```bash
pip install -r requirements.txt
python verify_submission.py
```

The check verifies the three packaged DrugCLIP adapter weights and runs all four modules on
temporary synthetic inputs. It does not reproduce or disclose the submitted candidate list.

## Modules

1. `module1_generate_filter.py`: applies the physicochemical and simplified pharmacophore
   rules used after BRICS and SMILES-LSTM generation.
2. `screen.py`: applies three M4 leave-one-target-out DrugCLIP projection heads and averages
   their cosine scores.
3. `module3_ensemble_dock.py`: combines Glide-PDB, Glide-BEmin, Glide-BEavg, Vina-PDB,
   Vina-BEmin and Vina-BEavg after conversion to protocol-specific rank percentiles.
4. `module4_fkg_eval.py`: exports the frozen four-context result based on apo, ACh-only,
   candidate-only and candidate+ACh simulations.

The original production scripts are kept in `src/production/`. They cover BRICS generation,
SMILES-LSTM generation, Glide/Vina docking, six-channel fusion and the full four-context
analysis. The top-level scripts are small inference interfaces around the same calculations.

## Unified inference entry

```bash
python predict.py \
  --library INPUT_LIBRARY.csv \
  --molecule-representations MOLECULES.npz \
  --pocket-representation M4_POCKET.npz \
  --docking-six-channel SIX_CHANNEL.csv \
  --stage4-summary STAGE4_SUMMARY.json \
  --stage4-id-map OPTIONAL_ID_MAP.csv \
  --output run/results.csv
```

Required identifiers must be consistent across the library, representation and docking files.
If Stage 4 uses historical display IDs, provide a two-column map named `candidate_id` and
`stage4_candidate_id`. Omit both Stage 4 options when four-context MD has not been run.

Build a local review page after inference:

```bash
python dashboard.py --results run/results.csv --output run/dashboard.html
```

## Expected input fields

- library CSV: `candidate_id` (or `molecule_id`) and `smiles`;
- molecule NPZ: `molecule_ids`, `molecule_representations`;
- pocket NPZ: `pocket_representations`;
- docking CSV: `candidate_id` plus raw score and `rankpct` columns for all six channels;
- Stage 4 JSON: output of the frozen prospective evaluation code.

## Training

`train.py` calls the repository's DrugCLIP triplet fine-tuning script. The three final adapter
weights used by `screen.py` are in `models/`; their SHA-256 values are listed in the model card
and checked by `verify_submission.py`.

## Full-run requirements

The compact smoke test does not invoke docking or MD. A full run needs AutoDock Vina,
licensed Schrodinger Glide, OpenMM, prepared receptor structures and receptor trajectories.
Those installations and large or licensed inputs are not bundled.

## Scope

DrugCLIP and docking measure binding compatibility. Four-context dynamics tests conditional
structural responses. These calculations do not establish PAM efficacy without functional
experiments.
