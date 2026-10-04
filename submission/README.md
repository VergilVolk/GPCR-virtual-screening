# M4 PAM virtual screening submission

This directory is the compact evaluation package for our Track 3 entry. It covers the full decision chain used in the project:

1. candidate chemistry checks;
2. M4 pocket retrieval with a DrugCLIP-derived adapter;
3. multi-conformation docking evidence;
4. four-context molecular-dynamics evidence.

The final output is a ranked list of computational candidates. It is not a claim that any listed molecule is an experimentally confirmed positive allosteric modulator (PAM).

## Quick start

Tested on Windows 11 and Python 3.11. CPU execution is sufficient for the packaged example.

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python predict.py --demo
```

The command writes `results/results.csv` and the four intermediate tables.

Each stage can also be run separately:

```bash
python module1_generate_filter.py --demo
python screen.py --demo
python module3_ensemble_dock.py --demo
python module4_fkg_eval.py --demo
```

## Inputs and outputs

| Stage | Input | Method | Output |
|---|---|---|---|
| Chemistry | candidate ID and SMILES | validity, PAINS, physicochemical and QED checks | `module1_chemistry.csv` |
| Pocket retrieval | frozen 512-dimensional molecule and M4-pocket representations | three M4 leave-one-target-out projection heads; cosine score | `module2_binding.csv` |
| Structure | frozen state-wise docking summary | ten-state Vina ensemble, pocket occupancy and contact gates | `module3_structure.csv` |
| Function | four receptor contexts: apo, ACh, candidate, ACh+candidate | interaction contrast `CP - P - C + A` and cross-replica direction agreement | `module4_function.csv` |
| Integration | outputs above | evidence-preserving merge | `results.csv` |

The example contains 12 candidates. The structure table was computed over ten receptor conformations. The prospective dynamics set contains three candidates, four contexts and three replicas per candidate (36 trajectories, 10 ns each). Raw trajectories are not duplicated in this compact package because of their size; the frozen summaries and full analysis code remain in the repository.

## Training and model files

`models/` contains three M4 leave-one-target-out DrugCLIP adapter checkpoints. They are the exact inference weights used by `screen.py`. The base representations were generated with the documented DrugCLIP checkpoint; redistribution of upstream weights follows the upstream license.

The adapter training entry is:

```bash
python train.py --repo-root .. --epochs 80 --dry-run
```

Remove `--dry-run` only in a full repository checkout and provide the training representations, base projection and contrast table through the corresponding command-line options. Target separation, scaffold controls and the three seeds are recorded in `logs/training_manifest.json`.

## Docking and MD reproduction

The packaged demo replays frozen docking and MD summaries so it can run on a CPU laptop. Full production runs require:

- AutoDock Vina 1.2.7 for state-wise docking;
- the receptor ensemble and prepared structures listed in `data/README.md`;
- OpenMM for membrane MD;
- substantially more time and storage than the compact demo.

Commercial Glide scores used in a separate public benchmark are not required by this submission and are not represented as locally generated results.

## Result fields

`results/results.csv` includes candidate ID, track, SMILES, chemistry checks, M4 binding score, docking energy, ensemble coverage, structural gate, functional-evidence status, model version and claim status. Blank functional fields mean that the candidate has not completed prospective four-context MD.

## Limits

- DrugCLIP and docking assess binding compatibility, not PAM efficacy.
- The four-context endpoint separated retrospective controls, but the prospective candidates did not yet show a robust, cross-replica cooperative signal.
- All final candidates require experimental pharmacology, counter-screening and safety testing.
