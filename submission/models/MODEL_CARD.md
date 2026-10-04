# Model card: M4 pocket retrieval adapter

## Purpose

The model ranks small molecules by compatibility with an M4 allosteric pocket representation. It is an early binding filter and is not a direct predictor of PAM efficacy.

## Architecture

The package contains three projection-head checkpoints. Each head maps frozen 512-dimensional molecule and pocket representations into a normalized 256-dimensional space. Molecule-pocket similarity is the cosine score, averaged over seeds and pocket conformations.

## Training protocol

- backbone representations: DrugCLIP-derived;
- adaptation: muscarinic/GPCR family transfer with M4 held out from fitting;
- seeds: 20260925, 20260926 and 20260927;
- inference model: equal-weight three-seed ensemble;
- validation: target-level and scaffold-aware checks recorded in the full repository.

## Inputs

NPZ files containing `molecule_ids`, `molecule_representations`, `pocket_ids` and `pocket_representations`.

## Checkpoint hashes

| Seed | SHA-256 |
|---|---|
| 20260925 | `a954c7c16fe857b5e3b02f1fbfdf7b2b1284d3559d86448b33736e16dd5d7398` |
| 20260926 | `3cc92da9f066a36f9d65b47b7859e69ee7a456d48047a3e1bbd509283d50d4b9` |
| 20260927 | `578ab8349fa8723f5609067e6511ebdfab8541dab8558e883bb7affeb7a8780c` |

## Output

A cosine binding score and rank for each molecule.

## Known limits

The score estimates pocket compatibility. It does not establish affinity, cooperativity, intrinsic agonism, selectivity, safety or clinical activity. Scores should be combined with structural and functional evidence and tested experimentally.
