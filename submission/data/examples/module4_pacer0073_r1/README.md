# Module 4 four-context MD input example

This sample comes from the actual PACER0073 Stage4 production simulations:
cluster 0, replica 1, paired seed 27101. It includes all four contexts:

- `cluster0__apo`: receptor alone.
- `cluster0__probe_only`: receptor and ACh.
- `PACER0073__candidate_no_probe`: receptor and candidate.
- `PACER0073__candidate_probe`: receptor, candidate and ACh.

Each folder contains an unchanged copy of its original `minimized.pdb` and
`trajectory_sample.dcd`. All atoms and their order, coordinates and periodic
boxes are retained. Coordinates are in angstrom and timestamps in picoseconds.
No alignment, interpolation, receptor cropping, simulation or fitting occurred.
Protein selection remains the existing frozen 270-residue mapping contract.

The excerpts use original zero-based frames `0, 5, 10, ..., 95`: the first 20
frames after the Stage4 frozen `frames[::5]` stride. Their timestamps are
approximately 10 to 960 ps, spaced by 50 ps. This is one nominal 1 ns analysis
block with 19 internal lags spanning 950 ps. It is an input example, not the
complete 10 ns trajectory or the three-replica evaluation. It does not establish
PAM activity or reproduce the final biological interpretation by itself.

`manifest.json` records the original relative paths, full original-file hashes,
sample hashes, atom counts, frame indices, timestamps and extraction checks.
The source-root field records historical provenance; verification of the sample
does not require that root to exist on the reviewing computer.

## Verify the packaged example

The utility uses NumPy and MDAnalysis 2.10.0. From the `submission` directory:

```bash
python -m pip install MDAnalysis==2.10.0
python src/examples/extract_four_context_sample.py --verify-only --output-dir data/examples/module4_pacer0073_r1
```

Verification checks all eight file hashes, atom/frame counts, finite coordinates
and timestamps. The expected result is `PASS`, four contexts, 20 frames each.
This command verifies MD inputs. To execute ViSNet and the frozen FKG computation,
install the environment in `docs/MODULE4_REPRODUCIBILITY.md`, then run
`python run_module4_example.py --device cpu` (or `--device cuda`).

## Regenerate from the original Stage4 data

```bash
python src/examples/extract_four_context_sample.py --source-root /path/to/PACER_STAGE4_MD_backup --output-dir /path/to/new/sample
```

The destination must not exist. Source files are opened read-only, hashed before
and after extraction, and selected coordinates and periodic boxes are compared
with the derived trajectories. The production data layout is unchanged.

The existing Stage4 `--smoke` command authenticates all 36 full production inputs.
This excerpt must not be renamed or substituted for those inputs. A dedicated
core-analysis example entry point is provided as `run_module4_example.py` and
does not bypass the full-production authentication gates.
