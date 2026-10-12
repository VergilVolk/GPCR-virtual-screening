# Module 4 frozen computation reproduction

The standalone entry point recomputes the actual input example from coordinates:
four MD contexts → unique frozen receptor mapping → strict-load ViSNet → C1-BS256
→ one STATE_MOTION/SIGNED_DRIFT block → frozen normalization and RFF → matched
contrasts → frozen graph diffusion → region vectors and a CSV summary.

All default input, checkpoint, upstream source and asset paths are relative to
the submission package. No historical project checkout, source repository,
external MD backup, historical cache or environment variable is required for
this example. The historical numerical implementations remain unchanged.

## Install the dedicated environment

Use Python 3.11. Module 4 uses PyTorch 2.6.0; the CPU result-replay environment in
`requirements.txt` uses PyTorch 2.8.0. Create separate virtual environments rather
than installing these two requirement files together.

From `submission/`, create a fresh environment and activate it:

```powershell
py -3.11 -m venv .venv-module4
.\.venv-module4\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

For Linux use `python3.11 -m venv .venv-module4` and
`source .venv-module4/bin/activate`. Choose ONE backend below, then install the
remaining scientific packages. The compiled PyG wheels must match PyTorch 2.6.0
and the chosen CPU/CUDA build. No source compilation or approximate neighbor
backend is required by the example.

CPU:

```bash
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
python -m pip install torch-cluster==1.6.3 torch-scatter==2.1.2 torch-sparse==0.6.18 --no-index -f https://data.pyg.org/whl/torch-2.6.0+cpu.html
python -m pip install -r requirements-module4.txt
python -m pip check
```

NVIDIA CUDA 12.6:

```bash
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu126
python -m pip install torch-cluster==1.6.3 torch-scatter==2.1.2 torch-sparse==0.6.18 --no-index -f https://data.pyg.org/whl/torch-2.6.0+cu126.html
python -m pip install -r requirements-module4.txt
python -m pip check
```

The GPU path requires an NVIDIA driver supporting the CUDA 12.6 PyTorch build.
The actual verified machine uses Windows, Python 3.11.16, an RTX 4070 Laptop GPU
with 8 GiB VRAM and driver 610.60. CUDA inference of the complete four-context
example took approximately 22 seconds on that machine; CPU execution on the
same machine took approximately 70 seconds. Runtime is hardware
dependent. Allow several GiB of RAM and at least 200 MiB free disk space for the
sample and computed arrays; CPU execution is slower. These are example resource
guidelines, not hardware minima for full production MD.

## Pinned upstream code and license

`vendor/geom2vec/` contains 54 original upstream files exported directly from
commit `371d642ec1061664f16e49fcac702d07fc8d0b51` of
https://github.com/dinner-group/geom2vec. The entry point loads this bundled
source explicitly and authenticates its bytes against `PINNED_SOURCE_MANIFEST.json`.
There is no dependency on an unpinned PyPI geom2vec release or a local clone.
The complete upstream LICENSE is retained, including notices for MIT and
GPL-derived portions. No upstream model code was edited.

The bundled ViSNet is an upstream pretrained encoder, not a project-trained PAM
classifier. Historical R2 calibration state is loaded, never refitted. The
sample's replica 1 is evaluation-only. The 270-residue ordering, readout widths,
normalization state, bandwidths, RFF weights, graph, regions and contrast formulas
are authenticated before computation.

## Run the complete sample computation

```bash
python run_module4_example.py --device cpu
# Or, with the CUDA environment:
python run_module4_example.py --device cuda
```

The output directory must not already exist. For a second run, choose another
directory with `--output-dir run/module4_example_second`. Input verification
rejects missing or modified files; model loading rejects mismatched weights;
the example requires the official torch-cluster neighbor backend.

Outputs:

- `run/module4_example/computed_arrays.npz`: four [20,270,256] BS256 arrays,
  descriptors, frozen RFF arrays, diffused contrasts and regional 512D vectors.
- `run/module4_example/module4_example_results.csv`: 90 regional summaries,
  covering two branches, five contrasts and nine frozen regions.
- `run/module4_example/run_receipt.json`: input/source/weight/asset identities,
  frozen residue selections, numerical replay checks, runtime and output hashes.

The three historical contrasts retain their exact definitions:
Delta_PAM = CP − P; Delta_AGO = C − A; Delta_INT = CP − P − C + A.
The two Stage4 descriptive contrasts are P − A and CP − C.

This example has only one correlated time block from replica 1. It demonstrates
the core computation, not the full three-replica result, repeat agreement,
independent-block inference, drug potency or PAM functionality. It does not
overwrite the final eight-candidate table or use the full Stage4 completion flag.

## Regression checks and full production entry points

```bash
python -m unittest test_module4_example
```

The sample descriptor adapter is checked for exact equality with the first block
of the formal Stage4 reduction, invalid/incomplete input rejection and the
zero-motion identity. Historical frozen asset authentication runs during every
real example computation.

CPU and CUDA both execute the declared official neighbor method, but their
compiled radius-graph backends can select different neighbors at the cap of 32.
Their regional numbers are not asserted to be interchangeable. Use the CUDA
environment to compare with the original CUDA production pathway; the CPU run
is a separately logged calculation example. Reference CSVs and run receipts
for both tested devices are provided under `results/` and `logs/`. No alternate
radius-graph fallback is installed. Deterministic replay is checked on the same
device using the frozen relative-RMSE criterion of 1e-4.

For full-production use, keep the existing 36-input authentication and three
independent replicas. The new sample entry point does not relax those gates.
Set `PACER_STAGE4_FROZEN_ROOT` to the absolute package root and use explicit
`--source-root` and `--geom2vec-source` arguments. The latter must be an upstream
Git checkout at the pinned commit, because the historical full Phase1 runner
also verifies Git provenance. Run its smoke before its full run, followed by
the existing Phase2a and Phase2b runners. This path still requires the full
production inputs; do not substitute the short sample for those inputs.
