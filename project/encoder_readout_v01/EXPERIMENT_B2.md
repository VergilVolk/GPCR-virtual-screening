# Experiment B2 v01

B2 changes the target to D = R - mu and fixes the masked-input centering leak.
It reuses A's authenticated loader/decomposition and B's split checks, model
construction, initialization, seeds and five-epoch fitting protocol. A and B
source files and archived outputs are unchanged.

```powershell
python -B -m unittest discover -s project/encoder_readout_v01/tests -v
python -B project/encoder_readout_v01/experiment_b2.py --input-root C:/projects/pacer_encoder_data/R2R3_full_v01 --output-root project/results/encoder_readout_B2_v01/run_001
python -B project/encoder_readout_v01/report_b2.py project/results/encoder_readout_B2_v01/run_001
```

The runner verifies B's archived source/output hashes and authenticates all 40
source NPZ files through A. It runs and captures the entire regression suite,
then prints/saves a frozen configuration before real fitting. A fresh result
directory is mandatory. An incomplete run cannot be overwritten or resumed.
After fitting, all nine models must pass full independent deterministic replay
and checkpoint reload checks before the one locked R3 evaluation stage.

## Numerical contracts

R is computed using A's FP64 full-frame decomposition, then rounded to FP32 at
the same boundary as B. mu is the FP64 accumulation of those R2-only values,
using B's exact window/frame order. It must equal B's saved predictor bitwise.
D is formed in FP64 for evaluation and cast FP32 for the unchanged FP32 training
loss. There is no scaling or R3 re-centering. Zero D plus mu exactly reproduces
the frozen predictor, including its baseline squared errors. R2 mean D must be
below 1e-12 in absolute value.

For every target i, the preprocessing explicitly gathers H[k != i], computes
its FP64 mean, and subtracts that mean from the visible regional H tokens.
It never subtracts H_i from a full-frame sum or adds a correction based on R_i.
The target cannot enter the centering reduction, regional values or padding.
Q is cast FP32 for the model. Full-frame canonical G is never a decoder or loss
input. A preserves the decomposition contract; G remains reproducible separately
from the immutable source archive.

The predefined dropout remains a readout-token stress test: remove the first
floor(fraction * visible region count) tokens in frozen membership order after
Q construction. Fractions remain 0, 0.25, 0.5. G_minus_i always uses all 269
non-target H residues. Thus dropout does not claim to remove those residues'
contributions to visible-frame centering. This follows the specified k != i
formula while retaining B's original readout dropout operation.

## Architecture and fitting

DynamicReadout inherits B's constructor and forward decoder; only pooling's
input interface changes to already-masked Q. Tests compare all initial state
bytes and parameter counts against B. Mean/static/attention retain 10,768 /
10,877 / 12,976 trainable parameters. Attention retains a shared 128-to-16 scorer,
tanh, 9 region queries of width 16, one head and no value projection.

All R2 frames and 109 target memberships are used. Replica x condition is the
trajectory unit, with all five ordered windows kept together. All R3 trajectories
are evaluation only. Seeds 17, 43, 101; batch size 25; Adam 0.001 with the same
betas/epsilon and zero weight decay; constant learning rate; exactly five epochs;
no shuffling, early stopping, architecture search or seed selection. Decreasing
epoch-five loss is reported as an optimization limitation, never extended.

The in-memory split cache contains H, R and per-target visible-only means. Q and
D are constructed on demand in chronological batches. No new R/Q/D archive is
written. Only small checkpoints, mu, JSON metrics, CSV tables, provenance and
logs are saved.

## Diagnostics and interpretation

Primary comparison is D MSE versus zero-D MSE = mean(D squared). Absolute
improvement is baseline minus model MSE; excess is the opposite sign. Relative
improvement is 1 - model MSE / baseline MSE. Normalized dynamic error is this
same ratio without subtraction from one; it is not channelwise whitening.
R3 energy remains relative to R2 mu, including any distribution shift.

Save every seed, comparator and dropout level for both splits; all five ordered
window blocks and all trajectory/window cells; all 128 per-channel MSE/baseline
energies; all nine per-region MSE, natural-log entropy and mean exp(entropy).
Zero denominators produce null normalized metrics. No channel or region result
selects a model. Full channels are exported for R3, including dropout cases.

Attention beating mean is insufficient if it cannot beat zero-D. Benefits that
reverse with seed/dropout/block are unstable. If all learned models lose to
zero-D, this objective has not shown incremental dynamic information and no
additional model complexity is justified by the result. Tests do not prove that
cached encoder features lack pre-existing contextual information; they prove
that changing masked H_i alone cannot change the Q values supplied for that task.
Q also uses a whole-visible-frame reference, not exclusively regional H.

R3 is historically inspected and reused; evaluation is descriptive engineering
validation. Correlated frames, windows and overlapping memberships are not
independent samples. No biological efficacy claim, biological contrast, kernel
calculation, new encoder inference, long MD or automatic commit is performed.
