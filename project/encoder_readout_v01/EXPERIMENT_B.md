# Experiment B v01: local residual readout

This standalone experiment reuses `experiment_a.preflight`, `read_input`,
`load_contract`, and `decompose`. It authenticates all 40 NPZ files before fitting
and reauthenticates the bytes used for decomposition. The authoritative manifest
and region graph are read through A; no source hash table is transcribed.
Experiment A, the frozen encoder, and archived results are unchanged.

Run from the repository root:

```powershell
python -B -m unittest discover -s project/encoder_readout_v01/tests -v
python -B project/encoder_readout_v01/experiment_b.py --input-root C:/projects/pacer_encoder_data/R2R3_full_v01 --output-root project/results/encoder_readout_B_v01/run_001
```

The output must be a new directory. Incomplete runs are retained and cannot be
resumed or overwritten. `frozen_config.json` is saved and printed before training;
its hash ties together provenance, the R3 start marker, and completion record.
The full configuration is the predeclaration. The command has no tuning flags.
All nine R2 fits and their independent deterministic replays finish before the
single R3 evaluation stage begins. R3 failures require investigation, not an
automatic rerun. No checkpoint or model variant is selected using R3.

## Fixed protocol

- Seeds: 17, 43, 101; all reported.
- R2: four condition trajectories, five windows each, all 2,000 frames.
- R3: four separate condition trajectories, five windows each, all 2,000 frames.
- Windows are ordered 0 through 4; within-window local frame IDs remain 0–99.
  No absolute trajectory times are inferred. No frames are shuffled.
- All 109 region memberships are masked individually on every frame. Frozen
  region order and residue order are retained, including overlapping memberships.
- FP64 A decomposition, then unscaled FP32 local features for training. No fitted
  scaling. The simple predictor is the FP64 accumulation of R2 FP32 residuals,
  independently for every residue/channel. Zero is also reported.
- Uniform mean, learned static membership logits, and single-head attention use
  the same target-index embedding and decoder family. Decoder weights start
  identically within each seed. Attention has only a shared 128-to-16 linear
  scorer, tanh, and one 16-component query per region. Values are unprojected R.
- Adam, learning rate 0.001, default betas/epsilon, no weight decay, constant
  schedule, five complete epochs, 25 consecutive frames per batch. No early
  stopping or epoch selection. The final checkpoint is used for all models.
- Decoder: target index embedding of width 8, concatenate with pooled R, linear
  136-to-32, tanh, linear 32-to-128. One shared decoder per comparator/seed.
- CPU FP32 deterministic Torch operations, two threads. Every fit is repeated
  from scratch; state bytes and all epoch losses must agree exactly. Checkpoint
  reload must preserve the state exactly.

G is computed by A and remains a separate, reproducible future branch. The local
model and loss receive only R. G is never concatenated with the decoder input or
used for fitting, stopping, or evaluation. H and its authenticated provenance
retain the information needed to regenerate `[P_local, G]`; neither another full
R archive nor another G archive is written. Only R2/R3 local arrays are cached in
memory, one split at a time. Model checkpoints and the simple predictor are small.

## Metrics and interpretation

MSE averages over frames, membership targets and channels. Normalized error is
MSE divided by target mean square on the same evaluation subset. Relative-simple
error is MSE divided by the R2-residue-mean predictor's MSE on that same subset;
improvement is one minus that ratio. Zero denominators produce JSON null.

Entropy is mean natural-log entropy of the masked-task weights. Effective count
is mean exp(entropy) per task, not exp(mean entropy). Mean and static models also
receive these diagnostics. Predictor baselines have no weights, so their entropy
and effective count are null.

Predeclared token-dropout stress tests remove 0%, 25%, or 50% (rounded down) of
the visible tokens, starting with the earliest in frozen membership order. The
target remains masked. These are deterministic sensitivity tests, not stochastic
training augmentation. Models are never refit. The two predictor baselines are
unchanged by dropout. Five ordered-window aggregate metrics and all 20
trajectory/window cells provide chronological-block sensitivity without
treating windows as independent samples. Full JSON contains every seed/model/
dropout/block, including negative results.

R3 has already been historically inspected. These results are descriptive
engineering validation, not fresh biological confirmation. Temporally correlated
frames, neighboring windows, overlapping regions and shared masked-task targets
are not independent experimental replicates. No uncertainty intervals based on
frame independence are used. A fixed short training budget can underfit all
models; it does not establish their best attainable performance.

The stipulated residual definition uses the full-frame mean. Thus visible R
tokens retain a shared-centering dependence on the original masked H token.
Tests establish exclusion of the target **R token** from its own pooling operation
and exclusion of the separate G branch. They do not establish independence from
masked H or claim to reconstruct H. No decomposition change is introduced to
remove this mathematical dependence.

There are no biological contrasts, kernel calculations, historical success
gates, encoder inference, or automatic commits in this experiment.
