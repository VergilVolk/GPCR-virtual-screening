# Retained production sources

- `brics_generate.py`, `lstm_generate.py`: candidate generation routes.
- `export_drugclip2023_m4_loto_adapters.py`: training and export entry for the packaged M4-held-out projections.
- `evaluate_drugclip_large_target_loso.py`, `finetune_drugclip_gpcr_screening.py`: model components used by that entry.
- `glide_runner.py`, `vina_runner.py`, `six_channel_fusion.py`: the Glide/Vina six-channel structural stage.
- `four_context_common.py`, `four_context_phase1.py`, `four_context_phase2a.py`, `four_context_phase2b.py`: four-context trajectory analysis.

These are retained project scripts rather than rewritten demonstration algorithms. External software and large prepared assets are still required for a full production rerun.
