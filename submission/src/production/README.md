# Retained production sources

- `brics_generate.py`, `lstm_generate.py`: candidate generation routes.
- `finetune_drugclip_gpcr_retrieval.py`: training code for the packaged GPCR-adapted projections.
- `finetune_drugclip_muscarinic_triplet.py`: shared projection definition used by the trainer.
- `glide_runner.py`, `vina_runner.py`, `six_channel_fusion.py`: the Glide/Vina six-channel structural stage.
- `four_context_common.py`, `four_context_phase1.py`, `four_context_phase2a.py`, `four_context_phase2b.py`: four-context trajectory analysis.

These are retained project scripts rather than rewritten demonstration algorithms. External software and large prepared assets are still required for a full production rerun.
