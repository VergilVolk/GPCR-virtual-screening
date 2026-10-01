from pathlib import Path
import sys

import torch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from finetune_drugclip_science2026_context_loto import ContextMetric  # noqa: E402


def test_deployable_triplet_checkpoints_reload_strictly():
    root = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "drugclip_science2026"
        / "litpcba_external_v01"
        / "deployable_triplet"
    )
    paths = sorted(root.glob("pacer_cgm_triplet_seed*.pt"))
    assert len(paths) == 3
    for path in paths:
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        model = ContextMetric(
            checkpoint["dim"], checkpoint["rank"], checkpoint["experts"], checkpoint["seed"]
        )
        model.load_state_dict(checkpoint["state_dict"], strict=True)
        assert sum(value.numel() for value in model.parameters()) == 2178
        assert all(torch.isfinite(value).all() for value in model.state_dict().values())


def test_deployable_projection_triplet_checkpoints_are_finite():
    root = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "drugclip_science2026"
        / "litpcba_external_v01"
        / "deployable_projection_triplet"
    )
    paths = sorted(root.glob("pacer_projection_triplet_seed*.pt"))
    assert len(paths) == 3
    for path in paths:
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        assert checkpoint["n_pairs"] == 8906
        assert len(checkpoint["training_targets"]) == 9
        for branch in ("mol_project", "pocket_project"):
            state = checkpoint[branch]
            assert state
            assert all(torch.is_tensor(value) and torch.isfinite(value).all() for value in state.values())
