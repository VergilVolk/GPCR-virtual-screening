from pathlib import Path
import sys

import pytest
import torch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from finetune_drugclip_science2026_context_loto import ContextMetric  # noqa: E402


def _optional_checkpoint_set(root: Path, pattern: str) -> list[Path]:
    """Return a complete optional checkpoint set or skip the external-asset test.

    Large development checkpoints are intentionally excluded from the public
    repository.  A missing directory is therefore not a code failure, while a
    partially copied checkpoint set remains an error.
    """
    paths = sorted(root.glob(pattern))
    if not paths:
        pytest.skip(f"optional external checkpoint bundle is not installed: {root}")
    assert len(paths) == 3, f"incomplete checkpoint bundle in {root}: found {len(paths)}/3"
    return paths


def test_deployable_triplet_checkpoints_reload_strictly():
    root = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "drugclip_science2026"
        / "litpcba_external_v01"
        / "deployable_triplet"
    )
    paths = _optional_checkpoint_set(root, "pacer_cgm_triplet_seed*.pt")
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
    paths = _optional_checkpoint_set(root, "pacer_projection_triplet_seed*.pt")
    for path in paths:
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        assert checkpoint["n_pairs"] == 8906
        assert len(checkpoint["training_targets"]) == 9
        for branch in ("mol_project", "pocket_project"):
            state = checkpoint[branch]
            assert state
            assert all(torch.is_tensor(value) and torch.isfinite(value).all() for value in state.values())
