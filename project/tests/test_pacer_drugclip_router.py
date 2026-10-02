import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from pacer_drugclip_router import RouterConfig, route_scores


def example_frame():
    return pd.DataFrame({
        "pair_id": ["m4_a", "m4_b", "b2_a", "b2_b"],
        "target": ["M4R", "M4R", "B2AR", "B2AR"],
        "score_2023_gpcr_loto": [0.9, 0.1, 0.9, 0.1],
        "score_2026_13t_famaug": [0.1, 0.9, 0.1, 0.9],
    })


def test_m4_uses_only_2023_rank():
    out = route_scores(example_frame()).set_index("pair_id")
    assert out.loc["m4_a", "pacer_binding_score"] == 1.0
    assert out.loc["m4_b", "pacer_binding_score"] == 0.5
    assert out.loc["m4_a", "route"] == "m4_2023_gpcr_loto"


def test_non_m4_uses_fixed_rank_fusion():
    out = route_scores(example_frame()).set_index("pair_id")
    np.testing.assert_allclose(out.loc[["b2_a", "b2_b"], "pacer_binding_score"], [0.75, 0.75])
    assert set(out.loc[["b2_a", "b2_b"], "route"]) == {"gpcr_fixed_rank_fusion"}


def test_claim_boundary_is_always_non_functional():
    out = route_scores(example_frame())
    assert not out.eligible_for_pam_claim.any()


def test_duplicate_pair_id_is_rejected():
    frame = example_frame()
    frame.loc[1, "pair_id"] = frame.loc[0, "pair_id"]
    with pytest.raises(ValueError, match="pair_id"):
        route_scores(frame)


def test_non_frozen_weights_are_marked():
    out = route_scores(example_frame(), RouterConfig(old_weight=0.7, new_weight=0.3))
    assert not out.frozen_protocol.any()


def test_m4_aliases_share_one_ranking_pool():
    frame = example_frame().iloc[:2].copy()
    frame.loc[0, "target"] = "CHRM4"
    frame.loc[1, "target"] = "M4"
    out = route_scores(frame).set_index("pair_id")
    assert out.loc["m4_a", "final_rank"] == 1
    assert out.loc["m4_b", "final_rank"] == 2
