import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from run_drugclip_generation_bakeoff import target_percentile_rank


def test_target_percentile_rank_is_target_local():
    table = pd.DataFrame({"target": ["A", "A", "B", "B"]})
    result = target_percentile_rank(table, np.asarray([1.0, 2.0, 100.0, 50.0]))
    np.testing.assert_allclose(result, [0.5, 1.0, 1.0, 0.5])


def test_equal_score_ties_receive_equal_rank():
    table = pd.DataFrame({"target": ["A", "A", "A"]})
    result = target_percentile_rank(table, np.asarray([1.0, 1.0, 2.0]))
    np.testing.assert_allclose(result, [0.5, 0.5, 1.0])
