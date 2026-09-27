import importlib.util
import unittest
from pathlib import Path

import numpy as np


MODULE = Path(__file__).parents[1] / "pacer_dc_training" / "pacer_factorial_kernel_graph.py"
SPEC = importlib.util.spec_from_file_location("pacer_fkg", MODULE)
FKG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FKG)


def contexts(seed=7, interaction=0.0, agonism=0.0, n=80, d=12):
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(n, d))
    probe = rng.normal(size=(n, d)) + np.r_[np.ones(3), np.zeros(d - 3)]
    candidate = rng.normal(size=(n, d)) + agonism * np.r_[np.zeros(3), np.ones(3), np.zeros(d - 6)]
    both = rng.normal(size=(n, d))
    both += np.r_[np.ones(3), np.zeros(d - 3)]
    both += agonism * np.r_[np.zeros(3), np.ones(3), np.zeros(d - 6)]
    both += interaction * np.r_[np.zeros(6), np.ones(3), np.zeros(d - 9)]
    return {"apo": base, "probe_only": probe, "candidate_no_probe": candidate, "candidate_probe": both}


class FactorialKernelTests(unittest.TestCase):
    def test_interaction_axis_detects_nonadditivity(self):
        signs = {"candidate_probe": 1, "probe_only": -1, "candidate_no_probe": -1, "apo": 1}
        null = FKG.signed_kernel_stat(contexts(interaction=0.0), signs)["unbiased_squared"]
        signal = FKG.signed_kernel_stat(contexts(interaction=2.0), signs)["unbiased_squared"]
        self.assertGreater(signal, null + 0.05)

    def test_agonism_axis_is_separate_from_interaction(self):
        groups = contexts(interaction=0.0, agonism=2.0)
        ago = FKG.signed_kernel_stat(groups, {"candidate_no_probe": 1, "apo": -1})["unbiased_squared"]
        interaction = FKG.signed_kernel_stat(
            groups, {"candidate_probe": 1, "probe_only": -1, "candidate_no_probe": -1, "apo": 1}
        )["unbiased_squared"]
        self.assertGreater(ago, interaction + 0.05)

    def test_frame_order_does_not_change_distribution_statistic(self):
        signs = {"candidate_probe": 1, "probe_only": -1, "candidate_no_probe": -1, "apo": 1}
        groups = contexts(interaction=1.5)
        original = FKG.signed_kernel_stat(groups, signs)["unbiased_squared"]
        reversed_groups = {name: values[::-1] for name, values in groups.items()}
        reversed_value = FKG.signed_kernel_stat(reversed_groups, signs)["unbiased_squared"]
        self.assertTrue(np.isclose(original, reversed_value, atol=1e-12))


if __name__ == "__main__":
    unittest.main()
