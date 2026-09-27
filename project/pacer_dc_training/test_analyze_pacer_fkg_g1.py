"""Standard-library unit tests for frozen PACER-FKG G1 diagnostics."""
import tempfile
import unittest
from pathlib import Path

from analyze_pacer_fkg_g1 import AXES, compute, rankdata, spearman, write_csv


class MathTests(unittest.TestCase):
    def test_rankdata_average_ties(self):
        self.assertEqual(rankdata([1, 2, 2, 4]), [1.0, 2.5, 2.5, 4.0])

    def test_spearman_ties(self):
        self.assertAlmostEqual(spearman([1, 1, 2], [2, 2, 3]), 1.0)

    def test_spearman_constant(self):
        self.assertIsNone(spearman([1, 1, 1], [1, 2, 3]))

    def test_spearman_inversion(self):
        self.assertAlmostEqual(spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1.0)

    def test_csv_escaping(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.csv"
            write_csv(p, [{"axis": "a,b", "rho": 0.5}], ["axis", "rho"])
            self.assertIn('"a,b"', p.read_text(encoding="utf-8-sig"))

    def test_computation_cross_checks_and_counts(self):
        regions = ["functional_region", "stable_core_control"]
        aggregates = {}
        for replica in (2, 3):
            for window in range(5):
                for region in regions:
                    row = {"n_residues": 3}
                    for axis in AXES:
                        v = 1.0 + window * 0.02 if region == "functional_region" else 0.2
                        row[f"{axis}_kernel_u2"] = v
                        row[f"{axis}_kernel_biased2"] = v + 0.01
                        row[f"{axis}_bandwidth"] = 10.0 + window
                        row[f"{axis}_graph_diffused_mean"] = v + 0.1
                        row[f"{axis}_linear"] = v + 0.2
                    aggregates[(replica, window, region)] = row
        summary = []
        for axis in AXES:
            for region in regions:
                val = 0.84 if region == "functional_region" else 0
                summary.append({"axis": axis, "region": region,
                                "replicas": {str(r): {"median_u2_minus_stable": val, "positive_windows": 5 if val else 0,
                                                       "median_graph_diffused_minus_stable": val} for r in (2, 3)},
                                "matched_window_spearman": 1.0 if val else None,
                                "matched_window_graph_diffused_spearman_descriptive": 1.0 if val else None,
                                "preregistered_specificity_gate": bool(val)})
        audit = {"cross_replica_summary": summary}
        rows, reps, comparisons = compute(audit, aggregates)
        self.assertEqual((len(rows), len(reps), len(comparisons)), (60, 12, 6))
        self.assertEqual(sum(x["gate_in_frozen_report"] for x in comparisons), 3)


if __name__ == "__main__":
    unittest.main()
