import unittest
import numpy as np

from pacer_mechanism_cv import cosine, extract_frame_features, factorial_contrast


class MechanismCVTests(unittest.TestCase):
    def setUp(self):
        self.mapping = {10: 0, 20: 1, 30: 2, 40: 3}
        self.config = {
            "regions": {"a": [10, 20], "b": [30, 40]},
            "ca_distance_pairs": [[10, 20], [20, 30]],
            "sidechain_distance_pairs": [[10, 20]],
            "region_centroid_pairs": [["a", "b"]],
        }
        arr = np.zeros((8, 4, 14, 3), dtype=float)
        for f in range(8):
            for r in range(4):
                base = np.array([r * 2.0, 0.1 * f, 0.0])
                arr[f, r, 0:5] = base + np.array([[0,0,0], [0.2,0,0], [0.4,0,0], [0.6,0,0], [0.2,0.3,0]])
        self.arr = arr

    def test_rigid_body_invariance(self):
        names, x = extract_frame_features(self.arr, "AAAA", self.mapping, self.config)
        shifted = self.arr.copy()
        shifted[:, :, :5, :] += np.array([10.0, -8.0, 2.0])
        _, y = extract_frame_features(shifted, "AAAA", self.mapping, self.config)
        self.assertEqual(len(names), x.shape[1])
        np.testing.assert_allclose(x, y, atol=1e-10)

    def test_additive_null_and_injected_interaction(self):
        zero = np.zeros(3); a = np.array([1., 0., 0.]); c = np.array([0., 2., 0.])
        null = factorial_contrast({"apo": zero, "probe_only": a,
                                   "candidate_no_probe": c, "candidate_probe": a + c},
                                  {"candidate_probe": 1, "probe_only": -1,
                                   "candidate_no_probe": -1, "apo": 1})
        np.testing.assert_allclose(null, 0.0)
        synergy = np.array([0., 0., 3.])
        hit = factorial_contrast({"apo": zero, "probe_only": a,
                                  "candidate_no_probe": c, "candidate_probe": a + c + synergy},
                                 {"candidate_probe": 1, "probe_only": -1,
                                  "candidate_no_probe": -1, "apo": 1})
        np.testing.assert_allclose(hit, synergy)
        self.assertAlmostEqual(cosine(hit, synergy), 1.0)


if __name__ == "__main__":
    unittest.main()
