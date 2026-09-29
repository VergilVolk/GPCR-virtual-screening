"""Synthetic contracts for the fixed local-only Experiment B."""
import inspect
from pathlib import Path
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import experiment_a as A
import experiment_b as B


class LocalTests(unittest.TestCase):
    def setUp(self):
        B.deterministic(17)
        self.regions = {'ordered': np.array([3, 1, 4]), 'overlap': np.array([0, 4])}
        self.local = torch.randn(4, 5, 6)

    def model(self, mode):
        return B.LocalReadout(self.regions, mode, channels=6, residues=5)

    def test_uniform_equals_mean(self):
        for mode in B.MODES:
            model = self.model(mode)
            pooled, weights = model.pool(self.local)
            expected = torch.stack([self.local[:, model.indices[k, model.valid[k]]].mean(1)
                                    for k in range(len(model.targets))], 1)
            torch.testing.assert_close(pooled, expected, rtol=1e-6, atol=1e-7)

    def test_weights_sum_and_dropout(self):
        model = self.model('attention')
        with torch.no_grad():
            model.query.normal_()
        for rate in (0.0, 0.25, 0.5):
            _, weights = model.pool(self.local, rate)
            torch.testing.assert_close(weights.sum(-1), torch.ones_like(weights[..., 0]))
            self.assertTrue((weights >= 0).all())
            self.assertTrue((weights[:, ~model.valid] == 0).all())

    def test_masked_target_unavailable(self):
        for mode in B.MODES:
            model = self.model(mode)
            for k, target in enumerate(model.targets):
                self.assertNotIn(int(target), model.indices[k].tolist())
                altered = self.local.clone()
                altered[:, target] += 12345
                before, _ = model.pool(self.local)
                after, _ = model.pool(altered)
                torch.testing.assert_close(before[:, k], after[:, k], rtol=0, atol=0)
                leaf = self.local.clone().requires_grad_()
                pooled, _ = model.pool(leaf)
                pooled[:, k].sum().backward()
                self.assertTrue((leaf.grad[:, target] == 0).all())

    def test_global_not_in_objective(self):
        model = self.model('mean')
        with self.assertRaises(TypeError):
            B.local_loss(model, self.local, torch.ones(4, 6))
        with self.assertRaises(TypeError):
            model(self.local, global_mean=torch.ones(4, 6))
        self.assertEqual(list(inspect.signature(B.local_loss).parameters), ['model', 'local'])
        h = np.arange(120, dtype=np.float64).reshape(4, 5, 6)
        a = A.decompose(h, self.regions)
        b = A.decompose(h + 1024, self.regions)
        ra = torch.tensor(a['local_residual'], dtype=torch.float32)
        rb = torch.tensor(b['local_residual'], dtype=torch.float32)
        torch.testing.assert_close(B.local_loss(model, ra), B.local_loss(model, rb), rtol=0, atol=0)

    def test_decomposition_correct(self):
        h = self.local.numpy()
        branches = A.decompose(h, self.regions)
        A.verify(h, branches, self.regions)
        np.testing.assert_allclose(branches['local_residual'] + branches['global_mean'][:, None], h, atol=1e-12)

    def test_all_predefined_seeds_deterministic(self):
        for seed in B.SEEDS:
            for mode in B.MODES:
                states = []
                for _ in range(2):
                    B.deterministic(seed)
                    model = self.model(mode)
                    opt = torch.optim.Adam(model.parameters(), lr=0.001)
                    for _ in range(2):
                        opt.zero_grad()
                        B.local_loss(model, self.local).backward()
                        opt.step()
                    states.append(B.state_hash(model))
                self.assertEqual(*states)

    def test_residue_order_and_membership(self):
        model = self.model('static')
        self.assertEqual(model.targets.tolist(), [3, 1, 4, 0, 4])
        self.assertEqual(model.indices[0].tolist(), [1, 4])
        self.assertEqual(model.indices[1].tolist(), [3, 4])
        self.assertEqual(model.membership[0].tolist(), [1, 2])
        self.assertEqual(model.membership[4, 0].item(), 3)
        self.assertEqual(model.region_ids.tolist(), [0, 0, 0, 1, 1])

    def test_malformed_nonfinite_rejected(self):
        model = self.model('mean')
        for local in (self.local[0], self.local[:, :-1], self.local.double(),
                      self.local * float('nan'), self.local * float('inf')):
            with self.assertRaises(ValueError):
                model.pool(local)
        for ix in ([1], [1, 1], [-1, 2], [1, 5], [1.5, 2.5]):
            with self.assertRaises(ValueError):
                B.LocalReadout({'bad': np.asarray(ix)}, 'mean', 6, 5)

    def test_decoder_capacity_and_initialization_identical(self):
        states = []
        for mode in B.MODES:
            B.deterministic(17)
            model = self.model(mode)
            states.append({k: v for k, v in model.state_dict().items() if k.startswith(('decoder.', 'identity_embedding.'))})
        for key in states[0]:
            for state in states[1:]:
                torch.testing.assert_close(state[key], states[0][key], rtol=0, atol=0)


class SplitTests(unittest.TestCase):
    def setUp(self):
        self.rows = [{'path': p} for p in sorted(A.expected_paths())]

    def test_split_leakage_rejected(self):
        split = B.split_rows(self.rows)
        self.assertEqual(len(split[2]), 20)
        self.assertEqual(len(split[3]), 20)
        self.rows[0]['split'] = 'R3'
        with self.assertRaises(ValueError):
            B.split_rows(self.rows)

    def test_trajectory_window_constraints(self):
        for path in ('replica_02/window_000/atom14/apo_w001.geom2vec.npz',
                     'replica_02/window_005/atom14/apo_w005.geom2vec.npz',
                     'replica_04/window_000/atom14/apo_w000.geom2vec.npz'):
            with self.assertRaises(ValueError):
                B.identity(path)
        with self.assertRaises(ValueError):
            B.split_rows(self.rows[:-1])
        with self.assertRaises(ValueError):
            B.split_rows(self.rows[:-1] + self.rows[:1])


if __name__ == '__main__':
    unittest.main()
