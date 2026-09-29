"""Synthetic leakage, baseline and architecture contracts for Experiment B2."""
from pathlib import Path
from contextlib import redirect_stdout
import io
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import experiment_a as A
import experiment_b as B
import experiment_b2 as C


class DynamicTests(unittest.TestCase):
    def setUp(self):
        B.deterministic(17)
        self.regions = {'first': np.array([3, 1, 4]), 'second': np.array([0, 4])}
        self.h = np.random.default_rng(17).normal(size=(4, 5, 6)).astype(np.float32)
        self.r = torch.from_numpy(A.decompose(self.h, self.regions)['local_residual'].astype(np.float32))
        self.data = [dict(row={'path':'replica_02/window_000/atom14/apo_w000.geom2vec.npz'}, r=self.r)]
        self.mu = C.fit_mu(self.data)

    def model(self, mode='attention'):
        return C.DynamicReadout(self.regions, mode, channels=6, residues=5)

    def q(self, model, h=None):
        return C.visible_q(self.h if h is None else h, model.targets.numpy(), model.indices.numpy())

    def test_mu_R2_only_and_R3_rejected(self):
        torch.testing.assert_close(self.mu, self.r.double().mean(0), rtol=0, atol=0)
        contaminated = dict(row={'path':'replica_03/window_000/atom14/apo_w000.geom2vec.npz'}, r=self.r*999)
        with self.assertRaises(ValueError):
            C.fit_mu(self.data+[contaminated])
        with self.assertRaises(ValueError):
            C.fit_mu([])

    def test_D_R2_mean_zero(self):
        model = self.model()
        d = C.dynamic_target(self.r, self.mu, model.targets)
        self.assertLess(d.mean(0).abs().max().item(), 1e-14)

    def test_zero_D_exact_frozen_mean_prediction_and_error(self):
        model = self.model()
        d = C.dynamic_target(self.r, self.mu, model.targets)
        prediction = C.reconstruct_r(torch.zeros_like(d), self.mu, model.targets)
        expected = self.mu[model.targets][None].expand_as(d)
        self.assertTrue(torch.equal(prediction, expected))
        self.assertTrue(torch.equal((prediction-self.r[:,model.targets].double()).square(), d.square()))

    def test_Q_uses_visible_H_directly(self):
        model = self.model()
        q = self.q(model)
        for k, target in enumerate(model.targets.tolist()):
            visible = [j for j in range(5) if j != target]
            center = self.h[:,visible].astype(np.float64).mean(1)
            expected = (self.h[:,model.indices[k]].astype(np.float64)-center[:,None]).astype(np.float32)
            np.testing.assert_array_equal(q[:,k].numpy(), expected)

    def test_masked_H_perturbation_cannot_change_Q_or_prediction(self):
        for mode in B.MODES:
            model = self.model(mode)
            q = self.q(model)
            for k, target in enumerate(model.targets.tolist()):
                changed = self.h.copy()
                changed[:,target] = 1e20
                altered = self.q(model, changed)
                self.assertTrue(torch.equal(q[:,k], altered[:,k]))
                self.assertTrue(torch.equal(model(q)[0][:,k], model(altered)[0][:,k]))

    def test_masked_target_absent_all_readouts_and_padding(self):
        for mode in B.MODES:
            model = self.model(mode)
            self.assertFalse((model.indices==model.targets[:,None]).any())
            invalid = model.indices.numpy().copy()
            invalid[0,0] = model.targets[0]
            with self.assertRaises(ValueError):
                C.visible_q(self.h, model.targets.numpy(), invalid)

    def test_G_not_decoder_or_loss_argument(self):
        model = self.model()
        q = self.q(model)
        d = C.dynamic_target(self.r,self.mu,model.targets)
        with self.assertRaises(TypeError):
            model(q, global_mean=torch.zeros(4,6))
        with self.assertRaises(TypeError):
            C.dynamic_loss(model,q,d,torch.zeros(4,6))
        # A uniform frame shift alters canonical G, while Q and D remain unchanged.
        h = np.arange(120,dtype=np.float32).reshape(4,5,6)
        np.testing.assert_array_equal(self.q(model,h).numpy(),self.q(model,h+1024).numpy())

    def test_uniform_attention_static_mean_equal_and_weights_sum(self):
        for mode in B.MODES:
            model = self.model(mode)
            q = self.q(model)
            for dropout in (0.,.25,.5):
                pooled,weights = model.pool(q,dropout)
                for k in range(len(model.targets)):
                    count = int(model.valid[k].sum())
                    expected = q[:,k,int(count*dropout):count].mean(1)
                    torch.testing.assert_close(pooled[:,k],expected,rtol=1e-6,atol=1e-7)
                torch.testing.assert_close(weights.sum(-1),torch.ones_like(weights[...,0]))

    def test_architecture_and_initialization_exactly_B(self):
        for mode in B.MODES:
            B.deterministic(17)
            old = B.LocalReadout(self.regions,mode,6,5)
            B.deterministic(17)
            new = self.model(mode)
            self.assertEqual(B.state_hash(old),B.state_hash(new))
            self.assertEqual(sum(p.numel() for p in old.parameters()),sum(p.numel() for p in new.parameters()))

    def test_split_and_trajectory_rejections_retained(self):
        rows = [{'path':p} for p in sorted(A.expected_paths())]
        self.assertEqual(len(B.split_rows(rows)[2]),20)
        rows[0]['split']='R3'
        with self.assertRaises(ValueError):
            B.split_rows(rows)
        with self.assertRaises(ValueError):
            B.identity('replica_02/window_000/atom14/apo_w001.geom2vec.npz')

    def test_all_seeds_replay(self):
        for seed in B.SEEDS:
            for mode in B.MODES:
                hashes=[]
                for _ in range(2):
                    B.deterministic(seed)
                    model=self.model(mode)
                    q=self.q(model)
                    d=C.dynamic_target(self.r,self.mu,model.targets)
                    opt=torch.optim.Adam(model.parameters(),lr=.001)
                    for _ in range(2):
                        opt.zero_grad()
                        C.dynamic_loss(model,q,d).backward()
                        opt.step()
                    hashes.append(B.state_hash(model))
                self.assertEqual(*hashes)

    def test_malformed_nonfinite_inputs_rejected(self):
        model=self.model()
        for h in (self.h[0], self.h.astype(int), self.h*np.nan, self.h[:0]):
            with self.assertRaises(ValueError):
                self.q(model,h)
        for q in (self.q(model).double(),self.q(model)*np.nan,self.q(model)[:0]):
            with self.assertRaises(ValueError):
                model(q)

    def test_evaluation_retains_all_channels_regions_and_baselines(self):
        regions={'tiny':np.array([0,3])}
        models={}
        for seed in B.SEEDS:
            for mode in B.MODES:
                B.deterministic(seed)
                models[(seed,mode)]=C.DynamicReadout(regions,mode)
        layout=next(iter(models.values()))
        h=np.random.default_rng(19).normal(size=(2,270,128)).astype(np.float32)
        r=torch.from_numpy(A.decompose(h,regions)['local_residual'].astype(np.float32))
        data=[dict(row={'path':path},h=h,r=r,centers=C.visible_centers(h,layout.targets.numpy()))
              for path in sorted(A.expected_paths()) if path.startswith('replica_02/')]
        mu=C.fit_mu(data)
        with redirect_stdout(io.StringIO()):
            result=C.evaluate(data,models,mu,list(regions))
        self.assertEqual(len(result['summary']),36)
        self.assertEqual(len(result['chronological_blocks']),180)
        self.assertEqual(len(result['trajectory_window_cells']),720)
        for row in result['summary']:
            self.assertEqual(len(row['per_channel_mse']),128)
            self.assertEqual(len(row['per_region_mse']),1)
            self.assertAlmostEqual(row['mse'],np.mean(row['per_channel_mse']),places=13)
            if row['model']=='zero_D':
                self.assertEqual(row['mse'],row['zero_d_mse'])
                self.assertEqual(row['relative_improvement'],0)
                self.assertEqual(row['normalized_dynamic_error'],1)


if __name__ == '__main__':
    unittest.main()
