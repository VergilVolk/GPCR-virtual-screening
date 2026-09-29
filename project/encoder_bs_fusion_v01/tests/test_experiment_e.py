from __future__ import annotations

import unittest

import torch

from project.encoder_bs_fusion_v01.experiment_e import FusionDecoder


class FusionTests(unittest.TestCase):
    def setUp(self):
        self.m = torch.randn(5, 128)
        self.b = torch.randn(5, 128)
        self.s = torch.randn(5, 128)

    def test_avg_is_exact(self):
        model = FusionDecoder("AVG", 17)
        z, gate = model.fuse(self.m, self.b, self.s)
        torch.testing.assert_close(z, (self.b + self.s) / 2, rtol=0, atol=0)
        self.assertIsNone(gate)

    def test_static_gate_is_channelwise_and_initially_half(self):
        model = FusionDecoder("STATIC", 17)
        z, gate = model.fuse(self.m, self.b, self.s)
        self.assertEqual(tuple(model.gate_logits.shape), (128,))
        torch.testing.assert_close(gate, torch.full_like(gate, 0.5), rtol=0, atol=0)
        torch.testing.assert_close(z, (self.b + self.s) / 2, rtol=0, atol=0)

    def test_dynamic_architecture_and_shape(self):
        model = FusionDecoder("DYNAMIC", 17)
        self.assertEqual(model.gate_in.in_features, 256)
        self.assertEqual(model.gate_in.out_features, 16)
        self.assertEqual(model.gate_out.in_features, 16)
        self.assertEqual(model.gate_out.out_features, 128)
        z, gate = model.fuse(self.m, self.b, self.s)
        self.assertEqual(tuple(z.shape), (5, 128))
        self.assertTrue(bool(((gate > 0) & (gate < 1)).all()))

    def test_decoder_initialization_matches_across_kinds(self):
        models = [FusionDecoder(kind, 43) for kind in ("M", "AVG", "STATIC", "DYNAMIC")]
        for model in models[1:]:
            torch.testing.assert_close(model.decoder.weight, models[0].decoder.weight, rtol=0, atol=0)
            torch.testing.assert_close(model.decoder.bias, models[0].decoder.bias, rtol=0, atol=0)


if __name__ == "__main__":
    unittest.main()

