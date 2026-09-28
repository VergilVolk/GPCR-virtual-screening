"""B2: fixed dynamic-target reconstruction with visible-only H centering."""
from __future__ import annotations

import argparse
import copy
import io
import json
from pathlib import Path
import sys
import time
import unittest

import numpy as np
import torch

import experiment_a as A
import experiment_b as B

VERSION = 'encoder_readout_B2_v01'
BASE = A.WORKSPACE / 'project/results' / VERSION
PREVIOUS = A.WORKSPACE / 'project/results/encoder_readout_B_v01/run_001'


class DynamicReadout(B.LocalReadout):
    """Identical trainable architecture; input is already masked, centered Q."""
    def pool(self, q, dropout=0.0):
        expected = (len(self.targets), self.indices.shape[1], self.channels)
        if q.ndim != 4 or len(q) == 0 or q.shape[1:] != expected or q.dtype != torch.float32 or not torch.isfinite(q).all():
            raise ValueError('Expected finite FP32 visible Q (frame,task,token,channel)')
        if dropout not in (0.0, 0.25, 0.5):
            raise ValueError('Undeclared dropout')
        valid = self.valid.clone()
        for k in range(len(valid)):
            count = int(valid[k].sum())
            valid[k, :int(count * dropout)] = False
        if self.mode == 'attention':
            scores = (torch.tanh(self.scorer(q)) * self.query[self.region_ids][None, :, None, :]).sum(-1)
        elif self.mode == 'static':
            scores = self.logits[self.membership][None].expand(len(q), -1, -1)
        else:
            scores = torch.zeros(q.shape[:-1], dtype=q.dtype)
        weights = scores.masked_fill(~valid[None], -torch.inf).softmax(-1)
        return (weights[..., None] * q).sum(-2), weights


def visible_centers(h, targets):
    """Gather only H[k != target]. Never subtract target from a full-frame sum."""
    h = np.asarray(h)
    targets = np.asarray(targets)
    if h.ndim != 3 or min(h.shape) < 1 or h.shape[1] < 2 or h.dtype.kind != 'f' or not np.isfinite(h).all():
        raise ValueError('Invalid H')
    if targets.ndim != 1 or not len(targets) or targets.dtype.kind not in 'iu' or np.any(targets < 0) or np.any(targets >= h.shape[1]):
        raise ValueError('Invalid target indices')
    x = np.array(h, dtype=np.float64, order='C', copy=True)
    centers = {}
    for target in dict.fromkeys(targets.tolist()):
        visible = np.arange(h.shape[1]) != target
        centers[target] = x[:, visible, :].mean(axis=1, dtype=np.float64)
    return np.stack([centers[int(target)] for target in targets], axis=1)


def q_from_centers(h, centers, targets, indices):
    targets, indices = np.asarray(targets), np.asarray(indices)
    if indices.ndim != 2 or indices.shape[0] != len(targets) or indices.dtype.kind not in 'iu' or np.any(indices < 0) or np.any(indices >= h.shape[1]):
        raise ValueError('Invalid visible region indices')
    if np.any(indices == targets[:, None]):
        raise ValueError('Masked target present in Q inputs')
    if centers.shape != (len(h), len(targets), h.shape[2]):
        raise ValueError('Invalid leave-one-out centers')
    q = np.asarray(h[:, indices, :], dtype=np.float64) - centers[:, :, None, :]
    if not np.isfinite(q).all():
        raise ValueError('Nonfinite Q')
    return torch.from_numpy(q.astype(np.float32))


def visible_q(h, targets, indices):
    return q_from_centers(h, visible_centers(h, targets), targets, indices)


def fit_mu(data):
    """Only R2 frames can contribute; R3 input is rejected, not ignored."""
    if not data or any(B.identity(item['row']['path'])[0] != 2 for item in data):
        raise ValueError('Residue means require R2-only data')
    return sum(item['r'].double().sum(0) for item in data) / sum(len(item['r']) for item in data)


def dynamic_target(r, mu, targets):
    return r[:, targets].double() - mu[targets][None]


def reconstruct_r(predicted_d, mu, targets):
    return predicted_d.double() + mu[targets][None]


def dynamic_loss(model, q, target_d):
    predicted, _ = model(q)
    if target_d.shape != predicted.shape or not torch.isfinite(target_d).all():
        raise ValueError('Invalid D target')
    return (predicted - target_d.float()).square().mean()


def load_windows(root, rows, sequence, regions, layout):
    data = []
    targets = layout.targets.numpy()
    for row in rows:
        arrays = A.read_input(root, row, sequence)
        h = np.array(arrays['residue_features'], dtype=np.float32, order='C', copy=True)
        r = torch.from_numpy(A.decompose(h, regions)['local_residual'].astype(np.float32))
        data.append(dict(row=row, h=h, r=r, centers=visible_centers(h, targets)))
    return data


def batches(data, layout):
    targets, indices = layout.targets.numpy(), layout.indices.numpy()
    for item in data:
        for start in range(0, len(item['h']), 25):
            sl = slice(start, start + 25)
            yield q_from_centers(item['h'][sl], item['centers'][sl], targets, indices), item['r'][sl]


def fit(data, mu, regions, mode, seed, cfg):
    B.deterministic(seed)
    model = DynamicReadout(regions, mode)
    settings = cfg['optimizer']
    opt = torch.optim.Adam(model.parameters(), lr=settings['lr'], betas=tuple(settings['betas']),
                           eps=settings['eps'], weight_decay=settings['weight_decay'])
    history = []
    for epoch in range(cfg['epochs']):
        total = 0.0
        for q, r in batches(data, model):
            opt.zero_grad(set_to_none=True)
            loss = dynamic_loss(model, q, dynamic_target(r, mu, model.targets))
            if not torch.isfinite(loss):
                raise ValueError('Nonfinite training loss')
            loss.backward()
            opt.step()
            total += float(loss.detach()) * len(q)
        history.append(total / sum(len(item['r']) for item in data))
    return model.eval(), history


def derived_metrics(mse, zero):
    return dict(mse=mse, zero_d_mse=zero, absolute_improvement=zero-mse,
                excess_mse=mse-zero, normalized_dynamic_error=mse/zero if zero else None,
                relative_improvement=1-mse/zero if zero else None)


@torch.no_grad()
def evaluate(data, models, mu, region_names):
    """One pass through each split, with all locked models/diagnostics together."""
    layout = next(iter(models.values()))
    cells = []
    region_ids = layout.region_ids.numpy()
    for item in data:
        target = dynamic_target(item['r'], mu, layout.targets)
        baseline_channel = target.square().mean((0, 1)).numpy()
        baseline_region = np.array([target[:, region_ids == j].square().mean().item() for j in range(len(region_names))])
        # Validate prediction and MSE equivalence without another model evaluation.
        zero_pred = reconstruct_r(torch.zeros_like(target), mu, layout.targets)
        if not torch.equal(zero_pred, mu[layout.targets][None].expand_as(target)):
            raise ValueError('Zero-D does not reproduce residue-mean predictor')
        if not torch.equal((zero_pred-item['r'][:, layout.targets].double()).square(), target.square()):
            raise ValueError('Zero-D baseline error differs from frozen mean baseline')
        keys = [(seed, mode, drop) for seed in B.SEEDS for mode in (*B.MODES, 'zero_D') for drop in (0.0, 0.25, 0.5)]
        totals = {key: dict(channel=np.zeros(128), entropy=0., effective=0.,
                           region=np.zeros(len(region_names)), region_entropy=np.zeros(len(region_names)),
                           region_effective=np.zeros(len(region_names))) for key in keys}
        for q, r in batches([item], layout):
            d = dynamic_target(r, mu, layout.targets)
            for seed, mode, drop in keys:
                if mode == 'zero_D':
                    continue
                pred, weights = models[(seed, mode)](q, drop)
                squared = (pred.double()-d).square()
                entropy = -(weights.double()*weights.double().clamp_min(1e-300).log()).sum(-1)
                entry = totals[(seed, mode, drop)]
                entry['channel'] += squared.mean((0, 1)).numpy()*len(q)
                entry['entropy'] += entropy.mean().item()*len(q)
                entry['effective'] += entropy.exp().mean().item()*len(q)
                for j in range(len(region_names)):
                    keep = region_ids == j
                    entry['region'][j] += squared[:, keep].mean().item()*len(q)
                    entry['region_entropy'][j] += entropy[:, keep].mean().item()*len(q)
                    entry['region_effective'][j] += entropy[:, keep].exp().mean().item()*len(q)
        replica, condition, window = B.identity(item['row']['path'])
        for seed, mode, drop in keys:
            entry = totals[(seed, mode, drop)]
            count = len(item['r'])
            if mode == 'zero_D':
                channel, regional = baseline_channel, baseline_region
                entropy = effective = region_entropy = region_effective = None
            else:
                channel, regional = entry['channel']/count, entry['region']/count
                entropy, effective = entry['entropy']/count, entry['effective']/count
                region_entropy = (entry['region_entropy']/count).tolist()
                region_effective = (entry['region_effective']/count).tolist()
            cells.append(dict(replica=replica, condition=condition, window=window, seed=seed, model=mode, dropout=drop,
                              **derived_metrics(float(channel.mean()), float(baseline_channel.mean())),
                              entropy=entropy, effective_count=effective, per_channel_mse=channel.tolist(),
                              per_channel_zero_d_mse=baseline_channel.tolist(), per_region_mse=regional.tolist(),
                              per_region_zero_d_mse=baseline_region.tolist(), per_region_entropy=region_entropy,
                              per_region_effective_count=region_effective))
        print(f"Evaluated R{replica} {condition} window {window}: all fixed models and dropout levels", flush=True)
    def aggregate(items):
        mse, zero = [float(np.mean([r[k] for r in items])) for k in ('mse', 'zero_d_mse')]
        out = derived_metrics(mse, zero)
        for key in ('entropy', 'effective_count', 'per_channel_mse', 'per_channel_zero_d_mse', 'per_region_mse',
                    'per_region_zero_d_mse', 'per_region_entropy', 'per_region_effective_count'):
            out[key] = np.mean([r[key] for r in items], axis=0).tolist() if items[0][key] is not None else None
        return out
    summary, blocks = [], []
    for seed in B.SEEDS:
        for mode in (*B.MODES, 'zero_D'):
            for drop in (0.0, 0.25, 0.5):
                key = dict(seed=seed, model=mode, dropout=drop)
                selected = [r for r in cells if all(r[k] == v for k, v in key.items())]
                summary.append({**key, **aggregate(selected)})
                for window in range(5):
                    blocks.append({**key, 'window': window, **aggregate([r for r in selected if r['window'] == window])})
    return dict(region_names=region_names, summary=summary, chronological_blocks=blocks, trajectory_window_cells=cells)


def configuration(regions, records):
    old = json.loads((PREVIOUS/'frozen_config.json').read_text(encoding='utf-8'))
    inherited = B.configuration(regions, records)
    unchanged = ('seeds', 'architecture', 'decoder', 'parameter_counts', 'split', 'trajectory_unit', 'masks',
                 'region_order', 'sampling', 'optimizer', 'schedule', 'epochs', 'stopping', 'runtime', 'replay', 'evaluation')
    if any(inherited[k] != old[k] for k in unchanged):
        raise ValueError('Experiment B fixed protocol mismatch')
    cfg = {k: copy.deepcopy(old[k]) for k in unchanged}
    cfg.update(version=VERSION, parent_config_sha256=A.file_sha256(PREVIOUS/'frozen_config.json'),
               comparators=[*B.MODES, 'zero_D'],
               baseline_trainable_parameters=0, frozen_mu_values=270*128,
               target='D = R - mu; canonical A FP64 decomposition then same FP32 R boundary as B; mu accumulated FP64 over R2 only; D FP64 evaluation/FP32 loss',
               preprocessing='No scaling. Per-target G_minus_i = FP64 mean of explicitly gathered H[k != i]; Q=H_visible-G_minus_i, cast FP32. No full-frame sum minus target shortcut.',
               inputs='Only region Q values with target removed reach readouts/decoder. No canonical G or target R/H enters Q construction.',
               dropout_centering='Same B readout-token dropout after Q construction. G_minus_i always excludes target only and uses all other 269 H residues; dropout does not recompute centering.',
               global_branch='Canonical G remains separately reproducible from original archive; never passed to decoder/loss.',
               means='Recompute from R2 only; require exact equality to B frozen R2_residue_mean.pt; no R3 fit or re-centering.',
               loss='Equal weight mean squared D error over frames,109 memberships,128 channels; same FP32 optimization as B',
               baseline='D_hat=0; R_hat=mu+D_hat exactly equals frozen residue-mean predictor',
               metrics={**copy.deepcopy(old['metrics']), 'mse':'mean squared D error',
                        'normalized_error':'MSE / same-subset mean(D^2); equals ratio to zero-D baseline',
                        'relative_simple':'1-MSE/zero_D_MSE is relative improvement; zero_D_MSE-MSE absolute improvement; reverse sign is excess',
                        'per_channel':'All 128 channel MSE and zero-D energies on both splits, every seed/dropout; no selection',
                        'per_region':'All nine region MSE, entropy and effective counts, every seed/dropout; no selection'},
               limitations=['R3 reused and historically inspected: descriptive engineering evaluation only',
                            'Temporal correlation, overlapping regions and target memberships are not independent samples',
                            'Epoch-five decreasing loss is an optimization limitation, never a reason to extend training',
                            'Q removes explicit masked-H centering dependence; precomputed encoder features may themselves contain contextual information',
                            'Q uses whole visible-frame centering, so this is not a strictly region-only preprocessing operation',
                            'No biological efficacy inference; no source encoder updates'])
    return cfg


def archived_snapshot():
    """Hashes only: never inspect historical result values."""
    roots = [A.OUTPUT_BASE, PREVIOUS]
    return {str(p): A.file_sha256(p) for root in roots if root.exists() for p in sorted(root.rglob('*')) if p.is_file()}


def verify_previous():
    audit = json.loads((PREVIOUS/'delivery_audit.json').read_text(encoding='utf-8'))
    for path, expected in audit['delivered_source_hashes'].items():
        if A.file_sha256(A.WORKSPACE/path) != expected:
            raise ValueError('Archived B source hash mismatch: '+path)
    for name, expected in audit['result_hashes'].items():
        if A.file_sha256(PREVIOUS/name) != expected:
            raise ValueError('Archived B result hash mismatch: '+name)


def run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if not BASE.resolve().is_relative_to(A.WORKSPACE) or output == BASE.resolve() or not output.is_relative_to(BASE.resolve()):
        raise ValueError('B2 output must be a new child of versioned result base')
    if output.exists():
        raise FileExistsError(output)
    if output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError('Input/output overlap')
    verify_previous()
    before_archives = archived_snapshot()
    manifest, graph, regions, _, records = A.preflight(root)
    groups = B.split_rows(records)
    cfg = configuration(regions, records)
    # Tests run before any fitting or R3 evaluation, and captured verbatim.
    log = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent/'tests'))
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    print(log.getvalue(), flush=True)
    if not result.wasSuccessful():
        raise ValueError('Regression tests failed; no real training started')
    output.mkdir(parents=True, exist_ok=False)
    (output/'tests.txt').write_text(log.getvalue(), encoding='utf-8')
    B.save_json(output/'frozen_config.json', cfg)
    config_hash = A.file_sha256(output/'frozen_config.json')
    print(json.dumps(cfg, indent=2), flush=True)
    print('FROZEN B2 CONFIG SHA256 '+config_hash, flush=True)
    source_hashes = {str(p.relative_to(A.WORKSPACE)): A.file_sha256(p)
                     for p in sorted(Path(__file__).parent.rglob('*')) if p.is_file() and p.suffix in ('.py', '.md')}
    B.save_json(output/'provenance.json', dict(config_sha256=config_hash, source_root=str(root), inputs=records,
               source_hashes=source_hashes, protected_archives=before_archives,
               manifest_path=str(A.MANIFEST), manifest_sha256=A.file_sha256(A.MANIFEST),
               graph_path=str(A.GRAPH), graph_sha256=A.file_sha256(A.GRAPH),
               tests_run=result.testsRun, python=sys.version, torch=torch.__version__, numpy=np.__version__,
               git_head=A.git_text('rev-parse','HEAD'), git_branch=A.git_text('branch','--show-current'),
               git_status=A.git_text('status','--short')))
    B.deterministic(B.SEEDS[0])
    layout = DynamicReadout(regions, 'mean')
    sequence = manifest['provenance']['sequence']
    training = load_windows(root, groups[2], sequence, regions, layout)
    mu = fit_mu(training)
    prior_mu = torch.load(PREVIOUS/'R2_residue_mean.pt', weights_only=True)
    if not torch.equal(mu, prior_mu):
        raise ValueError('R2 means differ from Experiment B frozen baseline')
    zero_mean = sum((item['r'].double()-mu).sum(0) for item in training)/sum(len(item['r']) for item in training)
    if zero_mean.abs().max() > 1e-12:
        raise ValueError('R2 D mean not zero')
    torch.save(mu, output/'mu_R2.pt')
    B.save_json(output/'preprocessing_checks.json', dict(mu_equal_B_bitwise=True,
               R2_D_mean_max_abs=zero_mean.abs().max().item(), fit_frames=2000, fit_replica=2))
    models, histories, replay = {}, {}, {}
    for seed in B.SEEDS:
        for mode in B.MODES:
            started = time.monotonic()
            model, history = fit(training, mu, regions, mode, seed, cfg)
            repeated, repeat_history = fit(training, mu, regions, mode, seed, cfg)
            if B.state_hash(model) != B.state_hash(repeated) or history != repeat_history:
                raise ValueError('B2 deterministic replay mismatch')
            name = f'{mode}_seed{seed}'
            torch.save(model.state_dict(), output/(name+'.pt'))
            restored = DynamicReadout(regions, mode)
            restored.load_state_dict(torch.load(output/(name+'.pt'), weights_only=True))
            if B.state_hash(restored) != B.state_hash(model):
                raise ValueError('Checkpoint roundtrip mismatch')
            models[(seed, mode)] = model
            histories[name] = history
            replay[name] = dict(bitwise_replay=True, loss_history_equal=True, checkpoint_roundtrip=True,
                                state_sha256=B.state_hash(model), checkpoint_sha256=A.file_sha256(output/(name+'.pt')))
            print(f'B2 R2 {name}: fit+replay PASS ({time.monotonic()-started:.1f}s)', flush=True)
    B.save_json(output/'training_history.json', histories)
    B.save_json(output/'reproducibility.json', replay)
    B.save_json(output/'R2_results.json', evaluate(training, models, mu, list(regions)))
    del training
    if A.file_sha256(output/'frozen_config.json') != config_hash or any(A.file_sha256(A.WORKSPACE/p) != h for p,h in source_hashes.items()):
        raise ValueError('Frozen configuration/source changed')
    B.save_json(output/'R3_evaluation_started.json', dict(config_sha256=config_hash, models_locked=replay))
    print('Beginning single B2 R3 evaluation stage', flush=True)
    heldout = load_windows(root, groups[3], sequence, regions, layout)
    B.save_json(output/'R3_results.json', evaluate(heldout, models, mu, list(regions)))
    if before_archives != archived_snapshot():
        raise ValueError('Protected A/B archive changed')
    if any(A.file_sha256(root/r['path']) != r['sha256'] for r in records):
        raise ValueError('Source NPZ changed')
    B.save_json(output/'completion.json', dict(status='COMPLETE', config_sha256=config_hash,
               tests_passed=result.testsRun, source_npz_unchanged=40, protected_A_B_archives_unchanged=True,
               exact_replays=9, R3_evaluation_stages=1))
    print('COMPLETE '+str(output), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    run(args.input_root, args.output_root)
