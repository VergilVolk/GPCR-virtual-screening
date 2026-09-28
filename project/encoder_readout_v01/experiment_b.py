"""Frozen local-only Experiment B. No encoder or biological evaluation imports."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
import torch
from torch import nn

import experiment_a as A

VERSION = "encoder_readout_B_v01"
BASE = A.WORKSPACE / "project/results" / VERSION
SEEDS = (17, 43, 101)
MODES = ("mean", "static", "attention")


def identity(path):
    match = re.fullmatch(r"replica_(02|03)/window_(00[0-4])/atom14/(.+)_w(00[0-4])\.geom2vec\.npz", path)
    if not match or match[2] != match[4] or match[3] not in A.CONTEXTS:
        raise ValueError("Invalid trajectory/window identity")
    return int(match[1]), match[3], int(match[2])


def split_rows(rows):
    if len(rows) != 40 or len({r['path'] for r in rows}) != 40:
        raise ValueError("Expected 40 distinct windows")
    groups = {2: [], 3: []}
    units = {}
    for row in rows:
        replica, condition, window = identity(row['path'])
        unit = (replica, condition)
        if row.get('split', 'R' + str(replica)) != 'R' + str(replica):
            raise ValueError("Trajectory crosses development/evaluation split")
        units.setdefault(unit, []).append(window)
        groups[replica].append(row)
    if set(units) != {(r, c) for r in (2, 3) for c in A.CONTEXTS}:
        raise ValueError("Missing trajectory")
    if any(sorted(w) != list(range(5)) for w in units.values()):
        raise ValueError("Window identity constraint violated")
    train_ids = {identity(r['path'])[:2] for r in groups[2]}
    eval_ids = {identity(r['path'])[:2] for r in groups[3]}
    if train_ids & eval_ids:
        raise ValueError("Trajectory leakage")
    return {k: sorted(v, key=lambda r: identity(r['path'])) for k, v in groups.items()}


def deterministic(seed):
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(seed)
    np.random.seed(seed)


class LocalReadout(nn.Module):
    """Only visible R values reach scores or sums; no global argument exists."""
    def __init__(self, regions, mode, channels=128, residues=270):
        super().__init__()
        if mode not in MODES:
            raise ValueError("Unknown comparator")
        self.mode = mode
        self.channels = channels
        self.residues = residues
        indices, targets, memberships, region_ids = [], [], [], []
        offset = 0
        sizes = []
        for rid, ix in enumerate(regions.values()):
            ix = np.asarray(ix)
            if ix.ndim != 1 or ix.dtype.kind not in 'iu' or len(ix) < 2 or len(set(ix.tolist())) != len(ix) or np.any(ix < 0) or np.any(ix >= residues):
                raise ValueError("Invalid region membership")
            sizes.append(len(ix))
            for j, target in enumerate(ix):
                visible = [k for k in range(len(ix)) if k != j]
                indices.append([int(ix[k]) for k in visible])
                memberships.append([offset + k for k in visible])
                targets.append(int(target))
                region_ids.append(rid)
            offset += len(ix)
        if not sizes:
            raise ValueError("Empty regions")
        width = max(map(len, indices))
        valid = [[True] * len(ix) + [False] * (width - len(ix)) for ix in indices]
        # Padding repeats a visible residue, never the masked target.
        for rows in (indices, memberships):
            for row in rows:
                row.extend([row[0]] * (width - len(row)))
        for name, value, dtype in [('indices', indices, torch.long), ('membership', memberships, torch.long),
                                   ('targets', targets, torch.long), ('region_ids', region_ids, torch.long),
                                   ('valid', valid, torch.bool)]:
            self.register_buffer(name, torch.tensor(value, dtype=dtype))
        # Create identical decoder initialization across all comparators.
        self.identity_embedding = nn.Embedding(residues, 8)
        self.decoder = nn.Sequential(nn.Linear(channels + 8, 32), nn.Tanh(), nn.Linear(32, channels))
        if mode == 'static':
            self.logits = nn.Parameter(torch.zeros(offset))
        if mode == 'attention':
            self.scorer = nn.Linear(channels, 16)
            self.query = nn.Parameter(torch.zeros(len(sizes), 16))

    def pool(self, local, dropout=0.0):
        if local.ndim != 3 or len(local) == 0 or local.shape[1:] != (self.residues, self.channels) or local.dtype != torch.float32 or not torch.isfinite(local).all():
            raise ValueError("Expected finite FP32 local residual tensor")
        if dropout not in (0.0, 0.25, 0.5):
            raise ValueError("Undeclared dropout")
        visible = local[:, self.indices, :]
        valid = self.valid.clone()
        # Fixed membership-order stress test, independent of seed, data and model.
        for k in range(len(valid)):
            count = int(valid[k].sum())
            valid[k, :int(count * dropout)] = False
        if self.mode == 'attention':
            scores = (torch.tanh(self.scorer(visible)) * self.query[self.region_ids][None, :, None, :]).sum(-1)
        elif self.mode == 'static':
            scores = self.logits[self.membership][None].expand(len(local), -1, -1)
        else:
            scores = torch.zeros(visible.shape[:-1], dtype=local.dtype)
        weights = scores.masked_fill(~valid[None], -torch.inf).softmax(-1)
        return (weights[..., None] * visible).sum(-2), weights

    def forward(self, local, dropout=0.0):
        pooled, weights = self.pool(local, dropout)
        ids = self.identity_embedding(self.targets)[None].expand(len(local), -1, -1)
        return self.decoder(torch.cat((pooled, ids), -1)), weights


def local_loss(model, local):
    predicted, _ = model(local)
    return (predicted - local[:, model.targets]).square().mean()


def configuration(regions, rows):
    split = split_rows(rows)
    counts = {}
    for mode in MODES:
        deterministic(SEEDS[0])
        m = LocalReadout(regions, mode)
        counts[mode] = {'total': sum(p.numel() for p in m.parameters()),
                        'decoder': sum(p.numel() for name, p in m.named_parameters() if name.startswith(('decoder.', 'identity_embedding.')))}
    return {
        'version': VERSION, 'seeds': list(SEEDS), 'comparators': [*MODES, 'R2_residue_mean', 'zero'],
        'architecture': {'heads': 1, 'scoring_hidden': 16, 'scorer': 'shared Linear(128,16)+tanh; region queries(9,16), initialized zero',
                         'value_projection': False, 'self_attention': False, 'static': '109 membership logits initialized zero'},
        'decoder': 'shared across tasks: target-index Embedding(270,8); concatenate P; Linear(136,32), tanh, Linear(32,128)',
        'parameter_counts': counts, 'simple_predictor_parameters': {'R2_residue_mean': 270 * 128, 'zero': 0},
        'preprocessing': 'A.decompose FP64; unscaled R cast FP32 for model; no fitted scaling; simple predictor R2-only FP64 residue mean',
        'global': 'G retained by A decomposition as separate branch; not passed to model/loss; recomputable from authenticated H; no R/G archive',
        'split': {f'R{k}': [r['path'] for r in v] for k, v in split.items()},
        'trajectory_unit': 'replica x condition; five windows grouped, never split',
        'masks': 'every region membership is a target on every frame; remove target before scoring; overlaps count once per membership',
        'region_order': {k: v.tolist() for k, v in regions.items()},
        'sampling': 'all frames; condition order then window 0..4 then frame 0..99; no shuffle; batch 25 frames',
        'optimizer': {'name': 'Adam', 'lr': 0.001, 'betas': [0.9, 0.999], 'eps': 1e-8, 'weight_decay': 0},
        'schedule': 'constant', 'epochs': 5, 'stopping': 'exactly 5 epochs; no early stopping or selection; final checkpoint only',
        'loss': 'equal-weight mean squared error across frames, 109 membership targets, 128 channels of R',
        'metrics': {'mse': 'mean squared R error', 'normalized_error': 'MSE / same-subset mean(R_target^2)',
                    'relative_simple': 'MSE / same-subset R2 residue-mean MSE; improvement=1-ratio',
                    'entropy': 'mean -sum(alpha log alpha), natural log', 'effective_count': 'mean exp(per-task entropy)',
                    'aggregation': 'all frames and memberships equally weighted; window/condition cells also reported, not independent replicates',
                    'blocks': 'five ordered window blocks; aggregate across conditions, plus per-trajectory/window cells',
                    'token_dropout': [0.0, 0.25, 0.5], 'dropout_rule': 'remove floor(fraction * visible_count) earliest visible membership tokens; keep target masked; no retraining'},
        'replay': 'full independent R2 training replay for every seed/comparator; exact state and loss-history equality; checkpoint roundtrip',
        'evaluation': 'one R3 stage after all R2 fits/replays complete; all frozen models and stress tests; never refit/reselect',
        'runtime': {'device': 'cpu', 'dtype': 'float32', 'threads': 2, 'deterministic_algorithms': True},
        'limitations': ['R3 historically inspected: descriptive engineering evaluation only',
                       'temporal correlation and overlapping memberships; no frame/window independence or biological efficacy inference',
                       'R is defined using full-frame mean; visible R carries shared centering dependence on masked H; direct target R is excluded',
                       'local frame IDs preserved; window order used without inventing absolute times']}


def save_json(path, obj):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def load_local(root, row, sequence, regions):
    arrays = A.read_input(root, row, sequence)
    branches = A.decompose(arrays['residue_features'], regions)
    # Deliberately narrow the training interface; G never leaves this boundary.
    return torch.from_numpy(branches['local_residual'].astype(np.float32))


def state_hash(model):
    digest = hashlib.sha256()
    for name, value in model.state_dict().items():
        digest.update(name.encode())
        digest.update(value.detach().cpu().numpy().tobytes())
    return digest.hexdigest()


def fit(data, regions, mode, seed, cfg):
    deterministic(seed)
    model = LocalReadout(regions, mode)
    opt = torch.optim.Adam(model.parameters(), lr=cfg['optimizer']['lr'])
    history = []
    for epoch in range(cfg['epochs']):
        total = 0.0
        for _, local in data:
            for batch in local.split(25):
                opt.zero_grad(set_to_none=True)
                loss = local_loss(model, batch)
                if not torch.isfinite(loss):
                    raise ValueError('Nonfinite training loss')
                loss.backward()
                opt.step()
                total += float(loss.detach()) * len(batch)
        history.append(total / sum(len(x) for _, x in data))
    return model.eval(), history


@torch.no_grad()
def evaluate(data, models, simple):
    cells = []
    targets = next(iter(models.values())).targets
    for row, local in data:
        rep, condition, window = identity(row['path'])
        target = local[:, targets]
        energy = float(target.double().square().mean())
        simple_mse = float((target.double() - simple[targets]).square().mean())
        for seed in SEEDS:
            for mode in (*MODES, 'R2_residue_mean', 'zero'):
                for dropout in (0.0, 0.25, 0.5):
                    mse, entropy, effective = 0.0, 0.0, 0.0
                    if mode in MODES:
                        model = models[(seed, mode)]
                        for batch in local.split(25):
                            pred, weights = model(batch, dropout)
                            mse += float((pred.double() - batch[:, targets].double()).square().mean()) * len(batch)
                            ent = -(weights.double() * weights.double().clamp_min(1e-300).log()).sum(-1)
                            entropy += float(ent.mean()) * len(batch)
                            effective += float(ent.exp().mean()) * len(batch)
                        mse /= len(local)
                        entropy /= len(local)
                        effective /= len(local)
                    else:
                        mse = simple_mse if mode == 'R2_residue_mean' else energy
                        entropy = effective = None
                    cells.append(dict(replica=rep, condition=condition, window=window, seed=seed, model=mode,
                                      dropout=dropout, mse=mse, target_energy=energy, simple_mse=simple_mse,
                                      entropy=entropy, effective_count=effective))
    def aggregate(items):
        out = {k: float(np.mean([x[k] for x in items])) for k in ('mse', 'target_energy', 'simple_mse')}
        out['normalized_error'] = out['mse'] / out['target_energy'] if out['target_energy'] else None
        out['relative_simple'] = out['mse'] / out['simple_mse'] if out['simple_mse'] else None
        out['improvement_over_simple'] = 1 - out['relative_simple'] if out['relative_simple'] is not None else None
        for k in ('entropy', 'effective_count'):
            out[k] = float(np.mean([x[k] for x in items])) if items[0][k] is not None else None
        return out
    summary, blocks = [], []
    for seed in SEEDS:
        for mode in (*MODES, 'R2_residue_mean', 'zero'):
            for dropout in (0.0, 0.25, 0.5):
                selected = [x for x in cells if (x['seed'], x['model'], x['dropout']) == (seed, mode, dropout)]
                key = dict(seed=seed, model=mode, dropout=dropout)
                summary.append({**key, **aggregate(selected)})
                for window in range(5):
                    blocks.append({**key, 'window': window, **aggregate([x for x in selected if x['window'] == window])})
    return dict(summary=summary, chronological_blocks=blocks, trajectory_window_cells=cells)


def run(root, output):
    output, root = Path(output).resolve(), Path(root).resolve()
    if not BASE.resolve().is_relative_to(A.WORKSPACE) or output == BASE.resolve() or not output.is_relative_to(BASE.resolve()):
        raise ValueError('Output must be a new child of Experiment B result base')
    if output.exists():
        raise FileExistsError(output)
    if output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError('Source/output overlap')
    manifest, graph, regions, batch, records = A.preflight(root)
    groups = split_rows(records)
    cfg = configuration(regions, records)
    output.mkdir(parents=True, exist_ok=False)
    save_json(output / 'frozen_config.json', cfg)
    config_hash = A.file_sha256(output / 'frozen_config.json')
    print(json.dumps(cfg, indent=2), flush=True)
    print('FROZEN CONFIG SHA256 ' + config_hash, flush=True)
    save_json(output / 'provenance.json', dict(config_sha256=config_hash, source_root=str(root),
              inputs=records, manifest_path=str(A.MANIFEST), manifest_sha256=A.file_sha256(A.MANIFEST),
              graph_path=str(A.GRAPH), graph_sha256=A.file_sha256(A.GRAPH),
              source_hashes={p.name: A.file_sha256(p) for p in Path(__file__).parent.glob('*.py')},
              git_head=A.git_text('rev-parse', 'HEAD'), git_branch=A.git_text('branch', '--show-current'),
              python=sys.version, torch=torch.__version__, numpy=np.__version__))
    sequence = manifest['provenance']['sequence']
    # In-memory R2 cache only (~276 MB); no residual archive written. G stays in A boundary.
    training = [(r, load_local(root, r, sequence, regions)) for r in groups[2]]
    simple = sum(x.double().sum(0) for _, x in training) / sum(len(x) for _, x in training)
    torch.save(simple, output / 'R2_residue_mean.pt')
    models, replay, histories = {}, {}, {}
    for seed in SEEDS:
        for mode in MODES:
            start = time.monotonic()
            model, history = fit(training, regions, mode, seed, cfg)
            repeated, repeated_history = fit(training, regions, mode, seed, cfg)
            if state_hash(model) != state_hash(repeated) or history != repeated_history:
                raise ValueError('Deterministic R2 replay failed')
            name = f'{mode}_seed{seed}'
            checkpoint = output / (name + '.pt')
            torch.save(model.state_dict(), checkpoint)
            restored = LocalReadout(regions, mode)
            restored.load_state_dict(torch.load(checkpoint, weights_only=True))
            if state_hash(restored) != state_hash(model):
                raise ValueError('Checkpoint roundtrip failed')
            replay[name] = dict(bitwise_replay=True, loss_history_equal=True, checkpoint_roundtrip=True,
                                state_sha256=state_hash(model), checkpoint_sha256=A.file_sha256(checkpoint))
            histories[name] = history
            models[(seed, mode)] = model
            print(f'R2 {name}: fit+replay PASS ({time.monotonic()-start:.1f}s)', flush=True)
    save_json(output / 'reproducibility.json', replay)
    save_json(output / 'training_history.json', histories)
    save_json(output / 'R2_results.json', evaluate(training, models, simple))
    del training
    if A.file_sha256(output / 'frozen_config.json') != config_hash:
        raise ValueError('Frozen configuration changed')
    # Exclusive marker: failed/incomplete held-out stages cannot be silently resumed.
    save_json(output / 'R3_evaluation_started.json', {'config_sha256': config_hash, 'models_locked': replay})
    print('Beginning single locked R3 evaluation stage', flush=True)
    heldout = [(r, load_local(root, r, sequence, regions)) for r in groups[3]]
    save_json(output / 'R3_results.json', evaluate(heldout, models, simple))
    save_json(output / 'completion.json', {'status': 'COMPLETE', 'config_sha256': config_hash,
              'interpretation': 'Local residual readout information preservation only; no biological efficacy claim'})
    print('COMPLETE ' + str(output), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args()
    run(args.input_root, args.output_root)
