#!/usr/bin/env python3
"""Multi-conformation GPCR pocket variant of the 13/20-target LOSO protocol.

Differences vs evaluate_drugclip_13target_loso.py (which stays frozen):
- All 10 GaMD cluster conformations per GPCR target enter the pocket matrix
  (40 GPCR + 16 external columns).
- Target-level score S[i, t] = max over target t's pocket columns; every loss
  (BCE, retrieval, preservation baseline unchanged on pockets) and every
  evaluation uses S. This integrates the published GaMD ensemble (our MD data)
  into fine-tuning instead of the cluster0-only pocket.
"""
from __future__ import annotations
import argparse, copy, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from finetune_drugclip_gpcr_screening import FastProjectionLoRA, binary_metrics
from finetune_drugclip_muscarinic_triplet import Proj

torch.set_num_threads(max(1, min(8, torch.get_num_threads())))


class Model(torch.nn.Module):
    def __init__(self, state, rank, seed):
        super().__init__()
        self.mol = FastProjectionLoRA(state['mol_project'], rank, seed)
        self.pocket = FastProjectionLoRA(state['pocket_project'], rank, seed + 17)


def weights(pt, y, n_slots):
    out = np.zeros(len(y), dtype=np.float32)
    for target in np.unique(pt):
        keep = pt == target
        pos = keep & (y == 1); neg = keep & (y == 0)
        if pos.sum() and neg.sum():
            out[pos] = 0.5 / pos.sum(); out[neg] = 0.5 / neg.sum()
        else:
            out[keep] = 1.0 / keep.sum()
    return out * len(y) / out.sum()


def fit(mol, pocket, pm, pt, y, seen, state, seed, args, groups, random=False):
    model = Model(state, args.rank, seed)
    with torch.no_grad():
        hm = model.mol.hidden(mol); hp = model.pocket.hidden(pocket)
        z0m = model.mol.base_from_hidden(hm); z0p = model.pocket.base_from_hidden(hp)
    y = np.asarray(y, dtype=int).copy()
    if random:
        rng = np.random.default_rng(seed)
        for target in seen:
            keep = np.flatnonzero(pt == target); y[keep] = rng.permutation(y[keep])
    pmt = torch.as_tensor(np.asarray(pm, dtype=np.int64).copy()); ptt = torch.as_tensor(np.asarray(pt, dtype=np.int64).copy())
    yt = torch.as_tensor(y, dtype=torch.float32); wt = torch.as_tensor(weights(pt, y, pocket.shape[0]))
    train_mol = torch.unique(pmt)
    remap = {t: i for i, t in enumerate(seen)}; active = {}
    for m, t, label in zip(pm, pt, y):
        if label: active.setdefault(int(m), set()).add(int(t))
    ret = [(m, remap[next(iter(ts))]) for m, ts in active.items() if len(ts) == 1 and next(iter(ts)) in remap]
    rm = torch.as_tensor([v[0] for v in ret]); rt = torch.as_tensor([v[1] for v in ret])
    seen_pockets = torch.as_tensor(np.concatenate([groups[t] for t in seen]))
    seg = torch.as_tensor(np.concatenate([[t] * len(groups[t]) for t in range(len(groups))]))
    seen_pockets = torch.as_tensor(np.concatenate([groups[t] for t in seen]))
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad()
        zm = model.mol.from_hidden(hm); zp = model.pocket.from_hidden(hp)
        score = zm @ zp.T                                   # n_mol x n_pockets
        if args.reduce == 'max':
            S = torch.full((score.shape[0], len(groups)), -1e9, dtype=score.dtype)
            S = S.scatter_reduce(1, seg.view(1, -1).expand_as(score), score, reduce='amax', include_self=False)
        else:
            S = torch.zeros(score.shape[0], len(groups), dtype=score.dtype)
            S = S.scatter_reduce(1, seg.view(1, -1).expand_as(score), score, reduce='mean', include_self=False)
        bce = (F.binary_cross_entropy_with_logits(S[pmt, ptt] / args.temperature, yt, reduction='none') * wt).mean()
        retrieval = F.cross_entropy(S[rm][:, seen] / args.temperature, rt)
        preserve = (1 - (zm[train_mol] * z0m[train_mol]).sum(1)).mean() + \
                   (1 - (zp[seen_pockets] * z0p[seen_pockets]).sum(1)).mean()
        loss = bce + args.retrieval_weight * retrieval + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.); opt.step()
    return model


def target_scores(fn_mol, fn_pocket, mol, pocket, groups, reduce='max'):
    with torch.inference_mode():
        zm = fn_mol(mol); zp = fn_pocket(pocket)
        score_all = zm @ zp.T                               # n_mol x n_pockets
        idx = torch.as_tensor(np.concatenate([groups[t] for t in range(len(groups))]))
        seg = torch.as_tensor(np.concatenate([[t] * len(groups[t]) for t in range(len(groups))]))
        score = score_all[:, idx]
        S = torch.full((score.shape[0], len(groups)), -1e9, dtype=score.dtype)
        if reduce == 'max':
            S = S.scatter_reduce(1, seg.view(1, -1).expand_as(score), score, reduce='amax', include_self=False)
        else:
            S = torch.zeros(score.shape[0], len(groups), dtype=score.dtype)
            S = S.scatter_reduce(1, seg.view(1, -1).expand_as(score), score, reduce='mean', include_self=False)
    return S.numpy()


def summary(table, scores, targets):
    per = {}
    for target in targets:
        keep = table.target.eq(target).to_numpy()
        per[target] = binary_metrics(table.loc[keep, 'label'].to_numpy(int), scores[keep])
    names = ['roc_auc', 'pr_auc', 'bedroc_alpha20', 'ef1pct', 'ef5pct']
    return {'macro': {n: float(np.mean([per[t][n] for t in targets])) for n in names}, 'per_target': per}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gpcr-representations', type=Path, required=True)
    p.add_argument('--gpcr-pairs', type=Path, required=True)
    p.add_argument('--external-representations', type=Path, required=True)
    p.add_argument('--external-pairs', type=Path, required=True)
    p.add_argument('--projection', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=int, default=20260926)
    p.add_argument('--rank', type=int, default=8)
    p.add_argument('--epochs', type=int, default=80)
    p.add_argument('--lr', type=float, default=5e-3)
    p.add_argument('--weight-decay', type=float, default=1e-3)
    p.add_argument('--temperature', type=float, default=.07)
    p.add_argument('--retrieval-weight', type=float, default=.25)
    p.add_argument('--preserve-weight', type=float, default=0.)
    p.add_argument('--reduce', choices=['max', 'mean'], default='max')
    p.add_argument('--eval-cluster0', action='store_true',
                   help='evaluate on cluster0-only pockets while training on the full GaMD ensemble')
    a = p.parse_args()
    ga = np.load(a.gpcr_representations, allow_pickle=False); ea = np.load(a.external_representations, allow_pickle=False)
    state = torch.load(a.projection, map_location='cpu')
    gp = pd.read_csv(a.gpcr_pairs); ep = pd.read_csv(a.external_pairs)
    gids = list(map(str, ga['molecule_ids'])); eids = list(map(str, ea['molecule_ids']))
    gmi = {v: i for i, v in enumerate(gids)}; emi = {v: i for i, v in enumerate(eids)}
    gp = gp[gp.canonical_smiles.astype(str).isin(gmi)].copy(); ep = ep[ep.canonical_smiles.astype(str).isin(emi)].copy()
    offset = len(gids)
    gp['mol_index'] = [gmi[v] for v in gp.canonical_smiles.astype(str)]; ep['mol_index'] = [offset + emi[v] for v in ep.canonical_smiles.astype(str)]
    gp_ids = list(map(str, ga['pocket_ids']))  # e.g. B2AR_cluster0..9 (GaMD clusters)
    gp_targets = ['B2AR', 'CCR2', 'M2R', 'M4R']
    groups = []
    for t in gp_targets:
        groups.append([gp_ids.index(v) for v in gp_ids if v.startswith(t + '_')])
    ext_pockets = list(map(str, ea['pocket_ids']))
    n_g_pockets = len(gp_ids)
    groups += [[n_g_pockets + i] for i in range(len(ext_pockets))]
    targets = gp_targets + ext_pockets
    ti = {v: i for i, v in enumerate(targets)}
    gp['target_index'] = [ti[v] for v in gp.target]; ep['target_index'] = [ti[v] for v in ep.target]
    table = pd.concat([gp, ep], ignore_index=True, sort=False)
    table['murcko_scaffold'] = [MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(v)) or v for v in table.canonical_smiles.astype(str)]
    mol = torch.as_tensor(np.concatenate([ga['molecule_representations'], ea['molecule_representations']]).astype(np.float32))
    pocket = torch.as_tensor(np.concatenate([ga['pocket_representations'], ea['pocket_representations']]).astype(np.float32))
    pm = table.mol_index.to_numpy(int); pt = table.target_index.to_numpy(int); labels = table.label.to_numpy(int)
    scaff = table.murcko_scaffold.astype(str).to_numpy()
    frozen_m = Proj(copy.deepcopy(state['mol_project'])).eval(); frozen_p = Proj(copy.deepcopy(state['pocket_project'])).eval()
    base_groups = [list(map(int, g)) for g in groups]
    if a.eval_cluster0:
        eval_groups = [[gp_ids.index(f'{t}_cluster0')] for t in gp_targets] + \
                      [[n_g_pockets + i] for i in range(len(ext_pockets))]
    else:
        eval_groups = base_groups
    S0 = target_scores(frozen_m, frozen_p, mol, pocket, eval_groups, a.reduce)
    official = S0[pm, pt]
    tuned = np.full(len(table), np.nan); random = np.full(len(table), np.nan); audit = {}
    for held, target in enumerate(targets):
        test = pt == held; held_scaff = set(scaff[test])
        train = (pt != held) & ~np.asarray([v in held_scaff for v in scaff])
        seen = [i for i in range(len(targets)) if i != held]
        audit[target] = {'train_pairs': int(train.sum()), 'test_pairs': int(test.sum()),
                         'purged_rows': int(((pt != held) & ~train).sum()), 'scaffold_overlap': 0}
        for randomized, out in [(False, tuned), (True, random)]:
            model = fit(mol, pocket, pm[train], pt[train], labels[train], seen, state, a.seed + held * 1009, a, base_groups, randomized)
            S = target_scores(lambda x: model.mol.from_hidden(model.mol.hidden(x)),
                              lambda x: model.pocket.from_hidden(model.pocket.hidden(x)),
                              mol, pocket, eval_groups, a.reduce)
            out[test] = S[pm[test], pt[test]]
    report = {'protocol': 'multi-conformation GaMD pocket LOSO (10 clusters/GPCR target, target score = max over conformations); '
                          'held-target scaffolds purged; held pocket absent from loss',
              'n_pairs': len(table), 'n_targets': len(targets), 'targets': targets,
              'metrics': {'official': summary(table, official, targets), 'bce_retrieval': summary(table, tuned, targets),
                          'random_label': summary(table, random, targets)},
              'training_audit': audit,
              'hyperparameters': {k: getattr(a, k) for k in ['rank', 'epochs', 'lr', 'temperature', 'retrieval_weight', 'preserve_weight', 'seed']},
              'claim_boundary': 'Retrospective target-transfer pilot with GaMD ensemble pockets; no efficacy claim.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    out = table[['target', 'label', 'canonical_smiles', 'murcko_scaffold']].copy()
    out['official'] = official; out['tuned'] = tuned; out['random'] = random
    out.to_csv(a.output.with_suffix('.predictions.csv'), index=False)
    print(json.dumps(report['metrics']['bce_retrieval']['macro'], indent=2))


if __name__ == '__main__':
    main()
