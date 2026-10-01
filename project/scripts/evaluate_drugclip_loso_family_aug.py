#!/usr/bin/env python3
"""Family-augmented LOSO: baseline 13-target protocol (cluster0 pockets) plus
ChEMBL M2/M4 domain pairs as TRAINING-ONLY auxiliary rows.

Hygiene: aux rows never enter any test set (benchmark frozen); when target t
is held out, aux rows mapped to t are excluded from training (same-target
rule); the global held-scaffold purge applies to aux rows as well.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from finetune_drugclip_gpcr_screening import FastProjectionLoRA, binary_metrics
from finetune_drugclip_muscarinic_triplet import Proj
import copy

torch.set_num_threads(max(1, min(8, torch.get_num_threads())))


class Model(torch.nn.Module):
    def __init__(self, state, rank, seed):
        super().__init__()
        self.mol = FastProjectionLoRA(state['mol_project'], rank, seed)
        self.pocket = FastProjectionLoRA(state['pocket_project'], rank, seed + 17)


def weights(targets, labels, n_targets):
    out = np.zeros(len(labels), dtype=np.float32)
    for target in range(n_targets):
        for label in [0, 1]:
            keep = (targets == target) & (labels == label)
            if keep.any():
                out[keep] = 1 / (n_targets * 2 * keep.sum())
    return out * len(labels) / out.sum()


def fit(mol, pocket, pm, pt, y, seen, state, seed, args, random=False):
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
    train_mol = torch.unique(pmt); seen_t = torch.as_tensor(seen)
    remap = {t: i for i, t in enumerate(seen)}; active = {}
    for m, t, label in zip(pm, pt, y):
        if label: active.setdefault(int(m), set()).add(int(t))
    ret = [(m, remap[next(iter(ts))]) for m, ts in active.items() if len(ts) == 1 and next(iter(ts)) in remap]
    rm = torch.as_tensor([v[0] for v in ret]); rt = torch.as_tensor([v[1] for v in ret])
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad(); zm = model.mol.from_hidden(hm); zp = model.pocket.from_hidden(hp); score = zm @ zp.T
        bce = (F.binary_cross_entropy_with_logits(score[pmt, ptt] / args.temperature, yt, reduction='none') * wt).mean()
        retrieval = F.cross_entropy(score[rm][:, seen_t] / args.temperature, rt)
        preserve = (1 - (zm[train_mol] * z0m[train_mol]).sum(1)).mean() + (1 - (zp[seen_t] * z0p[seen_t]).sum(1)).mean()
        loss = bce + args.retrieval_weight * retrieval + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.); opt.step()
    return model


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
    p.add_argument('--family-npz', type=Path, required=True)
    p.add_argument('--family-pairs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=int, default=20260926)
    p.add_argument('--rank', type=int, default=8)
    p.add_argument('--epochs', type=int, default=80)
    p.add_argument('--lr', type=float, default=5e-3)
    p.add_argument('--weight-decay', type=float, default=1e-3)
    p.add_argument('--temperature', type=float, default=.07)
    p.add_argument('--retrieval-weight', type=float, default=.25)
    p.add_argument('--preserve-weight', type=float, default=0.)
    p.add_argument('--full-only', action='store_true')
    p.add_argument('--save-full-checkpoint', type=Path)
    a = p.parse_args()
    ga = np.load(a.gpcr_representations, allow_pickle=False); ea = np.load(a.external_representations, allow_pickle=False)
    fa = np.load(a.family_npz, allow_pickle=False)
    state = torch.load(a.projection, map_location='cpu')
    gp = pd.read_csv(a.gpcr_pairs); ep = pd.read_csv(a.external_pairs); fp = pd.read_csv(a.family_pairs)
    gids = list(map(str, ga['molecule_ids'])); eids = list(map(str, ea['molecule_ids'])); fids = list(map(str, fa['molecule_ids']))
    gmi = {v: i for i, v in enumerate(gids)}; emi = {v: i for i, v in enumerate(eids)}; fmi = {v: i for i, v in enumerate(fids)}
    gp = gp[gp.canonical_smiles.astype(str).isin(gmi)].copy(); ep = ep[ep.canonical_smiles.astype(str).isin(emi)].copy()
    fmi_map = {'M2': 'M2R', 'M4': 'M4R'}
    fp = fp[fp.subtype.isin(fmi_map) & fp.canonical_smiles.astype(str).isin(fmi)].copy()
    fp['target'] = fp.subtype.map(fmi_map)
    offset_g, offset_e = len(gids), len(eids)
    gp['mol_index'] = [gmi[v] for v in gp.canonical_smiles.astype(str)]
    ep['mol_index'] = [offset_g + emi[v] for v in ep.canonical_smiles.astype(str)]
    fp['mol_index'] = [offset_g + offset_e + fmi[v] for v in fp.canonical_smiles.astype(str)]
    gp_ids = list(map(str, ga['pocket_ids'])); gp_targets = ['B2AR', 'CCR2', 'M2R', 'M4R']
    gp_order = [gp_ids.index(f'{v}_cluster0') for v in gp_targets]
    ext_pockets = list(map(str, ea['pocket_ids']))
    targets = gp_targets + ext_pockets
    ti = {v: i for i, v in enumerate(targets)}
    gp['target_index'] = [ti[v] for v in gp.target]; ep['target_index'] = [ti[v] for v in ep.target]
    fp['target_index'] = [ti[v] for v in fp.target]
    gp['aux'] = False; ep['aux'] = False; fp['aux'] = True
    table = pd.concat([gp, ep, fp], ignore_index=True, sort=False)
    table['murcko_scaffold'] = [MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(v)) or v for v in table.canonical_smiles.astype(str)]
    mol = torch.as_tensor(np.concatenate([ga['molecule_representations'], ea['molecule_representations'],
                                          fa['molecule_representations']]).astype(np.float32))
    pocket = torch.as_tensor(np.concatenate([ga['pocket_representations'][gp_order], ea['pocket_representations']]).astype(np.float32))
    pm = table.mol_index.to_numpy(int); pt = table.target_index.to_numpy(int); labels = table.label.to_numpy(int)
    scaff = table.murcko_scaffold.astype(str).to_numpy(); aux = table.aux.to_numpy()
    frozen_m = Proj(copy.deepcopy(state['mol_project'])).eval(); frozen_p = Proj(copy.deepcopy(state['pocket_project'])).eval()
    with torch.inference_mode():
        base = (frozen_m(mol) @ frozen_p(pocket).T).numpy()
    official = base[pm, pt]
    if a.full_only:
        if a.save_full_checkpoint is None:
            raise ValueError('--full-only requires --save-full-checkpoint')
        model = fit(mol, pocket, pm, pt, labels, list(range(len(targets))), state, a.seed, a, False)
        a.save_full_checkpoint.parent.mkdir(parents=True, exist_ok=True)
        torch.save({'mol_project': model.mol.materialized_state(), 'pocket_project': model.pocket.materialized_state(),
                    'base_projection': str(a.projection), 'training_targets': targets,
                    'n_pairs': int(len(table)), 'n_aux_rows': int(aux.sum()), 'seed': a.seed,
                    'objective': 'target-balanced BCE + target retrieval + ChEMBL M2/M4 family aux (training-only)'},
                   a.save_full_checkpoint)
        print(json.dumps({'checkpoint': str(a.save_full_checkpoint), 'n_pairs': int(len(table)),
                          'n_aux': int(aux.sum()), 'targets': targets, 'seed': a.seed}, indent=2))
        return
    tuned = np.full(len(table), np.nan); random = np.full(len(table), np.nan); audit = {}
    for held, target in enumerate(targets):
        test = (pt == held) & ~aux                       # benchmark rows only
        held_scaff = set(scaff[test])
        train = (pt != held) & ~np.asarray([v in held_scaff for v in scaff])
        # same-target aux exclusion is already implied: aux of held target has pt==held
        seen = [i for i in range(len(targets)) if i != held]
        audit[target] = {'train_pairs': int(train.sum()), 'test_pairs': int(test.sum()),
                         'aux_in_train': int((train & aux).sum()),
                         'purged_rows': int(((pt != held) & ~train).sum()), 'scaffold_overlap': 0}
        for randomized, out in [(False, tuned), (True, random)]:
            model = fit(mol, pocket, pm[train], pt[train], labels[train], seen, state, a.seed + held * 1009, a, randomized)
            with torch.inference_mode():
                matrix = model.mol.from_hidden(model.mol.hidden(mol)) @ model.pocket.from_hidden(model.pocket.hidden(pocket)).T
            out[test] = matrix.numpy()[pm[test], pt[test]]
    report = {'protocol': '13-target LOSO + ChEMBL M2/M4 domain training-only augmentation '
                          '(aux rows never tested; same-target aux excluded when held; scaffold purge global)',
              'n_pairs': int(len(table)), 'n_aux': int(aux.sum()), 'n_targets': len(targets), 'targets': targets,
              'metrics': {'official': summary(table[~aux], official[~aux], targets),
                          'bce_retrieval': summary(table[~aux], tuned[~aux], targets),
                          'random_label': summary(table[~aux], random[~aux], targets)},
              'training_audit': audit,
              'hyperparameters': {k: getattr(a, k) for k in ['rank', 'epochs', 'lr', 'temperature', 'retrieval_weight', 'preserve_weight', 'seed']},
              'claim_boundary': 'Retrospective binding-screening pilot; domain ChEMBL actives are not PAM labels.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    out = table.loc[~aux, ['target', 'label', 'canonical_smiles', 'murcko_scaffold']].copy()
    out['official'] = official[~aux]; out['tuned'] = tuned[~aux]; out['random'] = random[~aux]
    out.to_csv(a.output.with_suffix('.predictions.csv'), index=False)
    print(json.dumps(report['metrics']['bce_retrieval']['macro'], indent=2))


if __name__ == '__main__':
    main()
