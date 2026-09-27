#!/usr/bin/env python3
"""PACER-FKG corrective audit: shared-kernel signed vector magnitude vs legacy U2.

Reuses the exact frozen preprocessing, R2 node-bandwidth calibration, RFF seed
and window-block mapping of G2-C. Never modifies legacy FKG/G1/G2B/G2C files.
The new statistic is a *bag-of-residue signed mean embedding*, not the old
joint-region U-statistic. No hypothesis test or automatic replacement gate.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

import analyze_pacer_fkg_g2c as g2c

AXES = tuple(g2c.SIGNS)
REPLICAS = (2, 3)
WINDOWS = tuple(range(5))


def cosine(a, b):
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.clip(np.dot(a, b) / denom, -1., 1.)) if denom > 1e-12 else None


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return None
    def ranks(x):
        order = np.argsort(x, kind='stable')
        v = x[order]; out = np.empty(len(x), dtype=float)
        left = 0
        while left < len(x):
            right = left + 1
            while right < len(x) and v[right] == v[left]: right += 1
            out[order[left:right]] = (left + right - 1) / 2.
            left = right
        return out
    return float(np.corrcoef(ranks(a), ranks(b))[0, 1])


def safe_write_csv(path, rows):
    if not rows: raise ValueError('No data to write')
    with path.open('w', encoding='utf-8-sig', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)


def quantiles(values):
    a = np.asarray(values, dtype=float)
    return [float(np.percentile(a, p)) for p in (2.5, 50, 97.5)]


def paired_block_sensitivity(nodeblocks, ix, stable_ix, *, n_draws, seed):
    """Descriptive common-window and aligned-block resampling within each replica.

    Same window/block draws apply to region and stable core, and to both
    replicas; each draw averages signed contrast vectors before measuring norms.
    """
    rng = np.random.default_rng(seed)
    gaps = {2: [], 3: []}; cosines = []
    for _ in range(n_draws):
        wi = rng.integers(0, 5, size=5)
        bi = rng.integers(0, nodeblocks[2].shape[1], size=(5, nodeblocks[2].shape[1]))
        region_vec = {}
        for rep in REPLICAS:
            # Advanced indexing preserves (selected windows, selected blocks, nodes, RFF).
            chosen = nodeblocks[rep][wi[:, None], bi]
            rv = chosen[:, :, ix, :].mean(axis=(0, 1, 2), dtype=np.float64)
            sv = chosen[:, :, stable_ix, :].mean(axis=(0, 1, 2), dtype=np.float64)
            region_vec[rep] = rv
            gaps[rep].append(float(np.linalg.norm(rv) - np.linalg.norm(sv)))
        c = cosine(region_vec[2], region_vec[3])
        if c is not None: cosines.append(c)
    return {**{f'r{r}_gap_block_q{p}': q for r in REPLICAS
               for p, q in zip(('025', '50', '975'), quantiles(gaps[r]))},
            'cosine_block_q025': quantiles(cosines)[0] if cosines else None,
            'cosine_block_q975': quantiles(cosines)[2] if cosines else None,
            'block_draws': n_draws}


def compute(groups, graph, original, *, draws, seed):
    """Compute within-common-RFF regional norms; do not promote legacy gate."""
    region_ix = {k: g2c.get_graph_region(graph, k) for k in sorted(graph['regions'])}
    stable_ix = region_ix['stable_core_control']
    old_idx = {(s['region'], s['axis']): s for s in original['cross_replica_summary']}
    rows, windows = [], []
    for aidx, (axis, nodeblocks) in enumerate(groups.items()):
        for ridx, (name, ix) in enumerate(region_ix.items()):
            vec = {rep: nodeblocks[rep][:, :, ix, :].mean(axis=2).mean(axis=1, dtype=np.float64)
                   for rep in REPLICAS}  # per-window vectors
            stable = {rep: nodeblocks[rep][:, :, stable_ix, :].mean(axis=2).mean(axis=1, dtype=np.float64)
                      for rep in REPLICAS}
            by_rep = {}
            for rep in REPLICAS:
                norm = np.linalg.norm(vec[rep], axis=-1)
                base = np.linalg.norm(stable[rep], axis=-1)
                gap = norm - base
                by_rep[rep] = (norm, base, gap)
                for w in WINDOWS:
                    windows.append({'region': name, 'axis': axis, 'replica': rep, 'window': w,
                                    'shared_rff_norm': float(norm[w]),
                                    'stable_shared_rff_norm': float(base[w]),
                                    'shared_norm_minus_stable': float(gap[w]),
                                    'shared_norm_squared': float(norm[w] ** 2),
                                    'stable_shared_norm_squared': float(base[w] ** 2)})
            v2, v3 = vec[2].mean(axis=0), vec[3].mean(axis=0)
            orig = old_idx[(name, axis)]
            old_pass = bool(orig.get('legacy_magnitude_screen', orig.get('preregistered_specificity_gate', False)))
            # Pairwise checks in a COMMON RFF coordinate system, not different region kernels.
            local = {
                'region': name, 'axis': axis, 'residues': len(ix),
                'legacy_magnitude_screen_superseded': old_pass,
                'r2_median_shared_norm_minus_stable': float(np.median(by_rep[2][2])),
                'r3_median_shared_norm_minus_stable': float(np.median(by_rep[3][2])),
                'r2_windows_positive': int(sum(by_rep[2][2] > 0)),
                'r3_windows_positive': int(sum(by_rep[3][2] > 0)),
                'matched_window_spearman_shared_norm_gap': spearman(by_rep[2][2], by_rep[3][2]),
                'r2_pooled_shared_norm': float(np.linalg.norm(v2)),
                'r3_pooled_shared_norm': float(np.linalg.norm(v3)),
                'r2_pooled_stable_norm': float(np.linalg.norm(stable[2].mean(axis=0))),
                'r3_pooled_stable_norm': float(np.linalg.norm(stable[3].mean(axis=0))),
                'pooled_cross_replica_direction_cosine': cosine(v2, v3),
                'matched_window_cosines': [cosine(vec[2][w], vec[3][w]) for w in WINDOWS],
                'new_qualification_gate': 'NOT_DEFINED',
            }
            local.update(paired_block_sensitivity(nodeblocks, ix, stable_ix,
                         n_draws=draws, seed=seed + aidx * 1001 + ridx * 31))
            rows.append(local)
    return rows, windows


def verify_g2c(rows, existing_json, tolerance):
    """Hard-stop if implementation fails to reproduce frozen G2-C no-graph cosines."""
    if existing_json is None: return {'status': 'NOT_PROVIDED', 'count': 0}
    data = json.loads(existing_json.read_text(encoding='utf-8'))
    if data.get('G2C_STATUS') != 'COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL':
        raise ValueError('Wrong G2-C audit status')
    expected = {(v['region'], v['axis']): v['no_graph_cosine']
                for v in data['region_axis_summary']}
    if len(expected) != 27: raise ValueError('G2-C no-graph reference is incomplete')
    errors = []
    for r in rows:
        key = r['region'], r['axis']
        a, b = r['pooled_cross_replica_direction_cosine'], expected[key]
        if a is None or b is None: raise ValueError(f'Undefined cosine {key}')
        diff = abs(a - b)
        errors.append(diff)
        if diff > tolerance:
            raise AssertionError(f'Frozen G2-C reproduction FAILED {key}: {a} vs {b}; diff {diff}')
    return {'status': 'PASS', 'count': len(errors), 'max_absolute_difference': max(errors),
            'reference_sha256': g2c.sha256(existing_json), 'tolerance': tolerance}


def report_markdown(report):
    rows = report['summary']
    lines = ['# PACER-FKG 原 kernel_u2 标量筛选勘误 v02', '',
             '**结论性质：新增估计目标的描述性复核，旧门禁已撤销；不授权新门禁、药理学判断或训练。**', '',
             '## 确认的问题', '',
             '- 历史 `kernel_u2` 是 U-statistic 的平方核范数估计，不是平方根；平方本身合法。',
             '- 历史区域与 stable-core 使用独立的 per-call median 带宽，原样相减并据此设置门禁，缺乏共同核尺度。',
             '- 旧图传播扩散的是非负 `sqrt(max(kernel_u2,0))`，并非有符号 RKHS 差分向量。',
             '- 相邻 MD 帧并非独立；旧无偏表述和窗口 Spearman 不应作为统计显著性。', '',
             '## v02 修复的估计目标', '',
             '沿用冻结 G2-C 的 R2 通道预处理、共享节点 RBF 带宽、256D RFF 与种子。',
             '先分别构造每个残基的有符号四上下文核均值差，再对区域内残基取平均；',
             '该向量的范数在同一个 RFF 空间中与稳定核心范数比较。',
             '**这是共享节点核下的 bag-of-residue 估计目标，不等于旧版 joint-region U-statistic。**',
             '不根据结果设置或调节新门禁。窗口/块区间仅作敏感性分析。', '',
             '## 旧筛选结果与修正描述量对照', '',
             '| 区域 | 轴 | 历史筛选（撤销） | R2 中位范数差 | R3 中位范数差 | 新 pooled 余弦 |',
             '|---|---|---|---:|---:|---:|']
    for r in rows:
        c = r['pooled_cross_replica_direction_cosine']
        lines.append(f'| {r["region"]} | {r["axis"]} | '
                     f'{"曾通过" if r["legacy_magnitude_screen_superseded"] else "未通过"} | '
                     f'{r["r2_median_shared_norm_minus_stable"]:.5f} | '
                     f'{r["r3_median_shared_norm_minus_stable"]:.5f} | '
                     f'{c:.3f} |' if c is not None else '| N/A |')
    lines += ['', '## 校验与使用限制', '',
              f'- 输入 SHA256：{report["provenance"]["input_hashes_verified"]}/40。',
              f'- G2-C 无图 27 项交叉回放：{report["g2c_crosscheck"]["status"]}。',
              '- G2-B 的区域专用核与本次共享节点核不同，不能直接比较两者余弦大小。',
              '- 原 FKG / G1 保留为 SUPERSEDED 历史数据；G2-B / G2-C 不因本次勘误而被改写。',
              '- 长程 MD 尚未纳入；不能据此声称模型已验证或存在 PAM 功效。', '']
    return '\n'.join(lines)


def run(args):
    if args.rff_dim != 256 or args.seed != 271828 or args.block_frames != 20:
        raise ValueError('Erratum must replay G2-C frozen rff_dim=256, seed=271828, block_frames=20')
    if args.output_root.exists():
        raise FileExistsError('Refusing to overwrite output; choose a fresh v02 output directory')
    # g2c.preflight does all 40 input hashes + frozen source/checkpoint graph checks.
    graph, old, med, scale, edges, prov = g2c.preflight(args)
    if args.g2c_audit is None:
        raise ValueError('G2-C no-graph cross-check is mandatory; specify --g2c-audit')
    if not args.g2c_audit.is_file():
        raise FileNotFoundError(args.g2c_audit)
    g2c_ref = json.loads(args.g2c_audit.read_text(encoding='utf-8'))
    ref_prov = g2c_ref.get('provenance', {})
    if ref_prov.get('fkg_audit_sha256') != prov['fkg_audit_sha256'] or ref_prov.get('g2b_audit_sha256') != prov['g2b_audit_sha256'] or ref_prov.get('graph_sha256') != prov['graph_sha256']:
        raise ValueError('G2-C reference provenance does not match frozen FKG/G2B/graph SHA256')
    settings = g2c_ref.get('calibration', {})
    # Catch accidental change in alpha/other irrelevant settings as well as true RFF settings.
    if settings.get('rff_dim') != args.rff_dim or settings.get('rff_seed') != args.seed + 771:
        raise ValueError('G2-C reference RFF settings mismatch')
    if args.preflight_only:
        print(json.dumps({'ERRATUM_PREFLIGHT': 'PASS', 'inputs_hash_verified': 40,
                          'legacy_audit_sha256': prov['fkg_audit_sha256'],
                          'g2b_audit_sha256': prov['g2b_audit_sha256'],
                          'g2c_audit_sha256': g2c.sha256(args.g2c_audit),
                          'output_untouched': str(args.output_root)}, indent=2)); return
    width, calibration_samples = g2c.fit_shared_width(args, med, scale)
    if abs(width - float(settings['shared_node_bandwidth'])) > 1e-7:
        raise AssertionError('Shared bandwidth cannot reproduce frozen G2-C calibration')
    groups = g2c.load_contrasts(args, med, scale, width)
    original = json.loads(args.fkg_audit.read_text(encoding='utf-8'))
    rows, windows = compute(groups, graph, original, draws=args.bootstrap, seed=args.seed)
    cross = verify_g2c(rows, args.g2c_audit, args.crosscheck_tolerance)
    if len(rows) != 27 or len(windows) != 270 or cross['status'] != 'PASS':
        raise ValueError('Unexpected row counts / cross-validation')
    report = {'ERRATUM_STATUS': 'COMPLETED_DESCRIPTIVE_V02',
              'legacy_scalar_gate': 'SUPERSEDED_NOT_VALIDATED',
              'new_gate': 'NOT_DEFINED',
              'estimand': 'bag-of-residue signed shared-node RFF contrast mean norm',
              'provenance': {**prov, 'input_hashes_verified': 40,
                             'g2c_audit_sha256': g2c.sha256(args.g2c_audit),
                             'original_legacy_audit_sha256': g2c.sha256(args.fkg_audit)},
              'settings': {'rff_dim': args.rff_dim, 'seed': args.seed, 'width': width,
                           'calibration_samples': calibration_samples,
                           'block_frames': args.block_frames, 'descriptive_draws': args.bootstrap},
              'g2c_crosscheck': cross, 'summary': rows,
              'limitations': ['different estimand from original per-region joint U-statistic',
                              'shared finite RFF approximation, not exact infinite RKHS',
                              'two replicas and serially correlated blocks, no p-values',
                              'no training/efficacy/independent validation claim']}
    args.output_root.mkdir(parents=True, exist_ok=False)
    (args.output_root / 'FKG_U2_ERRATUM_AUDIT_v02.json').write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    safe_write_csv(args.output_root / 'FKG_U2_ERRATUM_REGION_AXIS_v02.csv',
                   [{k: v for k, v in row.items() if not isinstance(v, list)} for row in rows])
    safe_write_csv(args.output_root / 'FKG_U2_ERRATUM_WINDOWS_v02.csv', windows)
    (args.output_root / 'FKG_U2_ERRATUM_REPORT_v02.md').write_text(
        report_markdown(report), encoding='utf-8')
    print(json.dumps({'ERRATUM_STATUS': report['ERRATUM_STATUS'],
                      'G2C_NO_GRAPH_REPLAY': cross['status'],
                      'max_replay_difference': cross['max_absolute_difference'],
                      'rows': len(rows), 'window_rows': len(windows),
                      'output': str(args.output_root)}, indent=2))


def cli():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', type=Path, required=True)
    p.add_argument('--graph', type=Path, required=True)
    p.add_argument('--fkg-audit', type=Path, required=True)
    p.add_argument('--g2b-audit', type=Path, required=True)
    p.add_argument('--g2c-audit', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    # g2c preflight requires these fields even when this script does not diffuse graphs.
    p.add_argument('--rff-dim', type=int, default=256)
    p.add_argument('--seed', type=int, default=271828)
    p.add_argument('--block-frames', type=int, default=20)
    p.add_argument('--bootstrap', type=int, default=200)
    p.add_argument('--n-null', type=int, default=12)
    p.add_argument('--alpha', type=float, default=.65)
    p.add_argument('--steps', type=int, default=20)
    p.add_argument('--null-swaps', type=int, default=12000)
    p.add_argument('--crosscheck-tolerance', type=float, default=1e-5)
    p.add_argument('--preflight-only', action='store_true')
    args = p.parse_args()
    run(args)


if __name__ == '__main__': cli()
