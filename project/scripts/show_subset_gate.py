import json
out = r'D:\CLC\project\results\drugclip_science2026\cgda_subset_gate_v01'
for side, label in (('old', 'OLD 2023 ckpt'), ('new', 'NEW Science-2026')):
    d = json.load(open(out + '\\%s_cgda_3seed_summary.json' % side, encoding='utf-8'))
    print('=' * 16, label, '=' * 16)
    m = d['macro']
    print('  pocket     ROC=%.4f  BEDROC=%.4f  EF0.5%%=%.2f' % (m['pocket']['roc_auc'], m['pocket']['bedroc_alpha80_5'], m['pocket']['ef0.005']))
    print('  reference  ROC=%.4f  BEDROC=%.4f  EF0.5%%=%.2f' % (m['reference']['roc_auc'], m['reference']['bedroc_alpha80_5'], m['reference']['ef0.005']))
    print('  CGDA       ROC=%.4f  BEDROC=%.4f  EF0.5%%=%.2f' % (m['cgda']['roc_auc'], m['cgda']['bedroc_alpha80_5'], m['cgda']['ef0.005']))
    print('  CGDA - reference 95%% CI:')
    for k in ('roc_auc', 'bedroc_alpha80_5', 'ef0.005', 'ef0.01', 'ef0.05'):
        v = d['delta_95ci']['reference'][k]
        sig = 'POS' if v[0] > 0 else ('NEG' if v[2] < 0 else 'n.s.')
        print('    %-20s [%+.4f, %+.4f]  %s' % (k, v[0], v[2], sig))
