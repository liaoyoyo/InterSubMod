import json
d = json.load(open('shuffle_control_results.json'))
R = d['regions']
npos = sum(1 for r in R if r['real_minus_shuffle_p95'] > 0)
deltas = sorted(r['real_minus_shuffle_p95'] for r in R)
seps = sorted(abs(r['real_auc_raw'] - 0.5) for r in R)
out = []
out.append("N regions = %d" % len(R))
out.append("frac real>shuffle_p95 = %.0f%% (%d/%d)" % (100*npos/len(R), npos, len(R)))
out.append("delta(real-shuffle_p95): min=%.3f median=%.3f max=%.3f" % (deltas[0], deltas[len(deltas)//2], deltas[-1]))
out.append("|raw_auc-0.5| separation: median=%.3f min=%.3f max=%.3f" % (seps[len(seps)//2], seps[0], seps[-1]))
out.append("shuffle_sym null median=%.3f max=%.3f (真 null 基準, 非 0.5)" % (d['shuffle_auc_sym']['median'], d['shuffle_auc_sym']['max']))
out.append("real_sym median=%.3f (對稱化膨脹值, 勿當 effect size)" % d['real_auc_sym']['median'])
out.append("")
out.append("STRONG signal (top3 delta):")
for r in sorted(R, key=lambda x: -x['real_minus_shuffle_p95'])[:3]:
    out.append("  %-24s HP1=%d HP2=%d cpg=%d real=%.3f p95=%.3f d=%+.3f" %
               (r['region'], r['n_hp1'], r['n_hp2'], r['n_cpg'], r['real_auc_sym'], r['shuffle_auc_sym_p95'], r['real_minus_shuffle_p95']))
out.append("NO signal (bottom3 delta):")
for r in sorted(R, key=lambda x: x['real_minus_shuffle_p95'])[:3]:
    out.append("  %-24s HP1=%d HP2=%d cpg=%d real=%.3f p95=%.3f d=%+.3f" %
               (r['region'], r['n_hp1'], r['n_hp2'], r['n_cpg'], r['real_auc_sym'], r['shuffle_auc_sym_p95'], r['real_minus_shuffle_p95']))
open('shuffle_analysis.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
