#!/usr/bin/env python3
"""
Clean, reproducible per-CpG methylation discrimination test between the two
subclone lineages, with p-values — resolves the discrepancy between the two
earlier discrimination .txt files (different beta definitions) and the report
prose p-values (which were not saved to a raw file).

SINGLE, EXPLICIT lineage rule (pos3-based, the principled local marker; the
3.x CpGs are physically adjacent to pos3):
    base(pos3) == 'A'  -> alpha   (groups 1/2)
    base(pos3) == 'G'  -> beta    (groups 3/4/5)
    else               -> skip
Same rule for HKU_MOD (haplotagged) and DORADO (raw, no HP), so cross-basecaller
replication is apples-to-apples.

Per CpG x basecaller we report: n_alpha, n_beta, meanM_alpha, meanM_beta, delta,
Mann-Whitney U p, label-permutation p (20000, seed=42), flip direction.
A CpG is "robust" if significant (perm p < 0.05) AND same flip direction in BOTH
basecallers (when both have n>=3 per lineage).

Usage: methyl_discrimination_pvalues.py <hku_matrix.json> <dorado_matrix.json>
Output: prints table; writes methyl_discrimination_pvalues.txt next to inputs.
"""
import sys, json
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu

CPG_ORDER = ["2.1", "2.2", "3.1", "3.2", "3.3", "3.4", "3.5", "4.1", "5.1", "5.2"]
N_PERM = 20000
RNG = np.random.default_rng(42)


def lineage(r):
    b3 = r["bases"].get("3")
    if b3 == "A":
        return "alpha"
    if b3 == "G":
        return "beta"
    return None


def dedup(reads):
    seen, out = set(), []
    for r in reads:
        if r["name"] in seen:
            continue
        seen.add(r["name"]); out.append(r)
    return out


def perm_p(a, b):
    """Two-sided label-permutation test, statistic = |mean(a) - mean(b)|."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    obs = abs(a.mean() - b.mean())
    pool = np.concatenate([a, b]); na = len(a)
    cnt = 0
    for _ in range(N_PERM):
        RNG.shuffle(pool)
        if abs(pool[:na].mean() - pool[na:].mean()) >= obs - 1e-12:
            cnt += 1
    return (cnt + 1) / (N_PERM + 1)


def analyze(path):
    j = json.load(open(path))
    label = j["label"]
    reads = dedup(j["tumor_reads"])
    by = defaultdict(lambda: defaultdict(list))
    lin_count = defaultdict(int)
    for r in reads:
        lin = lineage(r)
        if lin is None:
            continue
        lin_count[lin] += 1
        for s in CPG_ORDER:
            m = r["meth"].get(s)
            if m and m.get("m_prob") is not None:
                by[s][lin].append(float(m["m_prob"]))
    res = {}
    for s in CPG_ORDER:
        a = by[s].get("alpha", []); b = by[s].get("beta", [])
        if len(a) >= 3 and len(b) >= 3:
            ma, mb = float(np.mean(a)), float(np.mean(b))
            try:
                mwu = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            except ValueError:
                mwu = 1.0
            pp = perm_p(a, b)
            flip = (ma >= 0.5) != (mb >= 0.5)
            res[s] = dict(na=len(a), nb=len(b), ma=ma, mb=mb, d=ma - mb,
                          mwu=mwu, pp=pp, flip=flip, ok=True)
        else:
            res[s] = dict(na=len(a), nb=len(b), ok=False)
    return label, dict(lin_count), res


def main():
    out = []
    def w(line=""):
        out.append(line); print(line)

    results = {}
    for p in sys.argv[1:]:
        label, lc, res = analyze(p)
        results[label] = res
        w(f"\n################ {label} ################")
        w(f"lineage rule: pos3==A -> alpha ; pos3==G -> beta   (single rule, both basecallers)")
        w(f"lineage read counts: " + ", ".join(f"{k}={v}" for k, v in sorted(lc.items())))
        w(f"perm: N={N_PERM} seed=42 ; MWU two-sided")
        w(f"{'CpG':>5} | {'nα':>3} {'nβ':>3} | {'meanα':>6} {'meanβ':>6} {'Δ':>6} | "
          f"{'MWU_p':>9} | {'perm_p':>8} | flip | sig(perm<.05)")
        w("-" * 84)
        for s in CPG_ORDER:
            r = res[s]
            if not r["ok"]:
                w(f"{s:>5} | {r['na']:>3} {r['nb']:>3} | {'(n<3, skip)':>30}")
                continue
            sig = "YES" if r["pp"] < 0.05 else "no"
            w(f"{s:>5} | {r['na']:>3} {r['nb']:>3} | {r['ma']:6.2f} {r['mb']:6.2f} {r['d']:+6.2f} | "
              f"{r['mwu']:9.2e} | {r['pp']:8.5f} |  {'F' if r['flip'] else '.'}   | {sig}")

    # robust = sig+same-flip in BOTH basecallers
    labels = list(results.keys())
    if len(labels) == 2:
        L1, L2 = labels
        w(f"\n################ ROBUST (sig+same-flip in BOTH {L1} & {L2}) ################")
        robust = []
        for s in CPG_ORDER:
            a = results[L1][s]; b = results[L2][s]
            if a["ok"] and b["ok"]:
                both_sig = a["pp"] < 0.05 and b["pp"] < 0.05
                same_flip = a["flip"] and b["flip"] and ((a["ma"] >= 0.5) == (b["ma"] >= 0.5))
                tag = "ROBUST" if (both_sig and same_flip) else (
                    "HKU-only" if a["pp"] < 0.05 and not b["pp"] < 0.05 else
                    "dir-mismatch" if both_sig and not same_flip else "weak")
                if both_sig and same_flip:
                    robust.append(s)
                w(f"  {s}: {L1} perm_p={a['pp']:.4f} flip={a['flip']} α={a['ma']:.2f} | "
                  f"{L2} perm_p={b['pp']:.4f} flip={b['flip']} α={b['ma']:.2f}  -> {tag}")
            else:
                w(f"  {s}: insufficient n in {'both' if not a['ok'] and not b['ok'] else (L1 if not a['ok'] else L2)}")
        w(f"\nROBUST cross-basecaller CpG set = {{{', '.join(robust)}}}  (n={len(robust)})")

    open("data/methyl_discrimination_pvalues.txt", "w").write("\n".join(out) + "\n")
    print("\n[written] data/methyl_discrimination_pvalues.txt")


if __name__ == "__main__":
    main()
