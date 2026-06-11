#!/usr/bin/env python3
"""
37_copy_partition_confirm.py — nail the copy-partition confound (3 tests).

(1) within-tag/copy-aware re-measure: for each T3, compare
      d_HP        = HP1 vs HP1-1            (original; conflates copy)
      d_within    = HP1-1: alt vs ref       (TRUE somatic-allele cis, copy-controlled)
      d_allele    = all ALT vs all REF      (ALLELE-axis; also conflates copy)
      d_copy      = HP1/ref vs HP1-1/ref    (pure copy/tag, same allele)
    -> do the 60 T3 collapse (|d_within| small) while d_HP / d_allele track d_copy?
(2) ALLELE-axis copy-confound: is d_allele ~ d_copy (confounded) or ~ d_within (clean)?
(3) per-read residual-copy permutation (BRCA2): is d_within above a random-split-of-HP1-1 null?
"""
import sys, glob, json, random
import numpy as np
sys.path.insert(0, "pipeline")
from lib.cis_asm_core import load_level1_plus, collapse_modtype, _beta_cpg, MIN_N
random.seed(20260603); np.random.seed(20260603)
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
recs = json.load(open(f"{ROOT}/genome_survey_v2/cis_scan_full.json"))
t3 = [r for r in recs if r.get("cis_tier") == "T3"]


def paired(A, B, n=MIN_N):
    sh = [c for c in set(A) & set(B) if A[c][1] >= n and B[c][1] >= n]
    return (round(float(np.mean([B[c][0] - A[c][0] for c in sh])), 3), len(sh)) if len(sh) >= 4 else (None, len(sh))


def all_tags(D, p):
    return {hp for (bam, hp, al) in D.get(p, {}) if bam == "tumor"}


def measures(D, p):
    b = lambda hps, al: _beta_cpg(D, p, "tumor", set(hps), set(al) if al else None)
    tags = all_tags(D, p)
    for g, s in (("1", "1-1"), ("2", "2-1")):
        d_hp, n = paired(b({g}, None), b({s}, None))
        if d_hp is not None:
            break
    d_within, nw = paired(b({s}, {"ref"}), b({s}, {"alt"}))
    # ALLELE-axis = ALT vs REF pooled across ALL tags (alt reads live in *-1 tags)
    d_allele, na = paired(b(tags, {"ref"}), b(tags, {"alt"}))
    d_copy, nc = paired(b({g}, {"ref"}), b({s}, {"ref"}))
    return {"d_HP": d_hp, "d_within": d_within, "d_allele": d_allele, "d_copy": d_copy,
            "n_within": nw, "g": g, "s": s}


print("=== (1)+(2) 60 T3: 4 measures (where computable) ===")
print(f"{'locus':18s}{'d_HP':>8}{'d_within':>10}{'d_allele':>10}{'d_copy':>9}{'n_within':>10}")
collapsed = 0; computable = 0; allele_confound = 0
rows = []
for r in t3:
    c, p = r["locus"].split(":")
    cache = glob.glob(f"{ROOT}/pipeline/cache/level1/{c}_{p}_*.tsv.gz")
    if not cache:
        continue
    Dp = load_level1_plus(cache[0]); D = collapse_modtype(Dp)
    m = measures(D, p)
    if m["d_within"] is None:
        continue
    computable += 1
    rows.append((r["locus"], m))
    # collapse = within-tag effect much smaller than HP-axis
    if abs(m["d_within"]) < 0.5 * abs(m["d_HP"] or 1):
        collapsed += 1
    # allele-axis confounded = d_allele tracks d_copy not d_within
    if m["d_allele"] is not None and m["d_copy"] is not None:
        if abs(m["d_allele"] - (m["d_copy"] or 0)) < abs(m["d_allele"] - (m["d_within"] or 0)):
            allele_confound += 1
for loc, m in sorted(rows, key=lambda x: -abs(x[1]["d_HP"] or 0))[:15]:
    def f(x): return f"{x:.3f}" if isinstance(x, float) else str(x)
    print(f"{loc:18s}{f(m['d_HP']):>8}{f(m['d_within']):>10}{f(m['d_allele']):>10}{f(m['d_copy']):>9}{m['n_within']:>10}")
print(f"\n  computable(both alleles in somatic tag): {computable}/{len(t3)}")
print(f"  COLLAPSED (|d_within| < 0.5|d_HP|): {collapsed}/{computable}  ← cis 效應遠小於 HP-axis = copy-partition")
print(f"  ALLELE-axis copy-confounded (d_allele~d_copy not d_within): {allele_confound}/{computable}")
# summary medians
import statistics as st
dw = [abs(m["d_within"]) for _, m in rows if m["d_within"] is not None]
dhp = [abs(m["d_HP"]) for _, m in rows if m["d_HP"] is not None]
dal = [abs(m["d_allele"]) for _, m in rows if m["d_allele"] is not None]
md = lambda x: st.median(x) if x else float("nan")
print(f"\n  median |d_HP|={md(dhp):.3f}  |d_allele|={md(dal):.3f} (n={len(dal)})  |d_within|={md(dw):.3f}")
if dal and dhp:
    print(f"  -> within-tag(真cis) 是 HP-axis 的 {100*md(dw)/md(dhp):.0f}%, ALLELE-axis 的 {100*md(dw)/md(dal):.0f}%")
# write machine-readable for the HTML report
json.dump({"rows": [{"locus": l, **{k: m[k] for k in ("d_HP","d_within","d_allele","d_copy","n_within")}} for l, m in rows],
           "computable": computable, "collapsed": collapsed, "n_t3": len(t3), "allele_confound": allele_confound},
          open(f"{ROOT}/genome_survey_v2/copy_partition_confirm.json", "w"), indent=1)

print("\n=== (3) BRCA2 per-read residual-copy permutation ===")
cache = glob.glob(f"{ROOT}/pipeline/cache/level1/chr13_32315128_*.tsv.gz")[0]
Dp = load_level1_plus(cache); D = collapse_modtype(Dp); p = "32315128"
# build HP1-1 read x cpg matrix + allele label
mat = {}; allele = {}
for (bam, hp, al), cpgd in D[p].items():
    if bam != "tumor" or hp != "1-1":
        continue
    for cpg, rd in cpgd.items():
        for rid, v in rd.items():
            mat.setdefault(rid, {})[cpg] = v
            allele[rid] = al
reads = [r for r in mat if len(mat[r]) >= 5]
def alt_ref_diff(labels):
    bycpg = {}
    for rid in reads:
        for cpg, v in mat[rid].items():
            bycpg.setdefault(cpg, {"alt": [], "ref": []})
            if labels[rid] in ("alt", "ref"):
                bycpg[cpg][labels[rid]].append(1 if v >= 0.5 else 0)
    diffs = [np.mean(d["alt"]) - np.mean(d["ref"]) for d in bycpg.values() if len(d["alt"]) >= MIN_N and len(d["ref"]) >= MIN_N]
    return np.mean(diffs) if len(diffs) >= 4 else None
obs = alt_ref_diff(allele)
alt_n = sum(1 for r in reads if allele[r] == "alt"); ref_n = sum(1 for r in reads if allele[r] == "ref")
null = []
labs = [allele[r] for r in reads]
for _ in range(2000):
    sh = labs[:]; random.shuffle(sh)
    nl = {reads[i]: sh[i] for i in range(len(reads))}
    d = alt_ref_diff(nl)
    if d is not None: null.append(d)
null = np.array(null)
pval = (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1) if obs is not None else None
print(f"  HP1-1 reads: alt={alt_n} ref={ref_n}")
print(f"  observed within-tag alt-ref = {obs:.4f}")
print(f"  random-split null: mean={null.mean():.4f} sd={null.std():.4f} 95%=[{np.percentile(null,2.5):.3f},{np.percentile(null,97.5):.3f}]")
print(f"  permutation p (obs vs random HP1-1 split) = {pval:.3f}")
verdict3 = ('within-tag 在 random-split null 內 = 非 allele-specific, 仍是 within-tag/copy 異質性'
            if pval and pval > 0.05 else 'within-tag 顯著 = 真有小 allele 效應')
print(f"  -> {verdict3}")
import json as _j
d = _j.load(open(f"{ROOT}/genome_survey_v2/copy_partition_confirm.json"))
d["brca2_perm"] = {"obs": round(float(obs), 4), "null_mean": round(float(null.mean()), 4),
                   "null_sd": round(float(null.std()), 4), "p": round(float(pval), 3),
                   "alt_n": alt_n, "ref_n": ref_n, "verdict": verdict3}
_j.dump(d, open(f"{ROOT}/genome_survey_v2/copy_partition_confirm.json", "w"), indent=1)
