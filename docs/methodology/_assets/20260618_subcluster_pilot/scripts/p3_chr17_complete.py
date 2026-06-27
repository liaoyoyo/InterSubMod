#!/usr/bin/env python3
"""[chr17:48360161 完整 4-sSNV 克隆樹] 加入 longphase 輸入裡被 SEQC2 誤判 FP 的 γ=48357368 →
完整 lineage: γ_subclone(sibling) / α_only(L1) / αβ_descendant(L2) / ancestral → 完整克隆樹 +
per-lineage 甲基 + 6 對 pairwise 2×2 + normal 確認。輸出 chr17_complete_data.json。"""
import json, csv
from collections import Counter, defaultdict
import numpy as np
import pysam
import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p2_snv_linkage as P2

A = P2.A
# 4 sSNV: γ(FP-誤判,longphase輸入), β1, α(祖先), β2
SNVS = [(48357368, "C", "T", "γ"), (48362515, "G", "A", "β1"), (48365089, "G", "C", "α"), (48365161, "T", "C", "β2")]
region = [c["region_dir"] for c in json.load(open(f"{A}/cis_candidates_resolved.json")) if c["chrom"] == "chr17" and c["pos"] == "48360161"][0]

# reads.tsv + methylation
rows = open(f"{region}/reads/reads.tsv").read().splitlines(); hdr = rows[0].split("\t")
ix = {k: hdr.index(k) for k in ("read_id", "read_name", "is_tumor", "hp")}
name2rid = {}; rid_meta = {}
for r in rows[1:]:
    c = r.split("\t"); name2rid[c[ix["read_name"]]] = c[ix["read_id"]]
    rid_meta[c[ix["read_id"]]] = {"is_tumor": c[ix["is_tumor"]], "hp": c[ix["hp"]]}
mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n")
cpgs = [int(x) for x in mr[0].split(",")[1:]]
meth = {}
for ln in mr[1:]:
    q = ln.split(","); meth[q[0]] = [None if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]

tb = pysam.AlignmentFile(P2.TBAM, "rb"); nb = pysam.AlignmentFile(P2.NBAM, "rb")
tcalls, thp = P2.per_read_allele(tb, "chr17", [(p, r, a) for p, r, a, _ in SNVS])
ncalls, _ = P2.per_read_allele(nb, "chr17", [(p, r, a) for p, r, a, _ in SNVS])
tb.close(); nb.close()
POS = {p: nm for p, _, _, nm in SNVS}


def lineage(g):
    γ, α, β2 = g.get(48357368), g.get(48365089), g.get(48365161)
    if γ == "ALT" and α != "ALT":
        return "γ_subclone(sibling)"
    if α == "ALT" and β2 == "ALT":
        return "L2_αβ(descendant)"
    if α == "ALT" and β2 != "ALT":
        return "L1_α_only(ancestor)"
    if γ != "ALT" and α != "ALT":
        return "ancestral(no_somatic)"
    return "complex"


reads = []
for nm, g in tcalls.items():
    rid = name2rid.get(nm)
    if rid is None or rid not in meth or rid_meta.get(rid, {}).get("is_tumor") != "1":
        continue
    reads.append({"rid": rid, "geno": {POS[p]: g.get(p) for p in POS}, "lineage": lineage(g), "hp": rid_meta[rid]["hp"], "meth": meth[rid]})
lc = Counter(r["lineage"] for r in reads)

# 6 對 pairwise 2×2
def tab2(p1, p2):
    t = Counter()
    for nm, g in tcalls.items():
        rid = name2rid.get(nm)
        if rid is None or rid_meta.get(rid, {}).get("is_tumor") != "1":
            continue
        if g.get(p1) in ("REF", "ALT") and g.get(p2) in ("REF", "ALT"):
            t[(g[p1], g[p2])] += 1
    return {"RR": t[("REF", "REF")], "RA": t[("REF", "ALT")], "AR": t[("ALT", "REF")], "AA": t[("ALT", "ALT")]}


pos_list = [p for p, _, _, _ in SNVS]
pairs = {}
for i in range(len(pos_list)):
    for j in range(i + 1, len(pos_list)):
        pairs[f"{POS[pos_list[i]]}×{POS[pos_list[j]]}"] = tab2(pos_list[i], pos_list[j])

# normal 確認 + per-lineage 甲基
snv_norm = {}
for p, ref, alt, nm in SNVS:
    nc = Counter(c.get(p) for c in ncalls.values() if p in c)
    tc = Counter(g.get(p) for g in tcalls.values() if p in g)
    snv_norm[nm] = {"pos": p, "ref": ref, "alt": alt, "tumor_REF": tc["REF"], "tumor_ALT": tc["ALT"],
                    "normal_REF": nc["REF"], "normal_ALT": nc["ALT"], "somatic": nc["ALT"] <= 1}


def lin_meth(lab):
    arr = np.array([[np.nan if v is None else v for v in r["meth"]] for r in reads if r["lineage"] == lab], float)
    if len(arr) == 0:
        return None, 0
    return [round(float(np.nanmean(arr[:, j])), 3) if np.isfinite(arr[:, j]).any() else None for j in range(len(cpgs))], len(arr)


LABS = ["ancestral(no_somatic)", "γ_subclone(sibling)", "L1_α_only(ancestor)", "L2_αβ(descendant)"]
lin_means = {lab: lin_meth(lab) for lab in LABS}
# γ vs α-branch 甲基差; L1 vs L2 甲基差
out = {"locus": "chr17:48360161", "snvs": [{"pos": p, "ref": r, "alt": a, "name": nm} for p, r, a, nm in SNVS],
       "lineage_counts": dict(lc), "cpgs": cpgs, "n_cpg": len(cpgs), "pairs_2x2": pairs, "snv_norm": snv_norm,
       "lineage_meth": {lab: {"mean": lin_means[lab][0], "n": lin_means[lab][1]} for lab in LABS},
       "reads": reads}
json.dump(out, open(f"{A}/chr17_complete_data.json", "w"), ensure_ascii=False)
print(json.dumps({"lineage_counts": dict(lc), "pairs_2x2": pairs,
                  "snv_normal確認": {nm: f"T:{v['tumor_REF']}/{v['tumor_ALT']} N:{v['normal_REF']}/{v['normal_ALT']} somatic={v['somatic']}" for nm, v in snv_norm.items()}},
                 ensure_ascii=False, indent=1))
print("[-> chr17_complete_data.json]")
