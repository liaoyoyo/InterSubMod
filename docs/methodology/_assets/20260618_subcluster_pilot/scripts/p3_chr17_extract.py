#!/usr/bin/env python3
"""[chr17:48360161 nested subclone 詳細抽取] per-read: 3 sSNV 基因型 + 甲基向量 + HP + tumor/normal →
按 lineage 分組(L0 89-REF / L1 89-ALT,161-REF 祖先 / L2 89-ALT,161-ALT 後代) → per-lineage per-CpG 甲基 +
L1-vs-L2 差異(非循環: lineage 由 sSNV 定義非甲基)。join BAM read_name ↔ reads.tsv read_id ↔ methylation.csv。
輸出 chr17_subclone_data.json。"""
import sys, json
from collections import defaultdict, Counter
import numpy as np
import pysam
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p2_snv_linkage as P2

A = P2.A
KEY = "chr17:48360161"
CHROM = "chr17"
SNVS = [(48362515, "G", "A"), (48365089, "G", "C"), (48365161, "T", "C")]  # β1, α(ancestral), β2
RD = {f"{c['chrom']}:{c['pos']}:{c['set']}": c["region_dir"] for c in json.load(open(f"{A}/cis_candidates_resolved.json"))}
region = RD[f"{KEY}:TP"]

# ---- reads.tsv: read_id ↔ read_name, is_tumor, hp ----
rid2name = {}; rid_meta = {}
rows = open(f"{region}/reads/reads.tsv").read().splitlines(); hdr = rows[0].split("\t")
ix = {k: hdr.index(k) for k in ("read_id", "read_name", "is_tumor", "hp")}
name2rid = {}
for r in rows[1:]:
    c = r.split("\t")
    rid, nm = c[ix["read_id"]], c[ix["read_name"]]
    rid2name[rid] = nm; name2rid[nm] = rid
    rid_meta[rid] = {"is_tumor": c[ix["is_tumor"]], "hp": c[ix["hp"]]}

# ---- methylation.csv: read_id → β vector ----
mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n")
cpgs = [int(x) for x in mr[0].split(",")[1:]]
meth = {}
for ln in mr[1:]:
    q = ln.split(",")
    meth[q[0]] = [None if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]

# ---- BAM per-read allele at 3 sSNV (tumor + normal) ----
tb = pysam.AlignmentFile(P2.TBAM, "rb"); nb = pysam.AlignmentFile(P2.NBAM, "rb")
tcalls, thp = P2.per_read_allele(tb, CHROM, SNVS)
ncalls, _ = P2.per_read_allele(nb, CHROM, SNVS)
tb.close(); nb.close()
P = [s[0] for s in SNVS]


def lineage(g):
    a89 = g.get(48365089); a161 = g.get(48365161)
    if a89 == "REF":
        return "L0_ancestral_root"          # 無 α
    if a89 == "ALT" and a161 == "REF":
        return "L1_alpha_only(ancestor)"     # α 祖先
    if a89 == "ALT" and a161 == "ALT":
        return "L2_alpha_beta(descendant)"   # α+β 後代
    return "other"


# ---- join: tumor reads with BAM genotype + methylation ----
reads = []
for nm, g in tcalls.items():
    rid = name2rid.get(nm)
    if rid is None or rid not in meth:
        continue
    if rid_meta.get(rid, {}).get("is_tumor") != "1":
        continue
    lin = lineage(g)
    reads.append({"rid": rid, "geno": {str(p): g.get(p) for p in P}, "lineage": lin,
                  "hp": rid_meta[rid]["hp"], "meth": meth[rid]})
lin_counts = Counter(r["lineage"] for r in reads)


# ---- per-lineage per-CpG 甲基 mean (tumor) ----
def lineage_cpg_mean(lin):
    arr = np.array([[np.nan if v is None else v for v in r["meth"]] for r in reads if r["lineage"] == lin], float)
    if len(arr) == 0:
        return [None] * len(cpgs), 0
    return [round(float(np.nanmean(arr[:, j])), 3) if np.isfinite(arr[:, j]).any() else None for j in range(len(cpgs))], len(arr)


L1_mean, n1 = lineage_cpg_mean("L1_alpha_only(ancestor)")
L2_mean, n2 = lineage_cpg_mean("L2_alpha_beta(descendant)")
# L1 vs L2 per-CpG 差異 (非循環: 由 sSNV 定義)
diff = []
for j in range(len(cpgs)):
    if L1_mean[j] is not None and L2_mean[j] is not None:
        diff.append({"cpg": cpgs[j], "L1": L1_mean[j], "L2": L2_mean[j], "dbeta": round(L2_mean[j] - L1_mean[j], 3)})
sig_diff = [d for d in diff if abs(d["dbeta"]) >= 0.2]

# ---- normal 甲基 (這些 CpG, 確認 tumor-specific) ----
normal_meth_n = sum(1 for rid, m in rid_meta.items() if m["is_tumor"] == "0" and rid in meth)
norm_arr = np.array([[np.nan if v is None else v for v in meth[rid]]
                     for rid, m in rid_meta.items() if m["is_tumor"] == "0" and rid in meth], float)
norm_mean = [round(float(np.nanmean(norm_arr[:, j])), 3) if len(norm_arr) and np.isfinite(norm_arr[:, j]).any() else None
             for j in range(len(cpgs))] if len(norm_arr) else [None] * len(cpgs)

# ---- normal sSNV 確認 (REF) ----
snv_norm = {}
for p, ref, alt in SNVS:
    nc = Counter(c.get(p) for c in ncalls.values() if p in c)
    snv_norm[p] = {"ref": ref, "alt": alt, "normal_REF": nc["REF"], "normal_ALT": nc["ALT"]}

lk = json.load(open(f"{A}/p2_linkage.json"))[KEY]
out = {"locus": KEY, "chrom": CHROM, "snvs": [{"pos": p, "ref": r, "alt": a} for p, r, a in SNVS],
       "cpgs": cpgs, "n_cpg": len(cpgs),
       "lineage_counts": dict(lin_counts), "n_tumor_reads_genotyped": len(reads),
       "pairs_2x2": lk["pairs"], "snv_stat": lk["snv_stat"], "snv_normal": snv_norm,
       "L1_n": n1, "L2_n": n2, "L1_cpg_mean": L1_mean, "L2_cpg_mean": L2_mean, "normal_cpg_mean": norm_mean,
       "n_normal_meth_reads": normal_meth_n,
       "L1_vs_L2_diff": diff, "n_sig_diff_cpg": len(sig_diff), "sig_diff_cpg": sig_diff[:30],
       "reads": reads}
json.dump(out, open(f"{A}/chr17_subclone_data.json", "w"), ensure_ascii=False)
print(json.dumps({"locus": KEY, "lineage_counts": dict(lin_counts), "n_cpg": len(cpgs),
                  "L1_n(ancestor)": n1, "L2_n(descendant)": n2, "n_sig_diff_cpg(L1vsL2,|db|>=0.2)": len(sig_diff),
                  "snv_normal": {p: v for p, v in snv_norm.items()}}, ensure_ascii=False, indent=1))
print("[-> chr17_subclone_data.json]")
