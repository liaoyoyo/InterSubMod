#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
§6 開放項補算(2026-07-04):40.4% 空向量區(n_sSNV>=2 但無 read 完整覆蓋 capped 位點集→populations 空)
中,有多少可用 partial-read / pairwise co-read 救回?
定義:
  - 空向量區 E = sm_region_integration.regions 中 n_sSNV>=2 且 populations 空(2887 = 40.4%,已對帳)。
  - 可救(recoverable)= 該區有 >=1 條 pairs_eps2 pairwise co-read(某 2 位點被同 read 覆蓋)。
    → 這是 partial-read soft-likelihood 的證據來源;complete-case genotyping 丟掉但 pairwise 仍在。
  - 不可救(truly single-locus)= 無任何 pair co-read(每 read <=1 位點),= Q2「單位點稀疏」極端。
分層:recoverable-any(>=1 pair)/ chain>=3loci(>=1 pair 連 >=3 不同位點=可拼 partial 樹)/ by CN / by n_sSNV。
🔴 誠實邊界:recoverable = 有 partial 證據的「上界」,非「完整 n_sSNV 拓樸可還原」。§13.0 compute batch。
輸出 partial_read_recoverable_fraction.json。
"""
import json, os, csv
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

R = json.load(open(f"{DATA}/sm_region_integration.json"))["regions"]

def has_vec(r):
    p = r.get("populations")
    return isinstance(p, (list, dict)) and len(p) > 0

# 全區 n_sSNV>=2(=全 7143);空向量集 E
mss = [r for r in R if r.get("n_sSNV", 0) >= 2]
E = [r for r in mss if not has_vec(r)]
E_idx = {id(r): r for r in E}

# 建 per-chrom 區間(非重疊 maximal chain);記 empty 標記與 region 物件
by_chrom = defaultdict(list)
for r in R:
    if r.get("n_sSNV", 0) < 2: continue
    ch = r["region"].split(":")[0]
    by_chrom[ch].append((r["start"], r["end"], r, not has_vec(r)))
for ch in by_chrom:
    by_chrom[ch].sort()

def find_region(ch, pa, pb):
    lst = by_chrom.get(ch)
    if not lst: return None
    lo = min(pa, pb); hi = max(pa, pb)
    # 線性(每 chrom 區數不大);找 start<=lo 且 hi<=end
    for st, en, r, isempty in lst:
        if st <= lo and hi <= en: return (r, isempty)
        if st > lo: break
    return None

# 掃 pairs_eps2:每 empty 區記錄命中的位點集 + pair 數
region_loci = defaultdict(set)      # id(r) -> set of positions covered by >=1 pair
region_npairs = Counter()           # id(r) -> n pairs
n_pair_rows = 0; n_mapped = 0; n_hit_empty = 0
with open(f"{DATA}/pairs_eps2_annotated.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        n_pair_rows += 1
        ch = row["chrom"]; pa = int(row["pos_a"]); pb = int(row["pos_b"])
        hit = find_region(ch, pa, pb)
        if hit is None: continue
        n_mapped += 1
        r, isempty = hit
        if isempty:
            n_hit_empty += 1
            region_loci[id(r)].update([pa, pb])
            region_npairs[id(r)] += 1

# 分層統計
n_E = len(E)
recoverable = [r for r in E if region_npairs[id(r)] >= 1]                 # >=1 pair
chain3 = [r for r in E if len(region_loci[id(r)]) >= 3]                   # pair 連 >=3 不同位點
truly_single = [r for r in E if region_npairs[id(r)] == 0]               # 無 pair = 每 read <=1 位點

def cn_break(rs): return dict(Counter(r.get("cn") for r in rs))
def nss_break(rs):
    c = Counter(min(r.get("n_sSNV", 0), 6) for r in rs)  # 6+ 併桶
    return {(f"{k}" if k < 6 else "6+"): c[k] for k in sorted(c)}

summary = {
    "sample": "HCC1395", "partial_flag": "single-sample HCC1395 / n_sSNV>=2 全區",
    "n_regions_total(n_sSNV>=2)": len(mss),
    "n_empty_vector(populations 空)": n_E,
    "empty_vector_frac": round(n_E / len(mss), 4),
    "pairs_eps2_rows": n_pair_rows, "pairs_mapped_to_region": n_mapped, "pairs_hitting_empty_region": n_hit_empty,
    "RECOVERABLE(>=1 pairwise co-read)": len(recoverable),
    "recoverable_frac_of_empty": round(len(recoverable) / n_E, 4),
    "recoverable_frac_of_all_mss": round(len(recoverable) / len(mss), 4),
    "chain>=3loci(可拼 partial 樹)": len(chain3),
    "chain3_frac_of_empty": round(len(chain3) / n_E, 4),
    "truly_single_locus(無任何 pair,不可救)": len(truly_single),
    "truly_single_frac_of_empty": round(len(truly_single) / n_E, 4),
    "recoverable_by_cn": cn_break(recoverable),
    "truly_single_by_cn": cn_break(truly_single),
    "recoverable_by_nsSNV": nss_break(recoverable),
    "truly_single_by_nsSNV": nss_break(truly_single),
    "note": "recoverable = 有 >=1 pairwise partial 證據的『上界』,非完整 n_sSNV 拓樸可還原;chain>=3loci 才可能拼出 partial 樹。",
}
json.dump({"summary": summary}, open(f"{DATA}/partial_read_recoverable_fraction.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(summary, ensure_ascii=False, indent=1))
