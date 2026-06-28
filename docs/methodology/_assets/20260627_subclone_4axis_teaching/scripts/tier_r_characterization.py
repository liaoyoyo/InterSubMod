#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#3 underpowered + #4 isolated 刻畫表(Tier-R 樹外但可刻畫):
- underpowered(有 partner 無共讀 link):CCF 分層(clonal/mid/low) + depth-limited(覆蓋可救) vs geometric(span>read,需長read)。
- isolated(read-span 內無 partner):caller VAF 分層 + 是否有 same-PS partner(Tier-PS 放單倍型)。
從 per_sSNV_census.tsv(每 sSNV: vaf/n_partners/n_realized_links/max_coread/cn_state)。
輸出 tier_r_characterization.json。compute batch(§13.0,純 census)。
"""
import json, os, csv
from collections import Counter
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MIN_COREAD = 6

def caller_af_index(loci):
    """從 caller VCF(FORMAT AF=tumor VAF)抽指定 loci 的 caller VAF。回 {locus: af}。"""
    by_chrom = {}
    for lo in loci:
        c, p = lo.split(":"); by_chrom.setdefault(c, set()).add(int(p))
    af = {}
    for c, ps in by_chrom.items():
        for src in ("tp", "fp"):
            path = f"{VCFD}/filtered_snv_{src}_{c}.vcf.gz"
            if not os.path.exists(path): continue
            try:
                tb = pysam.TabixFile(path)
                lo_hi = (min(ps), max(ps))
                for ln in tb.fetch(c, lo_hi[0] - 1, lo_hi[1] + 1):
                    f = ln.split("\t"); pos = int(f[1])
                    if pos not in ps: continue
                    fmt = f[8].split(":"); smp = f[9].split(":")
                    if "AF" in fmt:
                        try: af[f"{c}:{pos}"] = float(smp[fmt.index("AF")])
                        except Exception: pass
            except Exception: pass
    return af

def ccf_tier(v):
    if v is None: return "no_vaf"
    if v >= 0.4: return "clonal(≥0.4)"
    if v >= 0.1: return "mid(0.1-0.4)"
    return "low(<0.1,subclonal/noise)"

rows = []
with open(f"{DATA}/per_sSNV_census.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        def num(x):
            try: return float(x)
            except Exception: return None
        rows.append({"locus": r["locus"], "src": r["src"], "cn": r.get("cn_state"),
                     "vaf": num(r.get("vaf")),
                     "npart": int(r.get("n_partners_le50k") or 0),
                     "nlink": int(r.get("n_realized_links") or 0),
                     "mc": int(r.get("max_coread") or 0)})

underpowered = [r for r in rows if r["npart"] > 0 and r["nlink"] == 0]
isolated = [r for r in rows if r["npart"] == 0]
linked = [r for r in rows if r["nlink"] > 0]
# isolated census 無共讀 VAF → 補 caller VCF AF(tumor VAF)
caller_af = caller_af_index([r["locus"] for r in isolated])
for r in isolated:
    if r["vaf"] is None: r["vaf"] = caller_af.get(r["locus"])

def bucket_summary(rs, name):
    n = len(rs)
    vaf_avail = sum(1 for r in rs if r["vaf"] is not None)
    ccf = Counter(ccf_tier(r["vaf"]) for r in rs)
    src = Counter(r["src"] for r in rs)
    cn = Counter(r["cn"] for r in rs)
    # depth vs geometric(僅 underpowered 有意義:有 partner)
    depth_resc_2x = sum(1 for r in rs if r["mc"] * 2 >= MIN_COREAD)
    depth_resc_4x = sum(1 for r in rs if r["mc"] * 4 >= MIN_COREAD)
    geometric = sum(1 for r in rs if r["mc"] == 0)  # 從未共讀=span>read
    return {"n": n, "vaf_available": vaf_avail, "vaf_avail_pct": round(100*vaf_avail/n, 1) if n else 0,
            "ccf_tier": dict(ccf), "by_source": dict(src), "by_cn": dict(cn),
            "depth_rescuable_2x": depth_resc_2x, "depth_rescuable_4x": depth_resc_4x,
            "geometric(max_coread=0,需長read)": geometric}

out = {
    "universe": len(rows),
    "underpowered": {**bucket_summary(underpowered, "underpowered"),
        "note": "有 partner 但無共讀 link(coread<6);CCF 可刻畫;~57% 深覆蓋(4x)可救成 link、~43% 幾何需長 read。"},
    "isolated": {**bucket_summary(isolated, "isolated"),
        "note": "read-span≤50kb 內無 partner(樹外);caller VAF 可放 clonal 譜刻畫;可能 same-PS>50kb partner(Tier-PS 放單倍型非巢狀)。geometric/depth 對無-partner 不適用。"},
    "linked_ref": {"n": len(linked), "note": "對照:可建樹(進拓樸)"},
    "handling": {
        "underpowered": "① caller CCF 刻畫放 clonal 譜(現可) ② depth-limited 加深覆蓋救 link(~57%@4x) ③ geometric 需長 read",
        "isolated": "① caller VAF 刻畫(clonal 譜) ② Tier-PS 放單倍型(非巢狀) ③ 單樣本無法建樹但非永久 dead",
    },
}
json.dump(out, open(f"{DATA}/tier_r_characterization.json", "w"), ensure_ascii=False, indent=1)
print("TIER-R CHARACTERIZATION DONE")
print(f"underpowered n={out['underpowered']['n']} vaf%={out['underpowered']['vaf_avail_pct']} ccf={out['underpowered']['ccf_tier']}")
print(f"  depth 2x/4x={out['underpowered']['depth_rescuable_2x']}/{out['underpowered']['depth_rescuable_4x']} geometric={out['underpowered']['geometric(max_coread=0,需長read)']}")
print(f"isolated n={out['isolated']['n']} vaf%={out['isolated']['vaf_avail_pct']} ccf={out['isolated']['ccf_tier']}")
print(f"  by_source isolated={out['isolated']['by_source']}")
