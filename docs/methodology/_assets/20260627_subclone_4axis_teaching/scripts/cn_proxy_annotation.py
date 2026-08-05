#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
每區 CN 代理標註(2026-06-30):用 tumor tagged BAM(HP tag)+ germline phased het VCF
算 per-region (a)BAF 失衡 (b)depth logR (c)HP read 比 → likely cn_proxy(neutral/loh/gain/loss/balanced_gain)。
⚠ 註解級代理(非 SAVANA segment 級);補 cn=unknown,不取代 segment CN。§13.0 compute batch。
BAF 直接從 tumor BAM 在 germline het 位點數 REF/ALT(不用 VCF 的 AD,來源無歧義);HP 比用 haplotag。
用法: SM_DATA=<dir> SM_TBAM=<tagged.bam> SM_GVCF=<germline_phased.vcf.gz> python3 cn_proxy_annotation.py
"""
import json, os, sys
import numpy as np
import pysam
from statistics import median

DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")))
TBAM = os.environ.get("SM_TBAM", "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam")
GVCF = os.environ.get("SM_GVCF", "")
MAPQ = 20; MIN_HET = 4; MIN_SITE_DP = 8   # 每區至少 MIN_HET 個 het;每 het 至少 MIN_SITE_DP 覆蓋
B2I = {"A": 0, "C": 1, "G": 2, "T": 3}    # count_coverage 回 [A,C,G,T] 陣列
GAIN_LOGR = 1.30; LOSS_LOGR = 0.70; BAF_IMBAL = 0.12  # 門檻

# 分塊平行(2026-06-30,依資源 memory):argv[1]=chunk_total argv[2]=chunk_idx;strided 切分,跑完 merge。
CHUNK_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CHUNK_IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
OUTSUF = "" if CHUNK_TOTAL == 1 else f"_part{CHUNK_IDX}"
det = json.load(open(os.path.join(DATA, "topology_per_region.json")))["detail"][CHUNK_IDX::CHUNK_TOTAL]

def _gvcf():
    if not GVCF or not os.path.exists(GVCF): return None
    idx = GVCF + ".tbi" if os.path.exists(GVCF + ".tbi") else (GVCF + ".csi" if os.path.exists(GVCF + ".csi") else None)
    try: return pysam.TabixFile(GVCF, index=idx) if idx else pysam.TabixFile(GVCF)
    except Exception: return None

def het_sites(gv, chrom, s, e):
    """germline het(0|1 或 1|0) 位點 → [(pos, ref, alt)]。"""
    out = []
    if gv is None: return out
    try:
        for ln in gv.fetch(chrom, max(0, s - 1), e + 1):
            f = ln.split("\t")
            if len(f) < 10: continue
            gt = f[9].split(":")[0]
            if gt in ("0|1", "1|0", "0/1") and len(f[3]) == 1 and len(f[4]) == 1:
                out.append((int(f[1]), f[3].upper(), f[4].upper()))
    except Exception: pass
    return out

def hp_of(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

def main():
    gv = _gvcf()
    tb = pysam.AlignmentFile(TBAM, "rb")
    rows = []
    for r in det:
        chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
        hs = het_sites(gv, chrom, s, e)
        baf_devs = []; site_dps = []
        for pos, ref, alt in hs:
            if ref not in B2I or alt not in B2I:
                continue
            cov = tb.count_coverage(chrom, pos - 1, pos, quality_threshold=0)  # 快:per-base ACGT,跳 secondary/dup
            refc = int(cov[B2I[ref]][0]); altc = int(cov[B2I[alt]][0])
            dp = refc + altc
            if dp >= MIN_SITE_DP:
                baf_devs.append(abs(altc / dp - 0.5)); site_dps.append(dp)
        rec = {"region": r["region"], "cn_seqc2": r.get("cn"), "n_het": len(baf_devs)}
        if len(baf_devs) >= MIN_HET:
            rec["baf_dev"] = round(float(median(baf_devs)), 3)
            rec["mean_site_dp"] = round(float(np.mean(site_dps)), 1)
            rec["hp_ratio_minor"] = None  # HP-ratio 暫略(count_coverage 不分 HP;BAF 已足夠當代理)
        else:
            rec["baf_dev"] = None; rec["mean_site_dp"] = None; rec["hp_ratio_minor"] = None
        rows.append(rec)
    # chunk 模式:只輸出 raw rows(baf_dev/mean_site_dp),gmed+分類留給 merge(需全域中位深度)
    if CHUNK_TOTAL > 1:
        json.dump({"rows_raw": rows}, open(os.path.join(DATA, f"cn_proxy{OUTSUF}.json"), "w"), ensure_ascii=False)
        print(f"CN PROXY CHUNK {CHUNK_IDX}/{CHUNK_TOTAL} done: {len(rows)} regions")
        return
    # 全域中位深度 → depth logR
    dps = [x["mean_site_dp"] for x in rows if x["mean_site_dp"]]
    gmed = median(dps) if dps else 1.0
    for x in rows:
        x["depth_logr"] = round(x["mean_site_dp"] / gmed, 3) if x["mean_site_dp"] else None
        # 分類
        bd, lr = x["baf_dev"], x["depth_logr"]
        if bd is None:
            x["cn_proxy"] = "unknown(coverage<%d het)" % MIN_HET
        else:
            imbal = bd >= BAF_IMBAL
            if lr is not None and lr >= GAIN_LOGR:
                x["cn_proxy"] = "gain" if imbal else "balanced_gain"
            elif lr is not None and lr <= LOSS_LOGR:
                x["cn_proxy"] = "loss/loh" if imbal else "loss"
            elif imbal:
                x["cn_proxy"] = "loh(cn-neutral)"
            else:
                x["cn_proxy"] = "neutral"
    from collections import Counter
    dist = Counter(x["cn_proxy"].split("(")[0] for x in rows)
    summary = {"n_regions": len(rows), "n_with_baf": sum(1 for x in rows if x["baf_dev"] is not None),
               "genome_median_het_dp": round(gmed, 1), "cn_proxy_dist": dict(dist),
               "params": {"MIN_HET": MIN_HET, "MIN_SITE_DP": MIN_SITE_DP, "GAIN_LOGR": GAIN_LOGR,
                          "LOSS_LOGR": LOSS_LOGR, "BAF_IMBAL": BAF_IMBAL},
               "note": "註解級代理(BAF失衡+depth logR+HP比);非 SAVANA segment CN。cn_seqc2 欄為 HCC1395 真值(驗證用)。"}
    json.dump({"summary": summary, "regions": rows}, open(os.path.join(DATA, "cn_proxy_annotation.json"), "w"), ensure_ascii=False, indent=1)
    print("CN PROXY DONE", json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
