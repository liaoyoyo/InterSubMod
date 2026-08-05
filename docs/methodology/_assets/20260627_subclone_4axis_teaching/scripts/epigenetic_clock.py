#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
樣本級 epigenetic-clock 代理(2026-06-30):normal vs tumor 全域 5mC + 甲基失序(disorder)。
無外部 clock 模型(Horvath/epiTOC 需 array 校準座標)→ 用兩個 model-free 信號:
 (1) 全域平均 β(tumor hypomethylation = 癌症標誌)
 (2) 甲基失序 frac_intermediate(0.25-0.75 的 CpG call 比例 = 隨機漂移/有絲分裂齡代理;越多=越多分裂)
固定窗(deterministic,可重現):chr1-22 每 ~30Mb 一個 20kb 窗。§13.0 compute batch。
用法: SM_TBAM=<tumor.bam> SM_NBAM=<normal.bam> SM_TAG=<sample> python3 epigenetic_clock.py
"""
import os, json, sys
import numpy as np
import pysam

TBAM = os.environ.get("SM_TBAM", "")
NBAM = os.environ.get("SM_NBAM", "")
TAG = os.environ.get("SM_TAG", "sample")
OUT = os.environ.get("SM_CLOCK_OUT", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", f"epigenetic_clock_{TAG}.json")))
MAPQ = 20; WIN = 20000; READ_CAP = 400
# 固定窗:chr1-22 每 30Mb 一窗(deterministic)
WINDOWS = []
for c in range(1, 23):
    for start in range(20_000_000, 240_000_000, 30_000_000):
        WINDOWS.append((f"chr{c}", start, start + WIN))


def scan(bam):
    """回 (mean_beta, frac_intermediate, n_calls, n_windows_used)。"""
    tb = pysam.AlignmentFile(bam, "rb")
    betas = []
    nwin = 0
    for chrom, s, e in WINDOWS:
        try:
            got = False; nread = 0
            for a in tb.fetch(chrom, s, e):
                if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ:
                    continue
                try:
                    mb = a.modified_bases
                except Exception:
                    mb = None
                if not mb:
                    continue
                for k, lst in mb.items():
                    if k[2] != "m":
                        continue
                    for _q, mlq in lst:
                        betas.append(mlq / 255.0)
                        got = True
                nread += 1
                if nread >= READ_CAP:
                    break
            if got:
                nwin += 1
        except (ValueError, OSError):
            continue
    if not betas:
        return None
    arr = np.array(betas)
    return {"mean_beta": round(float(arr.mean()), 4),
            "frac_intermediate": round(float(((arr >= 0.25) & (arr <= 0.75)).mean()), 4),
            "n_calls": len(arr), "n_windows_used": nwin}


def main():
    tum = scan(TBAM) if TBAM and os.path.exists(TBAM) else None
    nor = scan(NBAM) if NBAM and os.path.exists(NBAM) else None
    out = {"sample": TAG, "tumor": tum, "normal": nor}
    if tum and nor:
        out["delta_mean_beta(tumor-normal)"] = round(tum["mean_beta"] - nor["mean_beta"], 4)
        out["delta_disorder(tumor-normal)"] = round(tum["frac_intermediate"] - nor["frac_intermediate"], 4)
        out["interpretation"] = ("tumor hypomethylation(delta_mean_beta<0=癌症標誌) + "
                                 "disorder 增加(delta_disorder>0=更多有絲分裂/漂移=表觀『較老』)。樣本級描述,非 subclone 排序。")
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"CLOCK {TAG}:", json.dumps({k: v for k, v in out.items() if k != "interpretation"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
