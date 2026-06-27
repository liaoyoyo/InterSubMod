#!/usr/bin/env python3
"""[SEQC2 CN/LOH 標註] 把全 34736 位點(本 session full run)用 SEQC2 NGS-benchmark CNV 真值標 CN 狀態+值+LOH。
來源(hg38, first-hand): /big8_disk/data/HCC1395/SEQC2/CNV/
  ngs_benchmark_cnvs_gain_loss_loh.bed (chrom start end state{gain/loss/loh}) — 660 區
  ngs_benchmark_cnv_gain_cn.bed / _loss_cn.bed (chrom start end CN值)
每位點 SNV pos = rec_pos + 5000 → 區間重疊 → cn_state(gain/loss/loh/neutral) + cn_value + is_loh。
輸出 phylo_cpp_wg_full_records_cn.json(records+cn) + cn_stats.json。純讀真值,不捏造。"""
import json, os
from bisect import bisect_right

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
CNV = "/big8_disk/data/HCC1395/SEQC2/CNV"


def load_bed(path, valcol=3, asint=False):
    """→ {chrom: (sorted_starts[], regions[(start,end,val)])}."""
    by = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        v = (round(float(f[valcol]), 2) if asint else f[valcol])   # CN 值可為小數(median CN, 如 4.5)
        by.setdefault(c, []).append((s, e, v))
    out = {}
    for c, regs in by.items():
        regs.sort()
        out[c] = ([r[0] for r in regs], regs)
    return out


def overlap(idx, chrom, pos):
    """回該 pos 落入的 region val, 無則 None。"""
    if chrom not in idx:
        return None
    starts, regs = idx[chrom]
    i = bisect_right(starts, pos) - 1
    # 往回找(區間可能不互斥, 線性回掃少量)
    while i >= 0 and i >= bisect_right(starts, pos) - 8:
        s, e, v = regs[i]
        if s <= pos < e:
            return v
        i -= 1
    # 保險: 線性掃該 chrom(660 區總量小)
    for s, e, v in regs:
        if s <= pos < e:
            return v
    return None


def main():
    state_idx = load_bed(f"{CNV}/ngs_benchmark_cnvs_gain_loss_loh.bed")
    gain_idx = load_bed(f"{CNV}/ngs_benchmark_cnv_gain_cn.bed", asint=True)
    loss_idx = load_bed(f"{CNV}/ngs_benchmark_cnv_loss_cn.bed", asint=True)
    rec = json.load(open(f"{A}/phylo_cpp_wg_full_records.json"))
    from collections import Counter
    st_count = Counter(); st_by_struct = Counter()
    for r in rec:
        snv = int(r["pos"]) + 5000
        state = overlap(state_idx, r["chrom"], snv) or "neutral"
        if state == "gain":
            cn = overlap(gain_idx, r["chrom"], snv)
        elif state == "loss":
            cn = overlap(loss_idx, r["chrom"], snv)
        else:
            cn = 2 if state == "neutral" else None   # loh = copy-neutral(允許 None→顯示 LOH)
        r["cn_state"] = state
        r["cn_value"] = cn
        r["is_loh"] = (state == "loh")
        st_count[state] += 1
        st_by_struct[(state, "multi" if r["coarse_ng"] >= 2 else "single")] += 1
    json.dump(rec, open(f"{A}/phylo_cpp_wg_full_records_cn.json", "w"))
    stats = {"n": len(rec), "cn_state_dist": dict(st_count),
             "cn_state_x_structure": {f"{s}_{k}": n for (s, k), n in sorted(st_by_struct.items())},
             "cn_value_dist": dict(Counter(r["cn_value"] for r in rec if r["cn_value"] is not None)),
             "source": "SEQC2 CNV ngs_benchmark (hg38, first-hand)"}
    json.dump(stats, open(f"{A}/cn_stats.json", "w"), indent=2)
    print(json.dumps(stats, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
