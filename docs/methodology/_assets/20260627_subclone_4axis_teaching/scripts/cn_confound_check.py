#!/usr/bin/env python3
"""驗證 sSNV 共現區是否受 CN 影響(HCC1395 SEQC2 CN 真值)。compute→file。
CN 對 VAF 的影響: gain=稀釋VAF(→假subclone) · loh/loss=抬高VAF(→假clonal,隱藏subclone) · neutral=乾淨(VAF≈CCF)。
"""
import json, glob, os
from collections import Counter, defaultdict

PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
MIN_COV, MIN_ALT = 6, 2
out = []
out.append(f"CNBED 存在: {os.path.exists(CNBED)}  ({CNBED})")

groups = []
for f in sorted(glob.glob(f"{PILOT}/mlhp_part_*.json")):
    for g in json.load(open(f))["groups"]:
        if g.get("n_sSNV", 0) >= 2:
            groups.append(g)
out.append(f"HCC1395 sSNV 共現區(≥2 sSNV) 總數: {len(groups)}")

# region 層 c + CN
def region_c(g):
    pbh = g.get("populations_by_hp", {}) or {}
    alt = set()
    for fam in ("1", "2", "3"):
        for gt in (pbh.get(fam, {}) or {}):
            if "A" in gt: alt.add((fam, gt))
    return len(alt)

# ① 全共現區 CN 分布
cn_all = Counter(g.get("cn") for g in groups)
tot = len(groups)
out.append("\n=== ① 全 sSNV 共現區的 CN 狀態分布 ===")
for k in ["neutral", "loh", "gain", "loss", None]:
    if cn_all.get(k): out.append(f"  {str(k):8}: {cn_all[k]:5d} ({100*cn_all[k]/tot:.1f}%)")
clean = cn_all.get("neutral", 0)
out.append(f"  → CN-neutral(VAF 可乾淨解讀)= {clean} ({100*clean/tot:.1f}%)")
out.append(f"  → CN-altered(gain/loh/loss, VAF 受混淆)= {tot-clean} ({100*(tot-clean)/tot:.1f}%)")
out.append(f"  🔴 其中 gain(稀釋VAF→可能假subclone)= {cn_all.get('gain',0)} ({100*cn_all.get('gain',0)/tot:.1f}%)")

# ② CN × c(clone 數)
out.append("\n=== ② CN 狀態 × c(region 層 clone 數) ===")
cn_by_c = defaultdict(Counter)
for g in groups:
    cn_by_c[min(region_c(g), 3)][g.get("cn")] += 1  # c: 0,1,2,3+
out.append(f"{'c':>4} | {'neutral':>8}{'loh':>7}{'gain':>7}{'loss':>6} | neutral%")
for c in [0, 1, 2, 3]:
    cc = cn_by_c[c]; n = sum(cc.values())
    lab = "3+" if c == 3 else str(c)
    out.append(f"{lab:>4} | {cc.get('neutral',0):>8}{cc.get('loh',0):>7}{cc.get('gain',0):>7}{cc.get('loss',0):>6} | {100*cc.get('neutral',0)/n:.0f}%" if n else f"{lab:>4} | (0)")

# ③ 巢狀 clone→subclone 區的 CN(關鍵:0.95→0.31 pattern 是否 CN 假象)
out.append("\n=== ③ 巢狀 clone→subclone 區(2-ALT 巢狀) 的 CN 狀態 ===")
out.append("   祖先VAF~1.0→衍生VAF~0.4 的 pattern:neutral=真subclone候選 · gain=可能稀釋假象")
nested_cn = Counter(); nested_gain_vaf = []; nested_neutral_vaf = []
for g in groups:
    positions = g.get("positions", [])
    pbh = g.get("populations_by_hp", {}) or {}
    cbh = g.get("col_coverage_by_hp", {}) or {}
    cn = g.get("cn")
    for fam in ("1", "2"):
        pops = pbh.get(fam, {}) or {}; cc = cbh.get(fam, {}) or {}
        alts = [gt for gt in pops if "A" in gt]
        if len(alts) != 2: continue
        a, b = alts
        sa = {i for i, ch in enumerate(a) if ch == "A"}; sb = {i for i, ch in enumerate(b) if ch == "A"}
        if not (sa < sb or sb < sa): continue
        parent, child = (a, b) if sa < sb else (b, a)
        der_idx = [i for i, ch in enumerate(child) if ch == "A" and parent[i] != "A"]
        dvs = []
        for i in der_idx:
            if i < len(positions):
                p = str(positions[i])
                if p in cc:
                    nr, na = cc[p]
                    if nr + na >= MIN_COV: dvs.append(na / (nr + na))
        dv = sum(dvs) / len(dvs) if dvs else None
        nested_cn[cn] += 1
        if dv is not None:
            if cn == "gain": nested_gain_vaf.append(dv)
            elif cn == "neutral": nested_neutral_vaf.append(dv)
ntot = sum(nested_cn.values())
for k in ["neutral", "loh", "gain", "loss"]:
    if nested_cn.get(k): out.append(f"  {k:8}: {nested_cn[k]:5d} ({100*nested_cn[k]/ntot:.1f}%)")
out.append(f"  → 巢狀區總 {ntot};CN-neutral(乾淨真subclone候選)= {nested_cn.get('neutral',0)} ({100*nested_cn.get('neutral',0)/ntot:.1f}%)")
out.append(f"  🔴 巢狀區在 gain(衍生VAF可能=拷貝稀釋非subclone)= {nested_cn.get('gain',0)} ({100*nested_cn.get('gain',0)/ntot:.1f}%)")
def med(xs): return round(sorted(xs)[len(xs)//2], 3) if xs else None
out.append(f"  衍生VAF 中位: neutral 區={med(nested_neutral_vaf)} (n={len(nested_neutral_vaf)}) vs gain 區={med(nested_gain_vaf)} (n={len(nested_gain_vaf)})")
out.append(f"     若兩者接近 → gain 的低VAF 未必來自稀釋;若 gain 明顯更低 → 稀釋跡象")

txt = "\n".join(out)
open("/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260710_cn_confound_out.txt", "w").write(txt)
print(txt)
