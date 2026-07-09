#!/usr/bin/env python3
"""VAF 峰值分群 + 2-ALT-群父子關係 pilot（compute→file，供讀回）。
c(Def B) = 除 REF 外 distinct ALT 組合群數（=該區 somatic clone 數）。
"""
import json, glob
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
MIN_COV, MIN_ALT = 6, 2

def wd(s): return PILOT if s == "HCC1395" else f"{MSROOT}/{s}"

def iter_groups(s):
    for f in sorted(glob.glob(f"{wd(s)}/mlhp_part_*.json")):
        for g in json.load(open(f))["groups"]:
            if g.get("n_sSNV", 0) >= 2:
                yield g

def peaks_from_hist(vafs, binw=0.05):
    """簡單直方圖 + local-maxima 峰值偵測（免 scipy 依賴）。"""
    import math
    if not vafs: return [], []
    nb = int(round(1.0 / binw))
    h = [0] * nb
    for v in vafs:
        b = min(nb - 1, int(v / binw))
        h[b] += 1
    centers = [(i + 0.5) * binw for i in range(nb)]
    # 平滑（3-bin moving avg）後找 local maxima，峰高 ≥ 總數 3%
    sm = [ (h[max(0,i-1)]+h[i]+h[min(nb-1,i+1)])/3 for i in range(nb) ]
    tot = sum(h); thr = 0.03 * tot
    pk = []
    for i in range(nb):
        left = sm[i-1] if i>0 else -1
        right = sm[i+1] if i<nb-1 else -1
        if sm[i] >= left and sm[i] >= right and h[i] >= thr:
            pk.append((round(centers[i],3), h[i]))
    return pk, h

out = []
out.append("="*100)
out.append("§A  c(Def B) = 除 REF 外 distinct ALT 組合群數 = 該區 somatic clone 數")
out.append("   c=0 germline-only(無 full-cov somatic 組合) · c=1 單 somatic clone · c≥2 有 subclonal 結構")
out.append("   同時給 min-2-read 穩健版(丟單-read genotype=可能 error;raw=最高可能上界)")
out.append("="*100)
out.append(f"{'樣本':16}{'c=0':>7}{'c=1':>7}{'c=2':>7}{'c=3':>7}{'c=4':>7}{'c≥5':>7}{'maxC':>6}  || min2穩健: c≥2佔raw c≥2 之 %")
grand_raw = Counter(); grand_r2 = Counter()
for s in SAMPLES:
    raw = Counter(); r2 = Counter()
    for g in iter_groups(s):
        pbh = g.get("populations_by_hp", {}) or {}
        for fam in ("1", "2"):
            pops = pbh.get(fam, {}) or {}
            if not pops: continue
            alt_raw = {gt for gt in pops if "A" in gt}
            alt_r2 = {gt for gt in pops if "A" in gt and pops[gt] >= 2}
            raw[len(alt_raw)] += 1; r2[len(alt_r2)] += 1
            grand_raw[len(alt_raw)] += 1; grand_r2[len(alt_r2)] += 1
    def gg(c, k): return c.get(k, 0)
    cge5 = sum(v for k, v in raw.items() if k >= 5)
    maxc = max(raw) if raw else 0
    raw_ge2 = sum(v for k, v in raw.items() if k >= 2)
    r2_ge2 = sum(v for k, v in r2.items() if k >= 2)
    pct = f"{100*r2_ge2/raw_ge2:.0f}%" if raw_ge2 else "-"
    out.append(f"{s:16}{gg(raw,0):>7}{gg(raw,1):>7}{gg(raw,2):>7}{gg(raw,3):>7}{gg(raw,4):>7}{cge5:>7}{maxc:>6}  || {r2_ge2}/{raw_ge2} = {pct}")
out.append(f"{'全樣本':16}{grand_raw.get(0,0):>7}{grand_raw.get(1,0):>7}{grand_raw.get(2,0):>7}{grand_raw.get(3,0):>7}{grand_raw.get(4,0):>7}{sum(v for k,v in grand_raw.items() if k>=5):>7}")
tot = sum(grand_raw.values())
out.append(f"   全樣本: c≥1(≥1 somatic clone) {100*sum(v for k,v in grand_raw.items() if k>=1)/tot:.0f}% · c≥2(有 subclonal 結構) {100*sum(v for k,v in grand_raw.items() if k>=2)/tot:.0f}%")

out.append("")
out.append("="*100)
out.append("§B  genome-wide somatic VAF 峰值(pooled;每 somatic 位點 within-family VAF,nALT≥2 & cov≥6)")
out.append("   峰值 = 可能的 clone/subclone 頻率群(SciClone/PyClone 式);峰在 VAF≈0.5/1.0=clonal,較低=subclonal")
out.append("   🔴 raw VAF 有 CN confound(gain 稀釋/LOH 抬高);未做 CN 校正=僅 HCC1395 有 SEQC2 可轉 CCF")
out.append("="*100)
for s in SAMPLES:
    vafs = []
    for g in iter_groups(s):
        cbh = g.get("col_coverage_by_hp", {}) or {}
        for fam in ("1", "2"):
            for pos, (nr, na) in (cbh.get(fam, {}) or {}).items():
                if na >= MIN_ALT and nr + na >= MIN_COV:
                    vafs.append(na / (nr + na))
    pk, _ = peaks_from_hist(vafs)
    pkstr = " · ".join(f"VAF~{c}({n})" for c, n in pk)
    out.append(f"{s:16} n_somatic_pos={len(vafs):>6}  峰值: {pkstr}")

out.append("")
out.append("="*100)
out.append("§C  2-ALT-群(c=2)區的父子關係:巢狀(clone→subclone) vs 姊妹(分支);祖先 VAF vs 衍生 VAF")
out.append("   巢狀 = 一 ALT 組合 ⊂ 另一(如 AR⊂AA);祖先位點(共有 A)VAF 應 > 衍生位點(獨有 A)VAF")
out.append("   驗證你的假設: 多個區是否一致呈『高 VAF 祖先 → 低 VAF 衍生』(=clone→subclone 巢狀)")
out.append("="*100)
def is_subset(a, b):
    """a 的 ALT 位點 ⊂ b 的 ALT 位點(genotype 同長字串)。"""
    if len(a) != len(b): return False
    sa = {i for i, ch in enumerate(a) if ch == "A"}
    sb = {i for i, ch in enumerate(b) if ch == "A"}
    return sa < sb  # 真子集
for s in SAMPLES:
    nested = 0; sister = 0
    anc_vafs = []; der_vafs = []; examples = []
    pairs = []  # (anc_vaf, der_vaf)
    for g in iter_groups(s):
        positions = g.get("positions", [])
        pbh = g.get("populations_by_hp", {}) or {}
        cbh = g.get("col_coverage_by_hp", {}) or {}
        for fam in ("1", "2"):
            pops = pbh.get(fam, {}) or {}
            cc = cbh.get(fam, {}) or {}
            alts = [gt for gt in pops if "A" in gt]
            if len(alts) != 2: continue
            a, b = alts
            # 判巢狀 or 姊妹
            if is_subset(a, b) or is_subset(b, a):
                nested += 1
                parent, child = (a, b) if is_subset(a, b) else (b, a)
                # 祖先位點 = parent 的 A;衍生位點 = child 有 parent 沒有的 A
                anc_idx = [i for i, ch in enumerate(parent) if ch == "A"]
                der_idx = [i for i, ch in enumerate(child) if ch == "A" and parent[i] != "A"]
                def vaf_at(idx):
                    vs = []
                    for i in idx:
                        if i < len(positions):
                            p = str(positions[i])
                            if p in cc:
                                nr, na = cc[p]
                                if nr + na >= MIN_COV: vs.append(na / (nr + na))
                    return sum(vs) / len(vs) if vs else None
                av = vaf_at(anc_idx); dv = vaf_at(der_idx)
                if av is not None and dv is not None:
                    anc_vafs.append(av); der_vafs.append(dv); pairs.append((av, dv))
                    if s == "HCC1395" and len(examples) < 6:
                        examples.append((g["chrom"], g["start"], round(av, 2), round(dv, 2)))
            else:
                sister += 1
    def med(xs): return sorted(xs)[len(xs)//2] if xs else float("nan")
    n_pair = len(pairs)
    anc_gt_der = sum(1 for av, dv in pairs if av > dv)
    out.append(f"{s:16} 巢狀(clone→subclone)={nested:>5} 姊妹(分支)={sister:>5}  | 有VAF對={n_pair:>4} 祖先VAF>衍生VAF={anc_gt_der}({100*anc_gt_der/n_pair:.0f}% if n_pair else 0)  祖先中位VAF={med(anc_vafs):.2f} 衍生中位VAF={med(der_vafs):.2f}")
    if s == "HCC1395" and examples:
        out.append("     HCC1395 範例(祖先VAF→衍生VAF): " + " · ".join(f"{c}:{st}({a}→{d})" for c, st, a, d in examples))

txt = "\n".join(str(x) for x in out)
open("/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260709_vaf_clone_pilot_out.txt", "w").write(txt)
print(txt)
