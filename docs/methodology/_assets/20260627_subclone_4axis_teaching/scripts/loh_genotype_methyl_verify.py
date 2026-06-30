#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LOH genotype-vector 甲基差異直接驗證(2026-06-30):驗用戶 LOH 假設。
LOH 區只剩一條 germline allele → germline-ASM 自動抵消(SAME-HP 必然) → genotype 向量間甲基差
  = somatic-cis(突變局部足跡) + clonal-lineage(殘差),不含 germline-ASM。
直接配對(用 BAM 自算 genotype,不靠 populations 欄,避免序對不上):
  - RR_vs_1ALT: 祖先群(全R) vs 單突變群(AR型);差異位點 = 該 ALT 位置
  - sib_1ALT  : 兩個不同單突變群(AR vs RA);差異位點 = 兩 ALT 位置
near/distal 拆解(near=差異位點 ±1000bp,與 §11 ordering pilot 一致):
  - 只 near 有 Δβ = somatic-cis 局部足跡(非 lineage)
  - distal 也有 Δβ = 超出突變足跡的 clonal/lineage 訊號
per pair: dbeta_near, dbeta_distal, perm_p_distal(2000× read-label 置換)。
CN 分層(loh 主測;neutral/gain 對照看 germline-ASM 有無差別)。seed 固定。§13.0 compute。
chunk: argv[1]=total argv[2]=idx。env SM_DATA/SM_TBAM/SM_VD/SM_VCF_MODE。
"""
import json, os, sys
import numpy as np
import pysam

CHUNK_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CHUNK_IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
OUTSUF = "" if CHUNK_TOTAL == 1 else f"_part{CHUNK_IDX}"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))
TBAM = os.environ.get("SM_TBAM", "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam")
VCFD = os.environ.get("SM_VD", "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup")
VCF_MODE = os.environ.get("SM_VCF_MODE", "perchrom")
MAPQ = 20; MINREAD_CPG = 3; MIN_PER_CLASS = 6; READ_CAP = 800
NEAR_BP = 1000; MIN_BAND = 4; N_PERM = 2000
RNG = np.random.default_rng(20260630)

det = json.load(open(os.path.join(DATA, "topology_per_region.json")))["detail"][CHUNK_IDX::CHUNK_TOTAL]
_TBX = {}
def _tabix(p):
    if p not in _TBX:
        if not os.path.exists(p): _TBX[p] = None
        else:
            idx = p+".tbi" if os.path.exists(p+".tbi") else (p+".csi" if os.path.exists(p+".csi") else None)
            try: _TBX[p] = pysam.TabixFile(p, index=idx) if idx else pysam.TabixFile(p)
            except Exception: _TBX[p] = None
    return _TBX[p]

def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}.vcf.gz" if VCF_MODE == "single" else f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        tbx = _tabix(p)
        if tbx is None: continue
        try:
            for ln in tbx.fetch(chrom, max(0, s-1), e+1):
                f = ln.split("\t"); pos = int(f[1])
                if pos not in m: m[pos] = (f[3].upper(), f[4].strip().upper())
        except Exception: pass
    return m

def read_betas(a):
    """回 {ref_pos: beta}。"""
    try: mb = a.modified_bases
    except Exception: return None
    if not mb: return None
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    out = {}
    for k, lst in mb.items():
        if k[2] != "m": continue
        for qp, ml in lst:
            rr = qr.get(qp)
            if rr is not None: out[rr] = ml/255.0
    return out

def geno(a, som):
    rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence; g = []
    for pos, ref, alt in som:
        q = rq.get(pos-1)
        if q is None or seq is None: return None
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    return None if "-" in g else "".join(g)

def band_means(betas, diff_pos):
    """read 的 CpG β 依離 diff_pos 距離拆 near(<=NEAR_BP)/distal。回 (near_mean, distal_mean) 可為 None。"""
    near = []; dist = []
    for rr, b in betas.items():
        p = rr + 1
        if min(abs(p - dp) for dp in diff_pos) <= NEAR_BP: near.append(b)
        else: dist.append(b)
    return (float(np.mean(near)) if near else None, float(np.mean(dist)) if dist else None)

def perm_p_distal(d1, d2):
    """兩群 per-read distal mean β 的 |Δ| 置換 p。"""
    a = np.array([x for x in d1 if x is not None]); b = np.array([x for x in d2 if x is not None])
    if len(a) < MIN_BAND or len(b) < MIN_BAND: return None, None
    obs = abs(a.mean() - b.mean())
    pool = np.concatenate([a, b]); n1 = len(a); cnt = 0
    idx = np.arange(len(pool))
    for _ in range(N_PERM):
        RNG.shuffle(idx)
        pa = pool[idx[:n1]]; pb = pool[idx[n1:]]
        if abs(pa.mean() - pb.mean()) >= obs: cnt += 1
    return round(float(obs), 4), (cnt + 1) / (N_PERM + 1)

def main():
    tb = pysam.AlignmentFile(TBAM, "rb")
    pairs = []
    for r in det:
        chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1]); cn = r.get("cn")
        ra = ref_alt_map(chrom, s, e)
        som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
        if len(som) < 1: continue
        pos = [p for p, _, _ in som]
        # 收 read: genotype + betas
        groups = {}  # g -> list of betas dict
        n = 0
        for a in tb.fetch(chrom, s, e+1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
            bt = read_betas(a)
            if not bt or len(bt) < MINREAD_CPG: continue
            g = geno(a, som)
            if g is None: continue
            groups.setdefault(g, []).append(bt); n += 1
            if n >= READ_CAP: break
        # 分類群
        founder = "R" * len(som)
        singles = {}  # idx -> g (exactly one A at idx)
        for g in groups:
            if g.count("A") == 1 and len(groups[g]) >= MIN_PER_CLASS:
                singles[g.index("A")] = g
        big = {g: v for g, v in groups.items() if len(v) >= MIN_PER_CLASS}
        def emit(g1, g2, diff_idx, ptype):
            dp = [pos[i] for i in diff_idx]
            b1 = [band_means(bt, dp) for bt in big[g1]]; b2 = [band_means(bt, dp) for bt in big[g2]]
            n1n = [x[0] for x in b1 if x[0] is not None]; n2n = [x[0] for x in b2 if x[0] is not None]
            n1d = [x[1] for x in b1 if x[1] is not None]; n2d = [x[1] for x in b2 if x[1] is not None]
            dbn = abs(np.mean(n1n) - np.mean(n2n)) if len(n1n) >= MIN_BAND and len(n2n) >= MIN_BAND else None
            obs, pp = perm_p_distal([x[1] for x in b1], [x[1] for x in b2])
            pairs.append({"region": r["region"], "cn": cn, "ptype": ptype,
                          "n1": len(big[g1]), "n2": len(big[g2]),
                          "dbeta_near": round(float(dbn), 4) if dbn is not None else None,
                          "n_distal1": len(n1d), "n_distal2": len(n2d),
                          "dbeta_distal": obs, "perm_p_distal": pp})
        # RR vs 單突變
        if founder in big:
            for i, g in singles.items():
                emit(founder, g, [i], "RR_vs_1ALT")
        # 單突變 vs 單突變 (sibling, AR vs RA)
        sk = sorted(singles)
        for ii in range(len(sk)):
            for jj in range(ii+1, len(sk)):
                emit(singles[sk[ii]], singles[sk[jj]], [sk[ii], sk[jj]], "sib_1ALT")
    if CHUNK_TOTAL > 1:
        json.dump({"pairs": pairs}, open(os.path.join(DATA, f"loh_gt_methyl{OUTSUF}.json"), "w"), ensure_ascii=False)
        print(f"LOH-GT CHUNK {CHUNK_IDX}/{CHUNK_TOTAL}: {len(pairs)} pairs"); return
    json.dump({"pairs": pairs}, open(os.path.join(DATA, "loh_genotype_methyl_verify.json"), "w"), ensure_ascii=False)
    print("DONE", len(pairs))

if __name__ == "__main__":
    main()
