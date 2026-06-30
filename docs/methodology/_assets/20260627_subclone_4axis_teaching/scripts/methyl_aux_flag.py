#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Methyl auxiliary differential-clustering FLAG(2026-06-30):
給『sSNV 定義的群』(germline=all-R vs mutated=≥1 ALT) → 檢定甲基是否提供差異分群(輔助驗證/標記用途)。
每區:兩群 read 的 mean β → |Δβ| 效應量 + N_PERM 次標籤置換經驗 p(對小 n 穩健,不依賴分布假設)。
merge:BH-FDR 跨區 → flag = (FDR<0.05 AND |Δβ|>=0.1);join founder balanced_acc(LOOCV 驗證)。
回答『甲基能否當輔助標記/驗證某群差異分群』:量化多少區甲基『真的、顯著地』提供差異分群訊號(非隨機)。
🔴 非循環但 confound 必標:群由 sSNV 定義(非甲基)→ 此檢定問『甲基是否 co-segregate genotype 軸』
   = cis-ASM 一致性(甲基足跡與突變共分離),非獨立 lineage 證據。是 consistency-flag,不是 subclone 確認器。
seed 固定(reproducible §7.2)。chunk: argv[1]=total argv[2]=idx。env SM_DATA/SM_TBAM/SM_VD/SM_VCF_MODE。§13.0 compute。
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
MAPQ = 20; MINREAD_CPG = 3; MIN_PER_CLASS = 6; READ_CAP = 600; N_PERM = 2000
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

def read_meth(a):
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

def perm_p(gm, lab, obs):
    """標籤置換經驗 p:H0 = 兩群 mean β 無差。回 (count(perm|Δ|>=obs)+1)/(N+1)。"""
    n = len(lab); cnt = 0
    idx = np.arange(n)
    for _ in range(N_PERM):
        RNG.shuffle(idx)
        pl = lab[idx]
        m0 = gm[pl == 0]; m1 = gm[pl == 1]
        if len(m0) == 0 or len(m1) == 0: continue
        if abs(m1.mean() - m0.mean()) >= obs: cnt += 1
    return (cnt + 1) / (N_PERM + 1)

def main():
    tb = pysam.AlignmentFile(TBAM, "rb")
    rows = []
    for r in det:
        chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
        ra = ref_alt_map(chrom, s, e)
        som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
        if len(som) < 1: continue
        gm = []; lab = []; n = 0
        for a in tb.fetch(chrom, s, e+1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
            m = read_meth(a)
            if not m: continue
            g = geno(a, som)
            if g is None: continue
            bv = list(m.values())
            if len(bv) < MINREAD_CPG: continue
            gm.append(float(np.mean(bv))); lab.append(0 if "A" not in g else 1); n += 1
            if n >= READ_CAP: break
        gm = np.array(gm); lab = np.array(lab)
        n0 = int((lab == 0).sum()); n1 = int((lab == 1).sum())
        if n0 < MIN_PER_CLASS or n1 < MIN_PER_CLASS: continue
        dbeta = abs(gm[lab == 0].mean() - gm[lab == 1].mean())
        pp = perm_p(gm, lab, dbeta)
        rows.append({"region": r["region"], "cn": r.get("cn"), "n": len(lab), "n_germ": n0, "n_mut": n1,
                     "dbeta_centroid": round(float(dbeta), 4), "perm_p": round(float(pp), 5)})
    if CHUNK_TOTAL > 1:
        json.dump({"rows": rows}, open(os.path.join(DATA, f"methyl_aux_flag{OUTSUF}.json"), "w"), ensure_ascii=False)
        print(f"AUX-FLAG CHUNK {CHUNK_IDX}/{CHUNK_TOTAL}: {len(rows)} regions"); return
    json.dump({"rows": rows}, open(os.path.join(DATA, "methyl_aux_flag.json"), "w"), ensure_ascii=False)
    print("DONE", len(rows))

if __name__ == "__main__":
    main()
