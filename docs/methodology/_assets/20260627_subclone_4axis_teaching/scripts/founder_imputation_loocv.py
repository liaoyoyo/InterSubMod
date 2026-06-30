#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Founder-imputation LOOCV(2026-06-30):驗證『甲基能否把 sSNV-uncovered read 補進遺傳定義的群』。
半監督非循環:群由 sSNV 定義(germline=all-R vs mutated=≥1 ALT),甲基當分類器,LOOCV 驗準確度。
每區:覆蓋全 sSNV 的 read 取 genotype label + 區內 mean β → LOOCV nearest-centroid(by mean β)
  → 報 accuracy / balanced_accuracy / majority_baseline / Δβ(群心距,=ASM強度)。
扣 majority baseline 才是真預測力(防 class imbalance 假高)。§13.0 compute。
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
MAPQ = 20; MINREAD_CPG = 3; MIN_PER_CLASS = 6; READ_CAP = 600

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

def main():
    tb = pysam.AlignmentFile(TBAM, "rb")
    rows = []
    for r in det:
        chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
        ra = ref_alt_map(chrom, s, e)
        som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
        if len(som) < 1: continue
        gm = []; lab = []  # mean β, label(0=germline all-R,1=mutated)
        n = 0
        for a in tb.fetch(chrom, s, e+1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
            m = read_meth(a)
            if not m: continue
            g = geno(a, som)
            if g is None: continue
            bv = list(m.values())
            if len(bv) < MINREAD_CPG: continue
            gm.append(float(np.mean(bv))); lab.append(0 if "A" not in g else 1)
            n += 1
            if n >= READ_CAP: break
        gm = np.array(gm); lab = np.array(lab)
        n0 = int((lab == 0).sum()); n1 = int((lab == 1).sum())
        if n0 < MIN_PER_CLASS or n1 < MIN_PER_CLASS: continue
        # LOOCV nearest-centroid(mean β)
        correct = 0; tp0 = tp1 = 0
        for i in range(len(lab)):
            c0 = gm[(lab == 0)]; c1 = gm[(lab == 1)]
            if i < len(lab) and lab[i] == 0: c0 = np.delete(c0, np.where(c0 == gm[i])[0][:1])
            else: c1 = np.delete(c1, np.where(c1 == gm[i])[0][:1])
            if len(c0) == 0 or len(c1) == 0: continue
            pred = 0 if abs(gm[i]-c0.mean()) <= abs(gm[i]-c1.mean()) else 1
            if pred == lab[i]:
                correct += 1
                if lab[i] == 0: tp0 += 1
                else: tp1 += 1
        acc = correct/len(lab)
        bal = 0.5*(tp0/max(n0,1) + tp1/max(n1,1))
        maj = max(n0, n1)/len(lab)
        dbeta = abs(gm[lab==0].mean() - gm[lab==1].mean())
        rows.append({"region": r["region"], "cn": r.get("cn"), "n": len(lab), "n_germ": n0, "n_mut": n1,
                     "accuracy": round(acc,3), "balanced_acc": round(bal,3), "majority_baseline": round(maj,3),
                     "acc_above_majority": round(acc-maj,3), "dbeta_centroid": round(float(dbeta),3)})
    if CHUNK_TOTAL > 1:
        json.dump({"rows": rows}, open(os.path.join(DATA, f"founder_imp{OUTSUF}.json"), "w"), ensure_ascii=False)
        print(f"FOUNDER-IMP CHUNK {CHUNK_IDX}/{CHUNK_TOTAL}: {len(rows)} regions"); return
    json.dump({"rows": rows}, open(os.path.join(DATA, "founder_imputation_loocv.json"), "w"), ensure_ascii=False)
    print("DONE", len(rows))

if __name__ == "__main__":
    main()
