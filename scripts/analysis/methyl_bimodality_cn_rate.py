#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CN-aware re-run of the within-genotype-group methylation bimodality test.

Reuses the parallel-session methyl_auxiliary GMM logic (2-comp BIC + peak diff>=0.2 +
subcomp>=4) but (a) takes cn from the SAVANA SM_CNBED (not topology's cn=unknown) and
(b) records EVERY tested group's (cn, is_bimodal) -> bimodality RATE per CN state.
This is the causal proof for methylation<->CN-gain (Q3): does CN-gain RAISE the
bimodality rate vs neutral? Plus the Q5 auxiliary flag: in hard-to-build tree regions
(ambiguous / underdetermined / incompatible), an ALT-cluster with cn-clean (neutral)
bimodality = L3 hidden-subclone / branching candidate (methylation-assisted tree hint).

Env: SM_TBAM, SM_VD (VCF dir), SM_TOPO (topology_per_region.json), SM_CNBED, SM_OUT.
Read-only on BAM. Output <SM_OUT>.
"""
import json
import os
import sys
from collections import defaultdict, Counter

import numpy as np
import pysam
from sklearn.mixture import GaussianMixture

TBAM = os.environ["SM_TBAM"]
VD = os.environ["SM_VD"]
TOPO = os.environ["SM_TOPO"]
CNBED = os.environ["SM_CNBED"]
OUT = os.environ["SM_OUT"]
MAPQ, MINREAD, MINGRP, SEP, MINSUB, READ_CAP = 20, 3, 12, 0.2, 4, 500
HARD = ("ambiguous", "underdetermined", "incompatible", "recurrence")
DET_CATS = ("A_determined", "A_ambiguous", "C_underdetermined", "incompatible",
            "recurrence_required", "B_pairwise", "E_subcube")


def norm_det(d):
    d = d or "other"
    for h in DET_CATS:
        if d.startswith(h):
            return h
    return "other"

_TBX = {}
def tabix(p):
    if p not in _TBX:
        if not os.path.exists(p):
            _TBX[p] = None
        else:
            idx = p + ".tbi" if os.path.exists(p + ".tbi") else (
                p + ".csi" if os.path.exists(p + ".csi") else None)
            try:
                _TBX[p] = pysam.TabixFile(p, index=idx) if idx else pysam.TabixFile(p)
            except Exception:
                _TBX[p] = None
    return _TBX[p]

def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        tbx = tabix(f"{VD}/filtered_snv_{src}.vcf.gz")
        if tbx is None:
            continue
        try:
            for ln in tbx.fetch(chrom, max(0, s - 1), e + 1):
                f = ln.split("\t"); pos = int(f[1])
                if pos not in m:
                    m[pos] = (f[3].upper(), f[4].strip().upper())
        except Exception:
            pass
    return m

def load_cn(bed):
    iv = defaultdict(list)
    for ln in open(bed):
        p = ln.split()
        if len(p) >= 4 and p[1].isdigit():
            iv[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in iv:
        iv[c].sort()
    return iv

def cn_at(iv, chrom, pos):
    for s, e, lab in iv.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"

def read_meth(a):
    try:
        mb = a.modified_bases
    except Exception:
        return None
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m":
                continue
            for qpos, mlq in lst:
                rr = qr.get(qpos)
                if rr is not None:
                    meth[rr] = mlq / 255.0
    return meth

def geno_str(a, som):
    rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence; g = []
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None:
            g.append("-"); continue
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    return "".join(g)

def hp_norm(a):
    if not a.has_tag("HP"):
        return 0
    try:
        return int(str(a.get_tag("HP")).split("-")[0])
    except Exception:
        return 0

def bimodal_test(vals):
    x = np.array(vals).reshape(-1, 1)
    if len(x) < MINGRP:
        return False, None
    try:
        g1 = GaussianMixture(1, covariance_type="full", random_state=0).fit(x)
        g2 = GaussianMixture(2, covariance_type="full", random_state=0, n_init=1).fit(x)
    except Exception:
        return False, None
    if g2.bic(x) >= g1.bic(x):
        return False, None
    m = g2.means_.flatten()
    if abs(m[0] - m[1]) < SEP:
        return False, None
    lab = g2.predict(x)
    c0, c1 = int((lab == 0).sum()), int((lab == 1).sum())
    if min(c0, c1) < MINSUB:
        return False, None
    return True, {"means": [round(float(m[0]), 3), round(float(m[1]), 3)], "sizes": [c0, c1], "labels": lab.tolist()}

def main():
    ct = int(os.environ.get("SM_CHUNK_TOTAL", "1"))
    ci = int(os.environ.get("SM_CHUNK_IDX", "0"))
    outp = OUT if ct == 1 else OUT.replace(".json", f"_part{ci}.json")
    det = json.load(open(TOPO))["detail"][ci::ct]
    iv = load_cn(CNBED)
    tb = pysam.AlignmentFile(TBAM, "rb")
    rate = defaultdict(lambda: [0, 0])   # cn -> [n_bimodal, n_tested]
    flag_by_cn = defaultdict(Counter)    # cn -> {hp/cn_gain/residual}
    q5 = []                              # Q5 hard-region ALT-cluster candidates
    det_ct = defaultdict(Counter)        # determinacy -> {total/none/hp/cn_gain/residual/q5}
    n_reg = 0
    for r in det:
        region = r["region"]; chrom = r["chrom"]; s = r["start"]
        e = int(region.split("-")[-1])
        ra = ref_alt_map(chrom, s, e)
        som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
        if len(som) < 2:
            continue
        n_reg += 1
        mycn = cn_at(iv, chrom, (s + e) // 2)
        hard = any(h in r.get("determinacy", "") for h in HARD)
        reg_status = "none"; reg_q5 = False
        _rank = {"none": 0, "hp": 1, "cn_gain": 2, "residual": 3}
        meth = {}; geno = {}; hp = {}; n = 0
        for a in tb.fetch(chrom, s, e + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ:
                continue
            m = read_meth(a)
            if m is None:
                continue
            g = geno_str(a, som)
            if "-" in g:
                continue
            meth[a.query_name] = m; geno[a.query_name] = g; hp[a.query_name] = hp_norm(a)
            n += 1
            if n >= READ_CAP:
                break
        grp = defaultdict(list)
        for rn, g in geno.items():
            grp[g].append(rn)
        for g, ids in grp.items():
            if len(ids) < MINGRP:
                continue
            mb = {rn: float(np.mean(list(meth[rn].values()))) for rn in ids
                  if len(meth[rn]) >= MINREAD}
            if len(mb) < MINGRP:
                continue
            rns = list(mb.keys()); vals = [mb[rn] for rn in rns]
            is_bi, info = bimodal_test(vals)
            rate[mycn][1] += 1
            if not is_bi:
                continue
            rate[mycn][0] += 1
            lab = info["labels"]
            def domhp(sel):
                c = Counter(hp[rns[i]] for i in range(len(rns)) if lab[i] == sel)
                tot = sum(c.values()) or 1
                for k in (1, 2, 3):
                    if c.get(k, 0) / tot >= 0.7:
                        return k
                return 0
            d0, d1 = domhp(0), domhp(1)
            hp_exp = (d0 in (1, 2) and d1 in (1, 2) and d0 != d1)
            cls = "hp" if hp_exp else ("cn_gain" if mycn == "gain" else "residual")
            flag_by_cn[mycn][cls] += 1
            if _rank[cls] > _rank[reg_status]:
                reg_status = cls
            # Q5: hard region + ALT cluster + cn-clean(neutral) + not HP -> L3 subclone/branch candidate
            if hard and (not hp_exp) and mycn == "neutral" and "A" in g:
                reg_q5 = True
                q5.append({"region": region, "determinacy": r.get("determinacy"),
                           "genotype": g, "n_reads": len(ids), "cn": mycn,
                           "means": info["means"], "sizes": info["sizes"],
                           "delta_beta": round(abs(info["means"][0] - info["means"][1]), 3),
                           "n_sSNV": r.get("n_sSNV"), "topology_type": r.get("topology_type")})
        dcat = norm_det(r.get("determinacy"))
        det_ct[dcat]["total"] += 1
        det_ct[dcat][reg_status] += 1
        if reg_q5:
            det_ct[dcat]["q5"] += 1
    tb.close()
    rate_out = {cn: {"n_bimodal": v[0], "n_tested": v[1],
                     "bimodal_rate": round(v[0] / v[1], 4) if v[1] else None}
                for cn, v in sorted(rate.items())}
    summary = {"n_regions_multiSNV": n_reg,
               "bimodal_rate_by_cn": rate_out,
               "flag_by_cn": {k: dict(v) for k, v in flag_by_cn.items()},
               "n_q5_hardregion_altcluster_candidates": len(q5),
               "params": {"MINGRP": MINGRP, "SEP": SEP, "MINSUB": MINSUB, "HARD": HARD}}
    # raw counts for chunk-merge
    raw = {cn: [v[0], v[1]] for cn, v in rate.items()}
    json.dump({"summary": summary, "raw_rate": raw,
               "flag_by_cn": {k: dict(v) for k, v in flag_by_cn.items()},
               "det_ct": {k: dict(v) for k, v in det_ct.items()},
               "q5_candidates": q5}, open(outp, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(rate_out, ensure_ascii=False))
    print(f"Q5 candidates: {len(q5)}  -> {outp}")


if __name__ == "__main__":
    main()
