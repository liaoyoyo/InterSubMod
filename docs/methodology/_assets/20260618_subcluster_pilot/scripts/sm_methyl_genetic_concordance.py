"""[甲基-遺傳一致性校正] 回答「能否用多-sSNV 的遺傳真值校正單-sSNV 的甲基子結構顯著性門檻」。
設計: 多-sSNV 區 = 有遺傳真值(sSNV 連鎖定義的 sub-lineage)。限同一 HP(控單倍型)→ 區內依 sSNV-combo 的
遺傳 sub-partition = subclone 真值。測「甲基能否 recover 此遺傳 partition」(PERMANOVA R² + label-shuffle null)。
若甲基顯著 recover 已知遺傳 subclone → 存在可轉移門檻給單 sSNV; 若不能(double-dip) → 不能外推。
argv: chroms_or_ALL out_path。
"""
import sys
import json
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import mannwhitneyu  # noqa (kept for parity)
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RD = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report/data"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
MAX_READS = 160
MIN_SHARED = 5     # 兩 read 至少共享 5 CpG 才算距離
N_PERM = 199


def read_meth_geno_hp(a, som):
    try:
        mb = a.modified_bases
    except Exception:
        return None
    pairs = a.get_aligned_pairs(matches_only=True)
    qr = {q: r for q, r in pairs}; rq = {r: q for q, r in pairs}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != 'm':
                continue
            for qpos, mlq in lst:
                r = qr.get(qpos)
                if r is not None:
                    meth[r] = mlq / 255.0
    seq = a.query_sequence
    if seq is None:
        return None
    g = []
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None:
            g.append("-"); continue
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    hp = a.get_tag("HP") if a.has_tag("HP") else None
    hpc = "HP1" if hp in (1, 11) else ("HP2" if hp in (2, 21) else "other")
    return meth, "".join(g), hpc


def permanova(D, labels, rng):
    """one-way distance-based PERMANOVA → (R2, pseudo_F, p)。D=方陣, labels=list。"""
    n = len(labels)
    labs = list(set(labels))
    g = len(labs)
    if g < 2 or n <= g:
        return None
    idx = {l: [i for i in range(n) if labels[i] == l] for l in labs}
    SST = (D ** 2).sum() / (2 * n)

    def within(lab_arr):
        ss = 0.0
        for l in labs:
            ii = [i for i in range(n) if lab_arr[i] == l]
            ng = len(ii)
            if ng < 2:
                continue
            sub = D[np.ix_(ii, ii)]
            ss += (sub ** 2).sum() / (2 * ng)
        return ss

    SSW = within(labels)
    SSA = SST - SSW
    if SSW <= 0 or SST <= 0:
        return None
    F = (SSA / (g - 1)) / (SSW / (n - g))
    R2 = SSA / SST
    # permutation null
    ge = 1
    larr = np.array(labels)
    for _ in range(N_PERM):
        perm = rng.permutation(larr)
        ssw = within(list(perm))
        ssa = SST - ssw
        fp = (ssa / (g - 1)) / (ssw / (n - g)) if ssw > 0 else 0
        if fp >= F:
            ge += 1
    p = ge / (N_PERM + 1)
    return round(float(R2), 4), round(float(F), 3), round(float(p), 4)


def main(chroms, out_path):
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    chset = set(chroms)
    tgt = [r for r in reg if (chroms == ["ALL"] or r["chrom"] in chset)
           and r["cn"] in ("loh", "neutral")
           and r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")
           and (r.get("n_populations") or 0) >= 2]
    rng = np.random.default_rng(0)
    n_tgt = len(tgt); n_testable = n_sig = 0
    records = []
    for idx, r in enumerate(tgt):
        chrom = r["chrom"]
        som = [(p, cen[f"{chrom}:{p}"]["ref"], cen[f"{chrom}:{p}"]["alt"])
               for p in range(r["start"], r["end"] + 1)
               if f"{chrom}:{p}" in cen and cen[f"{chrom}:{p}"].get("somatic")]
        if len(som) < 2:
            continue
        rec = {"region": r["region"], "shape": r["tree_shape"], "cn": r["cn"],
               "span": r["end"] - r["start"], "n_reads": 0, "n_geno_grp": 0, "R2": "", "pF": "", "p": "", "status": ""}
        if (r["end"] - r["start"]) > 80000:   # 無 read 能跨完所有 som（max read ~76kb）→ 單分子測不了
            rec["status"] = "region_too_wide(>80kb)"; records.append(rec); continue
        reads = {}
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_geno_hp(a, som)
            if res is None or "-" in res[1]:
                continue
            reads[a.query_name] = res
            if len(reads) >= MAX_READS:
                break
        # genetic partition = geno-string 群 (LOH 區: = subclone 真值, germline ASM 已被 LOH 排除)
        grp = defaultdict(list)
        for rn, v in reads.items():
            grp[v[1]].append(rn)
        grps = {g: ids for g, ids in grp.items() if len(ids) >= 3}
        rec.update({"n_reads": len(reads), "n_geno_grp": len(grps)})
        if len(grps) < 2:
            rec["status"] = "single_geno_grp"; records.append(rec); continue
        sub = reads
        # build read list + labels + methylation distance (CpG 限 region 窗內→密集共享)
        lo, hi = r["start"] - 1000, r["end"] + 1000
        wm = {rn: {c: b for c, b in sub[rn][0].items() if lo <= c <= hi} for rn in
              [x for ids in grps.values() for x in ids]}
        names = []; labels = []
        for g, ids in grps.items():
            for rn in ids:
                if len(wm[rn]) >= MIN_SHARED:
                    names.append(rn); labels.append(g)
        # 重判群數（窗內過濾後）
        from collections import Counter as _C
        if len(set(labels)) < 2 or min(_C(labels).values()) < 3:
            rec["status"] = "sparse_in_window"; records.append(rec); continue
        n = len(names)
        D = np.zeros((n, n))
        for i in range(n):
            mi = wm[names[i]]
            for j in range(i + 1, n):
                mj = wm[names[j]]
                shared = set(mi) & set(mj)
                if len(shared) < MIN_SHARED:
                    D[i, j] = D[j, i] = np.nan
                else:
                    d = np.mean([abs(mi[c] - mj[c]) for c in shared])
                    D[i, j] = D[j, i] = d
        if np.isnan(D).any():
            # drop reads with too many nan (keep complete submatrix greedily): simple → skip if >20% nan
            frac = np.isnan(D).sum() / (n * n - n)
            if frac > 0.2:
                rec["status"] = "sparse_shared_cpg"; records.append(rec); continue
            D = np.nan_to_num(D, nan=np.nanmean(D))
        n_testable += 1
        res = permanova(D, labels, rng)
        if res is None:
            rec["status"] = "permanova_NA"; records.append(rec); continue
        R2, F, p = res
        rec.update({"R2": R2, "pF": F, "p": p, "status": "tested"})
        if p < 0.05:
            n_sig += 1
        records.append(rec)
        if idx % 100 == 0:
            print(f"{idx}/{n_tgt} testable={n_testable} sig={n_sig}", flush=True)
    tb.close()
    tested = [r for r in records if r["status"] == "tested"]
    r2s = [r["R2"] for r in tested]
    sig = [r for r in tested if r["p"] < 0.05]

    def cn_block(cn):
        st = [r for r in tested if r["cn"] == cn]
        sg = [r for r in st if r["p"] < 0.05]
        return {"tested": len(st), "recover_sig": len(sg),
                "rate": round(len(sg) / len(st), 4) if st else None,
                "R2_median": round(float(np.median([r["R2"] for r in st])), 4) if st else None}
    out = {
        "question": "甲基能否 recover 多-sSNV 遺傳 partition → 校正門檻給單 sSNV",
        "design": ("geno-string 群 = 遺傳 partition (LOH 區 = subclone 真值, germline ASM 已被 LOH 排除); "
                   "PERMANOVA 甲基距離~遺傳label + shuffle null; region>80kb 跳(單分子跨不完)"),
        "params": {"max_reads": MAX_READS, "min_shared_cpg": MIN_SHARED, "n_perm": N_PERM},
        "n_target_multi_sSNV": n_tgt,
        "n_too_wide_skipped": sum(1 for r in records if r["status"] == "region_too_wide(>80kb)"),
        "n_with_genetic_partition(>=2 geno grp)": sum(1 for r in records if r["n_geno_grp"] >= 2),
        "n_testable_permanova": len(tested),
        "n_methyl_recovers_genetic(p<0.05)": len(sig),
        "recover_rate_of_testable": round(len(sig) / len(tested), 4) if tested else None,
        "recover_by_cn": {"loh": cn_block("loh"), "neutral": cn_block("neutral")},
        "R2_distribution": {"median": round(float(np.median(r2s)), 4) if r2s else None,
                            "p75": round(float(np.percentile(r2s, 75)), 4) if r2s else None,
                            "max": round(float(np.max(r2s)), 4) if r2s else None,
                            "median_sig": round(float(np.median([r["R2"] for r in sig])), 4) if sig else None},
        "calibration_verdict": ("若 recover_rate 高+R2 顯著 → 甲基 track 遺傳 subclone → 可設門檻轉移單 sSNV; "
                                "若 recover_rate 低/R2~null → 甲基子結構 = double-dip, 不能外推單 sSNV 偵測 subclone"),
    }
    json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=1)
    tsv = out_path.replace(".json", "_perregion.tsv")
    if records:
        cols = list(records[0].keys())
        with open(tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for rc in records:
                f.write("\t".join(str(rc[c]) for c in cols) + "\n")
    print(f"DONE: target={n_tgt} testable={len(tested)} sig(recover)={len(sig)} "
          f"rate={out['recover_rate_of_testable']} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
