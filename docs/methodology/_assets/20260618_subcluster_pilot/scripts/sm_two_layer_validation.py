"""[兩層大規模驗證] 驗證「甲基資訊分兩層」假設在全基因組是否穩健（非只 49 corroborated 小 N）。
對全部 740 CN-clean 結構區，單次 fetch，同時算三軸：
 (A) HP-axis germline ASM = HP1{1,1-1,11} vs HP2{2,2-1,21} 每 CpG MWU+BH-FDR。**完全非循環**（germline 軸,無 somatic）
     → 大規模驗證「haplotype 層 = germline 等位特異甲基化」是否真實且普遍。
 (B) genotype-axis somatic = 兩大基因型 population per-CpG MWU+FDR（= 既有 corroboration,應重現 49）。
 (C) 兩軸分解 = genotype-sig CpG 中,有多少同時 HP-sig(germline-explained) vs HP-not-sig(subclone 殘差候選)。
判別「資訊是否有用」: germline ASM 普及率 + Δβ effect size + 對 phasing 的可用性(HP1/HP2 可分性)。
argv: chroms out_path。HP 比對用 str()（修 2026-06-29 str/int bug）。
"""
import sys
import json
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2
HP_MIN = 4        # 每 HP 至少幾條 read 才評估 germline ASM
CPG_MIN = 3       # 每 CpG 每組至少幾點才測


def bh(p):
    p = np.asarray(p)
    n = len(p)
    if n == 0:
        return np.array([])
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


def hp_germ(v):
    """germline haplotype class（str 比對,兼容 longphase-s 字串與 -to 整數）。"""
    s = str(v)
    if s in ("1", "1-1", "11"):
        return 1
    if s in ("2", "2-1", "21"):
        return 2
    return 0


def read_meth_geno_hp(a, som):
    try:
        mb = a.modified_bases
    except Exception:
        return None
    pairs = a.get_aligned_pairs(matches_only=True)
    qr = {q: r for q, r in pairs}
    rq = {r: q for q, r in pairs}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != 'm':
                continue
            for qpos, mlq in lst:
                r = qr.get(qpos)
                if r is not None:
                    meth[r] = mlq / 255.0
    geno = []
    seq = a.query_sequence
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None:
            geno.append("-")
            continue
        b = seq[q].upper()
        geno.append("A" if b == alt else ("R" if b == ref else "-"))
    hp = a.get_tag("HP") if a.has_tag("HP") else None
    return meth, "".join(geno), hp


def per_cpg_test(group_a, group_b, reads_meth):
    """回 (n_cpg_tested, n_sig, max_dbeta, sig_cpg_set)。group=read_name list。"""
    cps = set()
    for rn in group_a + group_b:
        cps |= set(reads_meth[rn].keys())
    pvals = []
    dbs = []
    cplist = []
    for cp in cps:
        av = [reads_meth[rn][cp] for rn in group_a if cp in reads_meth[rn]]
        bv = [reads_meth[rn][cp] for rn in group_b if cp in reads_meth[rn]]
        if len(av) >= CPG_MIN and len(bv) >= CPG_MIN:
            try:
                _, p = mannwhitneyu(av, bv, alternative="two-sided")
            except ValueError:
                continue
            pvals.append(p)
            dbs.append(abs(np.mean(av) - np.mean(bv)))
            cplist.append(cp)
    if not pvals:
        return 0, 0, 0.0, set()
    q = bh(pvals)
    sig = {cplist[k] for k in range(len(pvals)) if q[k] < 0.05 and dbs[k] >= DBETA_MIN}
    return len(pvals), len(sig), float(max(dbs)) if dbs else 0.0, sig


def main(chroms, out_path):
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    chset = set(chroms)
    tgt = [r for r in reg if (chroms == ["ALL"] or r["chrom"] in chset)
           and r["cn"] in ("loh", "neutral")
           and r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")
           and (r.get("n_populations") or 0) >= 2]
    records = []
    for idx, r in enumerate(tgt):
        chrom = r["chrom"]
        som = [(p, cen[f"{chrom}:{p}"]["ref"], cen[f"{chrom}:{p}"]["alt"])
               for p in range(r["start"], r["end"] + 1)
               if f"{chrom}:{p}" in cen and cen[f"{chrom}:{p}"].get("somatic")]
        if len(som) < 2:
            continue
        rec = {"region": r["region"], "chrom": chrom, "shape": r["tree_shape"], "cn": r["cn"],
               "n_reads": 0, "hp1_n": 0, "hp2_n": 0, "hp_evaluable": 0,
               "n_cpg_hp": 0, "n_sig_hp": 0, "max_dbeta_hp": 0.0, "germline_asm": 0,
               "geno_corrob": 0, "n_sig_geno": 0, "geno_sig_hp_explained": 0, "geno_sig_residual": 0,
               "subclone_candidate": 0}
        meth_all = {}
        hp_all = {}
        geno_full = {}
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_geno_hp(a, som)
            if res is None:
                continue
            meth, geno, hp = res
            rn = a.query_name
            meth_all[rn] = meth
            hp_all[rn] = hp_germ(hp)
            if "-" not in geno:
                geno_full[rn] = geno
        rec["n_reads"] = len(meth_all)
        # (A) HP-axis germline ASM（所有有 HP 的 read,不需 genotype 覆蓋 → 非循環）
        hp1 = [rn for rn in meth_all if hp_all[rn] == 1]
        hp2 = [rn for rn in meth_all if hp_all[rn] == 2]
        rec["hp1_n"], rec["hp2_n"] = len(hp1), len(hp2)
        hp_sig = set()
        if len(hp1) >= HP_MIN and len(hp2) >= HP_MIN:
            rec["hp_evaluable"] = 1
            ncpg, nsig, mdb, hp_sig = per_cpg_test(hp1, hp2, meth_all)
            rec["n_cpg_hp"], rec["n_sig_hp"], rec["max_dbeta_hp"] = ncpg, nsig, round(mdb, 3)
            rec["germline_asm"] = 1 if nsig > 0 else 0
        # (B) genotype-axis somatic（重現 corroboration）
        pop = defaultdict(list)
        for rn, g in geno_full.items():
            pop[g].append(rn)
        pops = sorted([(g, ids) for g, ids in pop.items() if len(ids) >= 3], key=lambda x: -len(x[1]))
        geno_sig = set()
        if len(pops) >= 2:
            ncpg2, nsig2, mdb2, geno_sig = per_cpg_test(pops[0][1], pops[1][1], meth_all)
            rec["n_sig_geno"] = nsig2
            rec["geno_corrob"] = 1 if nsig2 > 0 else 0
        # (C) 兩軸分解: genotype-sig CpG 中有多少 HP 也 sig(germline-explained) vs 殘差
        if geno_sig:
            expl = len(geno_sig & hp_sig)
            resid = len(geno_sig - hp_sig)
            rec["geno_sig_hp_explained"] = expl
            rec["geno_sig_residual"] = resid
            # subclone 候選 = corroborated + hp 可評估 + 過半 geno-sig CpG 為殘差(HP 無法解釋)
            if rec["hp_evaluable"] and resid >= max(1, len(geno_sig) // 2):
                rec["subclone_candidate"] = 1
        records.append(rec)
        if idx % 100 == 0:
            print(f"[{chroms[0] if chroms!=['ALL'] else 'ALL'}] {idx}/{len(tgt)} "
                  f"germ_asm={sum(x['germline_asm'] for x in records)} "
                  f"corrob={sum(x['geno_corrob'] for x in records)}", flush=True)
    tb.close()

    # ===== 整體統計 =====
    hp_eval = [r for r in records if r["hp_evaluable"]]
    germ = [r for r in hp_eval if r["germline_asm"]]
    corrob = [r for r in records if r["geno_corrob"]]
    cand = [r for r in records if r["subclone_candidate"]]
    corrob_hpeval = [r for r in corrob if r["hp_evaluable"]]
    corrob_explained = [r for r in corrob_hpeval if r["geno_sig_hp_explained"] >= max(1, r["n_sig_geno"] // 2)]

    def by_cn(lst, cn):
        return sum(1 for r in lst if r["cn"] == cn)

    out = {
        "design": "全 740 區三軸: (A)HP-axis germline ASM 非循環 (B)genotype-axis somatic (C)分解",
        "params": {"hp_min": HP_MIN, "cpg_min": CPG_MIN, "dbeta_min": DBETA_MIN, "fdr": 0.05},
        "n_regions": len(records),
        "haplotype_layer": {
            "hp_evaluable": len(hp_eval),
            "germline_asm_regions": len(germ),
            "germline_asm_prevalence_of_evaluable": round(len(germ) / len(hp_eval), 4) if hp_eval else None,
            "germ_asm_by_cn": {"loh": by_cn(germ, "loh"), "neutral": by_cn(germ, "neutral")},
            "hpeval_by_cn": {"loh": by_cn(hp_eval, "loh"), "neutral": by_cn(hp_eval, "neutral")},
            "max_dbeta_hp_median": round(float(np.median([r["max_dbeta_hp"] for r in germ])), 3) if germ else None,
            "n_sig_hp_median": float(np.median([r["n_sig_hp"] for r in germ])) if germ else None,
        },
        "somatic_axis": {
            "geno_corroborated": len(corrob),
            "corrob_hp_evaluable": len(corrob_hpeval),
            "corrob_germline_explained": len(corrob_explained),
            "corrob_explained_frac": round(len(corrob_explained) / len(corrob_hpeval), 4) if corrob_hpeval else None,
        },
        "two_layer_split": {
            "subclone_candidates": len(cand),
            "candidate_by_cn": {"loh": by_cn(cand, "loh"), "neutral": by_cn(cand, "neutral")},
            "candidate_regions": [r["region"] for r in cand],
        },
        "cross_tab": {
            "germ_asm_AND_corrob": sum(1 for r in records if r["germline_asm"] and r["geno_corrob"]),
            "germ_asm_NOT_corrob": sum(1 for r in records if r["germline_asm"] and not r["geno_corrob"]),
            "corrob_NOT_germ_asm": sum(1 for r in records if r["geno_corrob"] and r["hp_evaluable"] and not r["germline_asm"]),
        },
    }
    json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=1)
    tsv = out_path.replace(".json", "_perregion.tsv")
    if records:
        cols = list(records[0].keys())
        with open(tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for rc in records:
                f.write("\t".join(str(rc[c]) for c in cols) + "\n")
    hl = out["haplotype_layer"]
    print(f"DONE: regions={len(records)} hp_eval={hl['hp_evaluable']} "
          f"germline_asm={hl['germline_asm_regions']} (prev={hl['germline_asm_prevalence_of_evaluable']}) "
          f"corrob={len(corrob)} candidates={len(cand)} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
