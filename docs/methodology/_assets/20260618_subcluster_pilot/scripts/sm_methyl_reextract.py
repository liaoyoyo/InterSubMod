"""[A1 甲基 genome-wide 重抽 — 補 L4 缺口] 直接從 tumor BAM 的 MM/ML(5mC) 重抽甲基, 對齊連鎖區域(非 ISM anchor 窗)。
已驗證: 對 chr17 與 ISM methylation.csv corr=1.000。
對每個 CN-clean 結構區(>=2 population): 單次 fetch → 每 read 同時取 (a)5mC 甲基 (b)som sSNV 基因型 (c)HP tag。
genotype-anchored: 兩大基因型 population 間 per-CpG MWU + BH-FDR(sig: q<0.05 且 |Δβ|>=0.2)。
cis-ASM control: sig CpG 在 HP1 vs HP2(germline 軸) 是否也區分 → 是=germline cis 非 subclone。
argv: chroms out_path。🔴 甲基 = corroborate 非 detect。
"""
import sys
import json
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2


def bh(p):
    p = np.asarray(p)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


def read_meth_and_geno(a, som):
    """單次: 回 (meth{refpos:beta}, geno_str over som, hp)。"""
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


def main(chroms, out_path):
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    cen = json.load(open(f"{A}/sm_linkage_genomewide.json"))["census"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    chset = set(chroms)
    tgt = [r for r in reg if r["chrom"] in chset and r["cn"] in ("loh", "neutral")
           and r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")
           and (r.get("n_populations") or 0) >= 2]
    n_tested = n_corrob = n_nopop = 0
    cis_explained = 0
    sig_counts = []
    examples = []
    for idx, r in enumerate(tgt):
        chrom = r["chrom"]
        som = [(p, cen[f"{chrom}:{p}"]["ref"], cen[f"{chrom}:{p}"]["alt"])
               for p in range(r["start"], r["end"] + 1)
               if f"{chrom}:{p}" in cen and cen[f"{chrom}:{p}"].get("somatic")]
        if len(som) < 2:
            continue
        poslist = [s[0] for s in som]
        reads = {}  # read_name -> (meth, geno, hp)
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_and_geno(a, som)
            if res is None:
                continue
            meth, geno, hp = res
            if "-" in geno:   # 只用覆蓋全部 som 的 read
                continue
            reads[a.query_name] = (meth, geno, hp)
        pop = defaultdict(list)
        for rn, (meth, geno, hp) in reads.items():
            pop[geno].append(rn)
        pops = [(g, ids) for g, ids in pop.items() if len(ids) >= 3]
        if len(pops) < 2:
            n_nopop += 1
            continue
        pops.sort(key=lambda x: -len(x[1]))
        idsA = set(pops[0][1]); idsB = set(pops[1][1])
        n_tested += 1
        # all CpG ref positions seen
        allcpg = set()
        for rn in (idsA | idsB):
            allcpg.update(reads[rn][0].keys())
        pv = []; cps = []; dbs = []
        for cp in allcpg:
            a_ = [reads[rn][0][cp] for rn in idsA if cp in reads[rn][0]]
            b_ = [reads[rn][0][cp] for rn in idsB if cp in reads[rn][0]]
            if len(a_) >= 3 and len(b_) >= 3:
                try:
                    _, p = mannwhitneyu(a_, b_, alternative="two-sided")
                except ValueError:
                    continue
                pv.append(p); cps.append(cp); dbs.append(abs(np.mean(a_) - np.mean(b_)))
        if not pv:
            continue
        q = bh(pv)
        sig = [(cps[k], dbs[k]) for k in range(len(pv)) if q[k] < 0.05 and dbs[k] >= DBETA_MIN]
        sig_counts.append(len(sig))
        if sig:
            n_corrob += 1
            # cis control: HP1({1,1-1}) vs HP2({2,2-1}) at sig CpGs
            hp1 = [rn for rn in (idsA | idsB) if reads[rn][2] in (1, 11)]
            hp2 = [rn for rn in (idsA | idsB) if reads[rn][2] in (2, 21)]
            cis = 0
            if len(hp1) >= 3 and len(hp2) >= 3:
                for cp, _ in sig:
                    a_ = [reads[rn][0][cp] for rn in hp1 if cp in reads[rn][0]]
                    b_ = [reads[rn][0][cp] for rn in hp2 if cp in reads[rn][0]]
                    if len(a_) >= 3 and len(b_) >= 3:
                        try:
                            _, p = mannwhitneyu(a_, b_, alternative="two-sided")
                            if p < 0.05 and abs(np.mean(a_) - np.mean(b_)) >= DBETA_MIN:
                                cis += 1
                        except ValueError:
                            pass
            if cis >= max(1, len(sig) // 2):
                cis_explained += 1
            if len(examples) < 25:
                examples.append({"region": r["region"], "shape": r["tree_shape"], "n_sig_cpg": len(sig),
                                 "popA": f"{pops[0][0]}({len(idsA)})", "popB": f"{pops[1][0]}({len(idsB)})", "cis_hit": cis})
        if idx % 100 == 0:
            print(f"[{chroms[0]}..] {idx}/{len(tgt)} tested={n_tested} corrob={n_corrob}", flush=True)
    tb.close()
    out = {"chroms": chroms, "n_target": len(tgt), "n_insufficient_pop": n_nopop, "n_tested": n_tested,
           "n_corroborated": n_corrob, "n_cis_explained": cis_explained,
           "sig_cpg_counts": sig_counts, "examples": examples}
    json.dump(out, open(out_path, "w"), ensure_ascii=False)
    print(f"DONE {chroms}: target={len(tgt)} tested={n_tested} corrob={n_corrob} cis_expl={cis_explained} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
