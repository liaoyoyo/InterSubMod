"""[HP tag 排解 cis-ASM 可行性] 回答「cis-ASM 能否用 HP tag 排解」。
正解設計 = within-HP 突變分層: 在 somatic 突變所在單倍型(som_HP)內, 比較「帶 ALT 的 read」vs「帶 REF 的 read」
→ germline 等位效應被 hold constant, 顯著=真 subclone(非 cis-ASM)。
量化「可行率」: 多少位點 2x2(som_HP x allele) 填滿能跑此測試; 失敗主因 LOH-單倍型 vs clonal-無REF。
輸入 _single_loci.jsonl (與 sm_single_locus_methyl 同集)。argv: chroms_or_ALL out_path。
"""
import sys
import json
from collections import Counter
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report/data"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2
WINDOW = 2500
MAX_READS = 160


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        idx = o[rank]; prev = min(prev, p[idx] * n / (rank + 1)); q[idx] = prev
    return q


def hpclass(hp):
    if hp in (1, 11):
        return "HP1"
    if hp in (2, 21):
        return "HP2"
    return "other"


def read_meth_allele_hp(a, pos, ref, alt):
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
    q = rq.get(pos - 1)
    if q is None or seq is None:
        return None
    b = seq[q].upper()
    allele = "A" if b == alt else ("R" if b == ref else "-")
    hp = a.get_tag("HP") if a.has_tag("HP") else None
    return meth, allele, hpclass(hp)


def per_cpg_sig(reads, idsA, idsB):
    allc = set()
    for rn in idsA + idsB:
        allc.update(reads[rn][0].keys())
    pv = []; dbs = []
    for cp in allc:
        a_ = [reads[rn][0][cp] for rn in idsA if cp in reads[rn][0]]
        b_ = [reads[rn][0][cp] for rn in idsB if cp in reads[rn][0]]
        if len(a_) >= 3 and len(b_) >= 3:
            try:
                _, p = mannwhitneyu(a_, b_, alternative="two-sided")
            except ValueError:
                continue
            pv.append(p); dbs.append(abs(np.mean(a_) - np.mean(b_)))
    if not pv:
        return 0, 0.0
    q = bh(pv)
    sig = sum(1 for k in range(len(pv)) if q[k] < 0.05 and dbs[k] >= DBETA_MIN)
    return sig, (max(dbs) if dbs else 0.0)


def main(chroms, out_path):
    rows = [json.loads(l) for l in open(f"{A}/_single_loci.jsonl")]
    if chroms != ["ALL"]:
        chset = set(chroms); rows = [r for r in rows if r["chrom"] in chset]
    tb = pysam.AlignmentFile(TBAM, "rb")
    n_tgt = len(rows)
    n_not_phaseable = n_loh_oneHP = n_clonal_noRef = n_lowcov = 0
    n_runnable = n_subclone_confirmed = n_cis_attributed = 0
    records = []
    for idx, r in enumerate(rows):
        chrom, pos, ref, alt = r["chrom"], int(r["pos"]), r["ref"], r["alt"]
        reads = {}
        for a in tb.fetch(chrom, max(0, pos - WINDOW), pos + WINDOW):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_allele_hp(a, pos, ref, alt)
            if res is None or res[1] == "-":
                continue
            reads[a.query_name] = res
            if len(reads) >= MAX_READS:
                break
        altreads = [rn for rn, v in reads.items() if v[1] == "A"]
        refreads = [rn for rn, v in reads.items() if v[1] == "R"]
        # somatic HP = mode HP-class among ALT reads (restricted to HP1/HP2)
        altc = Counter(reads[rn][2] for rn in altreads if reads[rn][2] in ("HP1", "HP2"))
        rec = {"locus": f"{chrom}:{pos}", "cn": r["cn"], "n_alt": len(altreads), "n_ref": len(refreads),
               "som_hp": "", "sameHP_alt": 0, "sameHP_ref": 0, "otherHP": 0,
               "status": "", "within_hp_sig": 0}
        if not altc:
            n_not_phaseable += 1; rec["status"] = "not_phaseable(ALT 無 HP1/HP2)"; records.append(rec); continue
        som_hp = altc.most_common(1)[0][0]
        other = "HP2" if som_hp == "HP1" else "HP1"
        sameHP_alt = [rn for rn in altreads if reads[rn][2] == som_hp]
        sameHP_ref = [rn for rn in refreads if reads[rn][2] == som_hp]
        otherHP = [rn for rn in (altreads + refreads) if reads[rn][2] == other]
        rec.update({"som_hp": som_hp, "sameHP_alt": len(sameHP_alt), "sameHP_ref": len(sameHP_ref),
                    "otherHP": len(otherHP)})
        # blocker classification
        if len(otherHP) < 3:
            n_loh_oneHP += 1; rec["status"] = "LOH/one-haplotype(otherHP<3)"
        if len(sameHP_ref) < 3:
            if len(otherHP) >= 3:
                n_clonal_noRef += 1; rec["status"] = rec["status"] or "clonal(sameHP 無 REF read)"
            elif rec["status"].startswith("LOH"):
                pass
        if len(sameHP_alt) >= 3 and len(sameHP_ref) >= 3:
            n_runnable += 1
            sig, mdb = per_cpg_sig(reads, sameHP_alt, sameHP_ref)
            rec["within_hp_sig"] = sig
            if sig:
                n_subclone_confirmed += 1; rec["status"] = "RUNNABLE→subclone-confirmed(within-HP sig)"
            else:
                n_cis_attributed += 1; rec["status"] = "RUNNABLE→null(within-HP 不顯著)"
        elif not rec["status"]:
            n_lowcov += 1; rec["status"] = "lowcov(sameHP cells<3)"
        records.append(rec)
        if idx % 200 == 0:
            print(f"{idx}/{n_tgt} runnable={n_runnable} confirmed={n_subclone_confirmed}", flush=True)
    tb.close()
    out = {"question": "cis-ASM 能否用 HP tag(within-HP 突變分層)排解?",
           "design": "som_HP 內 ALT-read vs REF-read 比甲基; 顯著=真 subclone(germline 等位 held constant)",
           "window": WINDOW, "n_target_single_sSNV_CNclean": n_tgt,
           "n_not_phaseable(ALT無HP1/HP2)": n_not_phaseable,
           "n_LOH_oneHaplotype(otherHP<3)": n_loh_oneHP,
           "n_clonal_no_sameHP_REF": n_clonal_noRef,
           "n_lowcov_other": n_lowcov,
           "n_RUNNABLE(within-HP test 可跑)": n_runnable,
           "runnable_rate": round(n_runnable / n_tgt, 4) if n_tgt else None,
           "of_runnable_subclone_confirmed": n_subclone_confirmed,
           "of_runnable_null(cis_or_no_signal)": n_cis_attributed,
           "verdict": ("HP 排解可行率 = runnable_rate; 失敗主因見 LOH/clonal 計數。"
                       "RUNNABLE 中 within-HP 顯著者 = HP 成功排解為真 subclone; 其餘為 cis/null。")}
    json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=1)
    tsv = out_path.replace(".json", "_perlocus.tsv")
    cols = list(records[0].keys()) if records else []
    with open(tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for rc in records:
            f.write("\t".join(str(rc[c]) for c in cols) + "\n")
    print(f"DONE: target={n_tgt} runnable={n_runnable} confirmed={n_subclone_confirmed} "
          f"LOH={n_loh_oneHP} clonal_noRef={n_clonal_noRef} not_phaseable={n_not_phaseable} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
