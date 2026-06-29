"""[單 sSNV 甲基外推測試] 回答「甲基能否外推到只有單個 sSNV(n_links=0)的位點」。
對每個 single(無連鎖)somatic 位點: 用該位點 REF-reads vs ALT-reads 做 genotype-anchored 甲基 per-CpG MWU+BH-FDR。
= 單位點 ASM。同時做 HP1-vs-HP2 cis-control（看單位點是否比多-sSNV 更可評估 cis）。
與 multi-sSNV(L4/L8) 直接可比。argv: chroms_or_ALL out_path。
🔴 單位點無連鎖錨 → REF/ALT 甲基差異 = subclone / cis-ASM / 直接效應 三者不可分（最 confounded）。
"""
import sys
import json
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report/data"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2
WINDOW = 2500
MAX_READS = 160  # 深 pileup subsample 上限（median cov=24，cap 幾乎不影響結果、砍掉慢尾）


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        idx = o[rank]; prev = min(prev, p[idx] * n / (rank + 1)); q[idx] = prev
    return q


def read_meth_allele(a, pos, ref, alt):
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
    return meth, allele, hp


def main(chroms, out_path):
    rows = [r for r in (json.loads(l) for l in open(f"{A}/_single_loci.jsonl"))]
    if chroms != ["ALL"]:
        chset = set(chroms); rows = [r for r in rows if r["chrom"] in chset]
    tb = pysam.AlignmentFile(TBAM, "rb")
    n_tgt = len(rows); n_testable = n_asm = n_nopop = 0
    hp_eval = 0; sig_counts = []; records = []
    for idx, r in enumerate(rows):
        chrom, pos, ref, alt = r["chrom"], int(r["pos"]), r["ref"], r["alt"]
        reads = {}
        for a in tb.fetch(chrom, max(0, pos - WINDOW), pos + WINDOW):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_allele(a, pos, ref, alt)
            if res is None or res[1] == "-":
                continue
            reads[a.query_name] = res
            if len(reads) >= MAX_READS:
                break
        idsA = [rn for rn, v in reads.items() if v[1] == "A"]
        idsR = [rn for rn, v in reads.items() if v[1] == "R"]
        rec = {"locus": f"{chrom}:{pos}", "cn": r["cn"], "n_alt": len(idsA), "n_ref": len(idsR),
               "n_cpg_tested": 0, "n_sig_cpg": 0, "max_dbeta": 0.0, "asm": 0, "hp_control_eval": 0, "cis_hit": 0}
        if len(idsA) < 3 or len(idsR) < 3:
            n_nopop += 1; records.append(rec); continue
        n_testable += 1
        allc = set()
        for rn in idsA + idsR:
            allc.update(reads[rn][0].keys())
        pv = []; cps = []; dbs = []
        for cp in allc:
            a_ = [reads[rn][0][cp] for rn in idsA if cp in reads[rn][0]]
            b_ = [reads[rn][0][cp] for rn in idsR if cp in reads[rn][0]]
            if len(a_) >= 3 and len(b_) >= 3:
                try:
                    _, p = mannwhitneyu(a_, b_, alternative="two-sided")
                except ValueError:
                    continue
                pv.append(p); cps.append(cp); dbs.append(abs(np.mean(a_) - np.mean(b_)))
        rec["n_cpg_tested"] = len(pv)
        if not pv:
            records.append(rec); continue
        q = bh(pv)
        rec["max_dbeta"] = round(float(max(dbs)), 3)
        sig = [(cps[k], dbs[k]) for k in range(len(pv)) if q[k] < 0.05 and dbs[k] >= DBETA_MIN]
        sig_counts.append(len(sig)); rec["n_sig_cpg"] = len(sig)
        if sig:
            n_asm += 1; rec["asm"] = 1
            # FIX(2026-06-29): HP:Z 回 str → 字串比對（原 int 比對恆 False）
            hp1 = [rn for rn in idsA + idsR if str(reads[rn][2]) in ("1", "1-1", "11")]
            hp2 = [rn for rn in idsA + idsR if str(reads[rn][2]) in ("2", "2-1", "21")]
            if len(hp1) >= 3 and len(hp2) >= 3:
                hp_eval += 1; rec["hp_control_eval"] = 1
                cis = 0
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
                rec["cis_hit"] = cis
        records.append(rec)
        if idx % 200 == 0:
            print(f"{idx}/{n_tgt} testable={n_testable} asm={n_asm}", flush=True)
    tb.close()
    out = {"scope": "single-sSNV (n_links=0) CN-clean somatic loci; REF-vs-ALT 單位點 ASM",
           "window": WINDOW, "n_target": n_tgt, "n_insufficient_pop": n_nopop, "n_testable": n_testable,
           "n_asm_positive": n_asm, "asm_rate_of_testable": round(n_asm / n_testable, 4) if n_testable else None,
           "hp_control_evaluable": hp_eval, "median_sig_cpg": float(np.median(sig_counts)) if sig_counts else 0}
    json.dump(out, open(out_path, "w"), ensure_ascii=False)
    tsv = out_path.replace(".json", "_perlocus.tsv")
    if records:
        cols = list(records[0].keys())
        with open(tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for rc in records:
                f.write("\t".join(str(rc[c]) for c in cols) + "\n")
    print(f"DONE: target={n_tgt} testable={n_testable} asm={n_asm} hp_eval={hp_eval} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
