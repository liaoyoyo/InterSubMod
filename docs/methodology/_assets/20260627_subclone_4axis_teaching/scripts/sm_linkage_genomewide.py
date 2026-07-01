#!/usr/bin/env python3
"""[S0-S2 全基因組 sSNV 單分子連鎖 — 完整枚舉 + 高效 single-pass]

S0 完整宇宙 : union TP∪FP per-chrom (filtered_snv_{tp,fp}_{chrom}.vcf.gz)，TP 優先；每 sSNV 標 source。
S1 高效枚舉 : 每 chrom sorted → 以「相鄰 gap > TIER_R」切 linkage group；group 內只配 dist<=TIER_R 的對
              (same-read 候選 = Tier-R)。group-based 讓每個 sSNV 只 pileup 一次(非每對重 pile = 效率核心)。
S2 共現分類 : per-group per-read 等位向量 → 每對 2×2 (RR/RA/AR/AA) → classify
              (mutual_excl/co_linked/nested_a_in_b/nested_b_in_a/independent/sparse) + HP 一致 gate
              (同 germline HP=subclone, 異 HP=allelic) + normal=REF somatic 確認 + PS 記錄。

完整輸出: 全 per-pair 連鎖地圖(非 top-40) + per-sSNV census(每 sSNV 一筆，供 S4 帳本) + tally + universe。
數字全落 sm_linkage_genomewide.json，下游 §13-A 注入，不手打。Tier-PS(same-PS) 由 sm_linkage_phaseset.py 另算。
"""
import sys
import os
import json
import gzip
import time
from collections import Counter, defaultdict
from bisect import bisect_right
import pysam

# per-sample 化(2026-06-29):env var 覆寫,預設 HCC1395(regression-safe)。
# SM_VCF_MODE=perchrom(預設,filtered_snv_{src}_{chrom}.vcf.gz) 或 single(filtered_snv_{src}.vcf.gz,tabix fetch)
VD = os.environ.get("SM_VD", "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup")
TBAM = os.environ.get("SM_TBAM", "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam")
NBAM = os.environ.get("SM_NBAM", "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam")
A = os.environ.get("SM_WORKDIR", "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot")
OUT = os.environ.get("SM_OUT", f"{A}/sm_linkage_genomewide.json")
VCF_MODE = os.environ.get("SM_VCF_MODE", "perchrom")

TIER_R = 50000      # same-read 枚舉上限 (read aligned-span p95~34kb, max~76kb → 50kb 涵蓋, 由 co-read 過濾)
COREAD_MIN = 6      # 連鎖「confirmed」最低共讀 (記錄 >=2 的所有對, powered = coread>=6)
MAPQ_MIN = 20
NORM_REF_MIN = 4    # normal somatic 確認: norm_ref>=NORM_REF_MIN 且 normal VAF<NORM_VAF_MAX
NORM_VAF_MAX = 0.05  # 容忍 normal 端低階噪聲 ALT (chr17 β1 真 somatic 但 normal 56/2=3.4%; 嚴格 ==0 會誤殺)
CHROMS = [f"chr{c}" for c in range(1, 23)]


def is_somatic(nref, nalt):
    if nref < NORM_REF_MIN:
        return False
    tot = nref + nalt
    return tot > 0 and (nalt / tot) < NORM_VAF_MAX


def _vcf_lines(src, chrom):
    """yield VCF 資料行。perchrom: 整檔讀;single: tabix fetch 該 chrom。"""
    if VCF_MODE == "single":
        path = f"{VD}/filtered_snv_{src}.vcf.gz"
        if not os.path.exists(path):
            return
        idx = path + ".tbi" if os.path.exists(path + ".tbi") else (path + ".csi" if os.path.exists(path + ".csi") else None)
        try:
            tb = pysam.TabixFile(path, index=idx) if idx else pysam.TabixFile(path)
            for ln in tb.fetch(chrom):
                yield ln
        except (ValueError, OSError):
            return  # 該 chrom 不在 VCF
    else:
        path = f"{VD}/filtered_snv_{src}_{chrom}.vcf.gz"
        if not os.path.exists(path):
            return
        with gzip.open(path, "rt") as fh:
            for ln in fh:
                if not ln.startswith("#"):
                    yield ln


def load_union(chrom):
    """回 sorted list of (pos, ref, alt, source). TP 優先 (同 pos 不覆蓋)。"""
    snvs = {}
    for src in ("tp", "fp"):
        for ln in _vcf_lines(src, chrom):
            p = ln.split("\t")
            pos = int(p[1])
            if pos not in snvs:
                snvs[pos] = (p[3].upper(), p[4].strip().upper(), src.upper())
    return sorted((pos, v[0], v[1], v[2]) for pos, v in snvs.items())


def per_read_calls(bam, chrom, group):
    """group=[(pos,ref,alt,src)]. 回 calls{rn:{pos:'REF'/'ALT'/'OTHER'}}, hp{rn:int}, ps{rn:int}.
    O2(2026-07-01): group span 發『單次』region pileup,只在 sSNV 位點 tabulate — 取代原逐位點 pileup。
    長 read 橫跨多 sSNV 時只解壓一次(非 k 次)。byte-identical 已驗證(original vs 本版,同輸入):
    397 group unit + chr22 端到端 7/7(pairs/census/tally/universe/計數)+ 多片切割合併 7/7;
    加速 per-unit 1.4-1.7x、超寬 group(chr8 503kb/33-sSNV)2.98-6.19x。原逐位點版見 git 歷史。"""
    calls = defaultdict(dict)
    hp = {}
    ps = {}
    refalt = {pos: (ref, alt) for pos, ref, alt, _ in group}
    posset = set(refalt)
    lo = min(posset)
    hi = max(posset)
    for col in bam.pileup(chrom, lo - 1, hi, truncate=True, min_base_quality=0, stepper="samtools"):
        pos = col.reference_pos + 1  # pileup 0-based reference_pos → 1-based sSNV pos
        if pos not in posset:
            continue
        ref, alt = refalt[pos]
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if aln.mapping_quality < MAPQ_MIN:
                continue
            rn = aln.query_name
            base = aln.query_sequence[pr.query_position].upper()
            calls[rn][pos] = "REF" if base == ref else ("ALT" if base == alt else "OTHER")
            if aln.has_tag("HP"):
                hp[rn] = aln.get_tag("HP")
            if aln.has_tag("PS"):
                ps[rn] = aln.get_tag("PS")
    return calls, hp, ps


def normal_status(nb, chrom, pos, ref, alt):
    nref = nalt = 0
    for col in nb.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="samtools"):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            if pr.alignment.mapping_quality < MAPQ_MIN:
                continue
            b = pr.alignment.query_sequence[pr.query_position].upper()
            if b == ref:
                nref += 1
            elif b == alt:
                nalt += 1
    return nref, nalt


def classify(rr, ra, ar, aa):
    if aa < 2:
        return "mutual_excl" if (ra + ar) >= 2 else "sparse"
    if ra == 0 and ar == 0:
        return "co_linked"
    if ar == 0:
        return "nested_a_in_b"
    if ra == 0:
        return "nested_b_in_a"
    return "independent"


def make_groups(snvs):
    """以相鄰 gap > TIER_R 切 group。回 list of groups (each = list of (pos,ref,alt,src))。"""
    groups = []
    cur = []
    for s in snvs:
        if cur and s[0] - cur[-1][0] > TIER_R:
            groups.append(cur)
            cur = []
        cur.append(s)
    if cur:
        groups.append(cur)
    return groups


def main(chroms=None, out_path=None):
    chroms = chroms or CHROMS
    out_path = out_path or OUT
    tb = pysam.AlignmentFile(TBAM, "rb")
    nb = pysam.AlignmentFile(NBAM, "rb")
    pairs = []
    census = {}          # "chrom:pos" -> dict
    tally = Counter()
    uni = {"n_tp": 0, "n_fp": 0, "n_union": 0, "per_chrom": {}}
    n_groups_multi = 0
    n_singletons = 0
    t0 = time.time()

    for chrom in chroms:
        snvs = load_union(chrom)
        ntp = sum(1 for s in snvs if s[3] == "TP")
        nfp = sum(1 for s in snvs if s[3] == "FP")
        uni["n_tp"] += ntp
        uni["n_fp"] += nfp
        uni["n_union"] += len(snvs)
        uni["per_chrom"][chrom] = {"tp": ntp, "fp": nfp, "union": len(snvs)}
        groups = make_groups(snvs)
        nstat = {}  # pos -> (nref, nalt) cache (only computed for linkable sSNV)

        for g in groups:
            poslist = [s[0] for s in g]
            # census skeleton for every sSNV (singleton or in-group)
            for s in g:
                key = f"{chrom}:{s[0]}"
                census[key] = {"chrom": chrom, "pos": s[0], "ref": s[1], "alt": s[2], "src": s[3],
                               "n_partners_le50k": 0, "n_realized_links": 0, "max_coread": 0,
                               "normal_ref": None, "normal_alt": None, "somatic": None,
                               "tumor_ref": None, "tumor_alt": None, "vaf": None,
                               "group_size": len(g)}
            if len(g) < 2:
                n_singletons += 1
                continue
            n_groups_multi += 1
            calls, hp, ps = per_read_calls(tb, chrom, g)
            # per-sSNV tumor depth/VAF (for S3 Foltz VAF-gradient)
            tdepth = defaultdict(lambda: [0, 0])
            for c in calls.values():
                for pp, v in c.items():
                    if v == "REF":
                        tdepth[pp][0] += 1
                    elif v == "ALT":
                        tdepth[pp][1] += 1
            for s in g:
                k = f"{chrom}:{s[0]}"
                tr, ta = tdepth.get(s[0], [0, 0])
                census[k]["tumor_ref"] = tr
                census[k]["tumor_alt"] = ta
                census[k]["vaf"] = round(ta / (tr + ta), 3) if (tr + ta) > 0 else None

            # enumerate pairs within TIER_R distance
            for i in range(len(g)):
                ai, aref, aalt, asrc = g[i]
                ja = i + 1
                while ja < len(g) and poslist[ja] - ai <= TIER_R:
                    bj, bref, balt, bsrc = g[ja]
                    census[f"{chrom}:{ai}"]["n_partners_le50k"] += 1
                    census[f"{chrom}:{bj}"]["n_partners_le50k"] += 1
                    # co-reads covering both with clean REF/ALT
                    multi = {rn: c for rn, c in calls.items()
                             if c.get(ai) in ("REF", "ALT") and c.get(bj) in ("REF", "ALT")}
                    coread = len(multi)
                    if coread >= 2:
                        tab = Counter((c[ai], c[bj]) for c in multi.values())
                        rr = tab[("REF", "REF")]; ra = tab[("REF", "ALT")]
                        ar = tab[("ALT", "REF")]; aa = tab[("ALT", "ALT")]
                        rel = classify(rr, ra, ar, aa)
                        hpa = Counter(hp[rn] for rn, c in calls.items() if c.get(ai) == "ALT" and rn in hp)
                        hpb = Counter(hp[rn] for rn, c in calls.items() if c.get(bj) == "ALT" and rn in hp)
                        da = hpa.most_common(1)[0][0] if hpa else None
                        db = hpb.most_common(1)[0][0] if hpb else None
                        same_hp = (da is not None and da == db)
                        # normal somatic confirm (cache)
                        for pp, rf, al in ((ai, aref, aalt), (bj, bref, balt)):
                            if pp not in nstat:
                                nstat[pp] = normal_status(nb, chrom, pp, rf, al)
                        na_r, na_a = nstat[ai]; nb_r, nb_a = nstat[bj]
                        som_a = is_somatic(na_r, na_a)
                        som_b = is_somatic(nb_r, nb_a)
                        for k, nr, nal, som in ((f"{chrom}:{ai}", na_r, na_a, som_a),
                                                (f"{chrom}:{bj}", nb_r, nb_a, som_b)):
                            census[k]["normal_ref"] = nr; census[k]["normal_alt"] = nal; census[k]["somatic"] = som
                            census[k]["max_coread"] = max(census[k]["max_coread"], coread)
                        powered = coread >= COREAD_MIN
                        if powered:
                            census[f"{chrom}:{ai}"]["n_realized_links"] += 1
                            census[f"{chrom}:{bj}"]["n_realized_links"] += 1
                        rec = {"chrom": chrom, "a": ai, "b": bj, "dist": bj - ai, "coread": coread,
                               "cells": {"RR": rr, "RA": ra, "AR": ar, "AA": aa}, "rel": rel,
                               "hp_a": da, "hp_b": db, "same_hp": same_hp,
                               "norm_a": [na_r, na_a], "norm_b": [nb_r, nb_a],
                               "somatic_a": som_a, "somatic_b": som_b, "src_a": asrc, "src_b": bsrc,
                               "powered": powered}
                        pairs.append(rec)
                        # tally only powered + both somatic-confirmed
                        if powered and som_a and som_b:
                            if rel in ("nested_a_in_b", "nested_b_in_a"):
                                tally["nested_sameHP" if same_hp else "nested_diffHP"] += 1
                            elif rel == "mutual_excl":
                                tally["sibling_sameHP" if same_hp else "allelic_diffHP"] += 1
                            else:
                                tally[rel] += 1
                    ja += 1
        print(f"[{chrom}] union={len(snvs)} groups_multi={n_groups_multi} singletons={n_singletons} "
              f"pairs={len(pairs)} t={time.time()-t0:.0f}s", flush=True)

    tb.close(); nb.close()
    out = {
        "params": {"TIER_R": TIER_R, "COREAD_MIN": COREAD_MIN, "MAPQ_MIN": MAPQ_MIN, "NORM_REF_MIN": NORM_REF_MIN},
        "universe": uni,
        "n_groups_multi": n_groups_multi, "n_singletons": n_singletons,
        "n_pairs_recorded": len(pairs),
        "tally_powered_somatic": dict(tally),
        "pairs": pairs,
        "census": census,
        "elapsed_sec": round(time.time() - t0, 1),
    }
    json.dump(out, open(out_path, "w"), ensure_ascii=False)
    print(f"DONE union={uni['n_union']} (TP={uni['n_tp']} FP={uni['n_fp']}) "
          f"pairs={len(pairs)} tally={dict(tally)} -> {out_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        main(sys.argv[1].split(","), sys.argv[2])  # chroms, out_path
    else:
        main()
