"""[甲基錨點招募力測試] 回答「甲基錨點能否招募更多 read（含未穿 sSNV / unphased / HP3）歸群」。
對每個有 germline-ASM 錨點的區域，測三件事：
 (1) HP 招募力: 用錨點 CpG 的甲基,leave-one-out 最近質心分類 HP-tagged read → 招募準確度(能不能只靠甲基把 read 分到 HP)。
 (2) 可招募池: 無 HP tag(unphased/HP3=0)但覆蓋錨點 CpG 的 read 數 → 甲基能救回多少讀。
 (3) 同 HP 內 subclone: 某 HP 內若有≥2 個 somatic 基因型子群,甲基能否再分開它們(within-HP subclone 招募,你問的關鍵鏈)。
🔑 判準: HP 招募準確度高=甲基能招募到 haplotype 層(phasing 用); within-HP subclone 可分=甲基能招募到 subclone 層。
argv: chroms out_path。HP 比對用 str()(修 2026-06-29 bug)。
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
HP_MIN = 4
CPG_MIN = 3
ANCHOR_MIN = 3   # read 至少覆蓋幾個錨點 CpG 才可被招募/分類


def bh(p):
    p = np.asarray(p); n = len(p)
    if n == 0:
        return np.array([])
    o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]; prev = min(prev, p[i] * n / (rank + 1)); q[i] = prev
    return q


def hp_germ(v):
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
    geno = []
    seq = a.query_sequence
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None:
            geno.append("-")
        else:
            b = seq[q].upper()
            geno.append("A" if b == alt else ("R" if b == ref else "-"))
    hp = a.get_tag("HP") if a.has_tag("HP") else None
    return meth, "".join(geno), hp_germ(hp)


def anchor_cpgs(g1, g2, meth):
    """g1,g2=read name lists. 回 germline-ASM sig CpG set(HP1 vs HP2 差異顯著)。"""
    cps = set()
    for rn in g1 + g2:
        cps |= set(meth[rn].keys())
    pv = []; db = []; cl = []
    for cp in cps:
        a_ = [meth[rn][cp] for rn in g1 if cp in meth[rn]]
        b_ = [meth[rn][cp] for rn in g2 if cp in meth[rn]]
        if len(a_) >= CPG_MIN and len(b_) >= CPG_MIN:
            try:
                _, p = mannwhitneyu(a_, b_, alternative="two-sided")
            except ValueError:
                continue
            pv.append(p); db.append(abs(np.mean(a_) - np.mean(b_))); cl.append(cp)
    if not pv:
        return set()
    q = bh(pv)
    return {cl[k] for k in range(len(pv)) if q[k] < 0.05 and db[k] >= DBETA_MIN}


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
        meth = {}; hp = {}; geno = {}
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_geno_hp(a, som)
            if res is None:
                continue
            m, g, h = res
            rn = a.query_name
            meth[rn] = m; hp[rn] = h
            if "-" not in g:
                geno[rn] = g
        hp1 = [rn for rn in meth if hp[rn] == 1]
        hp2 = [rn for rn in meth if hp[rn] == 2]
        hp0 = [rn for rn in meth if hp[rn] == 0]   # unphased / HP3
        rec = {"region": r["region"], "cn": r["cn"], "shape": r["tree_shape"],
               "hp1_n": len(hp1), "hp2_n": len(hp2), "hp0_n": len(hp0),
               "n_anchor": 0, "loo_tested": 0, "loo_correct": 0, "loo_acc": "",
               "recruitable_hp0": 0, "within_hp_subclone_tested": 0, "within_hp_subclone_sep": 0}
        if len(hp1) < HP_MIN or len(hp2) < HP_MIN:
            records.append(rec); continue
        anc = anchor_cpgs(hp1, hp2, meth)
        rec["n_anchor"] = len(anc)
        if not anc:
            records.append(rec); continue
        anc = list(anc)
        # (1) HP 招募 LOO: 最近質心(錨點甲基均值)
        def prof(rn):
            return {c: meth[rn][c] for c in anc if c in meth[rn]}
        labeled = [(rn, 1) for rn in hp1] + [(rn, 2) for rn in hp2]
        for rn, true_hp in labeled:
            pr = prof(rn)
            if len(pr) < ANCHOR_MIN:
                continue
            # centroids from the OTHER reads (leave-one-out)
            c1 = {}; c2 = {}
            for c in anc:
                v1 = [meth[x][c] for x in hp1 if x != rn and c in meth[x]]
                v2 = [meth[x][c] for x in hp2 if x != rn and c in meth[x]]
                if v1:
                    c1[c] = np.mean(v1)
                if v2:
                    c2[c] = np.mean(v2)
            shared1 = [c for c in pr if c in c1]; shared2 = [c for c in pr if c in c2]
            if len(shared1) < ANCHOR_MIN or len(shared2) < ANCHOR_MIN:
                continue
            d1 = np.mean([abs(pr[c] - c1[c]) for c in shared1])
            d2 = np.mean([abs(pr[c] - c2[c]) for c in shared2])
            pred = 1 if d1 <= d2 else 2
            rec["loo_tested"] += 1
            rec["loo_correct"] += int(pred == true_hp)
        if rec["loo_tested"] > 0:
            rec["loo_acc"] = round(rec["loo_correct"] / rec["loo_tested"], 3)
        # (2) 可招募池: hp0 read 覆蓋≥ANCHOR_MIN 錨點
        rec["recruitable_hp0"] = sum(1 for rn in hp0 if len([c for c in anc if c in meth[rn]]) >= ANCHOR_MIN)
        # (3) within-HP subclone: 某 HP 內 ≥2 geno 子群, 甲基能否分開
        for hpg, hpreads in [(1, hp1), (2, hp2)]:
            sub = defaultdict(list)
            for rn in hpreads:
                if rn in geno:
                    sub[geno[rn]].append(rn)
            subs = [ids for ids in sub.values() if len(ids) >= 3]
            if len(subs) >= 2:
                rec["within_hp_subclone_tested"] += 1
                sep = anchor_cpgs(subs[0], subs[1], meth)
                if sep:
                    rec["within_hp_subclone_sep"] += 1
        records.append(rec)
        if idx % 100 == 0:
            print(f"[{chroms[0] if chroms!=['ALL'] else 'ALL'}] {idx}/{len(tgt)} anchored={sum(1 for x in records if x['n_anchor']>0)}", flush=True)
    tb.close()

    anchored = [r for r in records if r["n_anchor"] > 0]
    loo = [r for r in anchored if r["loo_tested"] >= 5]
    tot_t = sum(r["loo_tested"] for r in loo); tot_c = sum(r["loo_correct"] for r in loo)
    wsub_t = sum(r["within_hp_subclone_tested"] for r in records)
    wsub_s = sum(r["within_hp_subclone_sep"] for r in records)
    out = {
        "design": "甲基錨點招募力: (1)HP LOO 最近質心分類準確度 (2)可招募 unphased/HP3 池 (3)within-HP subclone 可分性",
        "params": {"hp_min": HP_MIN, "anchor_min": ANCHOR_MIN, "dbeta_min": DBETA_MIN},
        "n_regions": len(records),
        "n_anchored_regions": len(anchored),
        "hp_recruitment": {
            "regions_with_LOO(>=5 reads)": len(loo),
            "total_reads_classified": tot_t,
            "total_correct": tot_c,
            "overall_accuracy": round(tot_c / tot_t, 4) if tot_t else None,
            "per_region_acc_median": round(float(np.median([r["loo_acc"] for r in loo if r["loo_acc"] != ""])), 3) if loo else None,
        },
        "recruitable_pool": {
            "total_unphased_hp3_recruitable": sum(r["recruitable_hp0"] for r in anchored),
            "total_hp0_in_anchored": sum(r["hp0_n"] for r in anchored),
            "total_hp_tagged_in_anchored": sum(r["hp1_n"] + r["hp2_n"] for r in anchored),
        },
        "within_hp_subclone": {
            "regions_tested": wsub_t,
            "regions_methyl_separable": wsub_s,
            "separable_frac": round(wsub_s / wsub_t, 4) if wsub_t else None,
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
    hr = out["hp_recruitment"]; rp = out["recruitable_pool"]; ws = out["within_hp_subclone"]
    print(f"DONE: anchored={len(anchored)} HP_acc={hr['overall_accuracy']} "
          f"recruitable_hp0={rp['total_unphased_hp3_recruitable']} "
          f"within_hp_subclone_sep={ws['regions_methyl_separable']}/{ws['regions_tested']} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
