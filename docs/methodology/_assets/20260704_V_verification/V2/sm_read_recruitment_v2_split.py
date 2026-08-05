"""[V2 unphase/HP3 分樹決策] 重寫自 sm_read_recruitment.py。
核心修正: 原 hp_germ() 把 HP tag "3"(HP3) 與無 tag(unphase) 都映射成 0 → hp0 池混淆。
本版把 hp0 拆成 hp3(HP tag 家族 3, 經過 ALT) 與 unphase(無 HP tag, REF-only 難處理),
並對每條 read 記錄是否帶 somatic ALT (PI 提出的判別軸: 經過 ALT → 可併入帶該突變的 lineage)。
argv: chroms out_path。所有數字來自實際 BAM 掃描。
"""
import sys, json
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ_MIN = 20
DBETA_MIN = 0.2
HP_MIN = 4
CPG_MIN = 3
ANCHOR_MIN = 3


def bh(p):
    p = np.asarray(p); n = len(p)
    if n == 0:
        return np.array([])
    o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]; prev = min(prev, p[i] * n / (rank + 1)); q[i] = prev
    return q


def hp_class(v):
    """把 HP tag 分成 hp1/hp2/hp3/unphase (原 hp_germ 把 3 與 no-tag 混為 0)。"""
    if v is None:
        return "unphase"
    s = str(v)
    if s in ("1", "1-1", "11"):
        return "hp1"
    if s in ("2", "2-1", "21"):
        return "hp2"
    if s in ("3", "3-1", "31") or s.startswith("3"):
        return "hp3"
    return "other"


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
    hp_raw = a.get_tag("HP") if a.has_tag("HP") else None
    return meth, "".join(geno), hp_class(hp_raw), hp_raw


def anchor_cpgs(g1, g2, meth):
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
    hp_census = Counter()          # 全部 raw HP tag 分佈 (跨 tgt regions, 過濾後 reads)
    for idx, r in enumerate(tgt):
        chrom = r["chrom"]
        som = [(p, cen[f"{chrom}:{p}"]["ref"], cen[f"{chrom}:{p}"]["alt"])
               for p in range(r["start"], r["end"] + 1)
               if f"{chrom}:{p}" in cen and cen[f"{chrom}:{p}"].get("somatic")]
        if len(som) < 2:
            continue
        meth = {}; cls = {}; geno = {}; hasalt = {}; nalt = {}
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_geno_hp(a, som)
            if res is None:
                continue
            m, g, klass, raw = res
            rn = a.query_name
            meth[rn] = m; cls[rn] = klass
            hp_census[str(raw)] += 1
            hasalt[rn] = ("A" in g)          # 帶任一 somatic ALT
            nalt[rn] = g.count("A")
            if "-" not in g:
                geno[rn] = g
        hp1 = [rn for rn in meth if cls[rn] == "hp1"]
        hp2 = [rn for rn in meth if cls[rn] == "hp2"]
        hp3 = [rn for rn in meth if cls[rn] == "hp3"]
        unp = [rn for rn in meth if cls[rn] == "unphase"]

        def alt_ct(pool):
            return sum(1 for rn in pool if hasalt.get(rn))

        rec = {"region": r["region"], "cn": r["cn"], "shape": r["tree_shape"],
               "hp1_n": len(hp1), "hp2_n": len(hp2), "hp3_n": len(hp3), "unphase_n": len(unp),
               "hp1_alt_n": alt_ct(hp1), "hp2_alt_n": alt_ct(hp2),
               "hp3_alt_n": alt_ct(hp3), "unphase_alt_n": alt_ct(unp),
               "n_anchor": 0, "loo_tested": 0, "loo_correct": 0, "loo_acc": "",
               "recruitable_hp3": 0, "recruitable_unphase": 0,
               "recruitable_hp3_alt": 0, "recruitable_unphase_alt": 0,
               "hp3_assign_hp1": 0, "hp3_assign_hp2": 0,
               "unphase_assign_hp1": 0, "unphase_assign_hp2": 0,
               "within_hp_subclone_tested": 0, "within_hp_subclone_sep": 0,
               "within_hp3_subclone_tested": 0, "within_hp3_subclone_sep": 0}
        if len(hp1) < HP_MIN or len(hp2) < HP_MIN:
            records.append(rec); continue
        anc = anchor_cpgs(hp1, hp2, meth)
        rec["n_anchor"] = len(anc)
        if not anc:
            records.append(rec); continue
        anc = list(anc)

        def prof(rn):
            return {c: meth[rn][c] for c in anc if c in meth[rn]}

        # (1) HP1/HP2 LOO 最近質心 (原始 96.3% 招募力, 不動)
        labeled = [(rn, 1) for rn in hp1] + [(rn, 2) for rn in hp2]
        for rn, true_hp in labeled:
            pr = prof(rn)
            if len(pr) < ANCHOR_MIN:
                continue
            c1 = {}; c2 = {}
            for cc in anc:
                v1 = [meth[x][cc] for x in hp1 if x != rn and cc in meth[x]]
                v2 = [meth[x][cc] for x in hp2 if x != rn and cc in meth[x]]
                if v1:
                    c1[cc] = np.mean(v1)
                if v2:
                    c2[cc] = np.mean(v2)
            shared1 = [cc for cc in pr if cc in c1]; shared2 = [cc for cc in pr if cc in c2]
            if len(shared1) < ANCHOR_MIN or len(shared2) < ANCHOR_MIN:
                continue
            d1 = np.mean([abs(pr[cc] - c1[cc]) for cc in shared1])
            d2 = np.mean([abs(pr[cc] - c2[cc]) for cc in shared2])
            pred = 1 if d1 <= d2 else 2
            rec["loo_tested"] += 1
            rec["loo_correct"] += int(pred == true_hp)
        if rec["loo_tested"] > 0:
            rec["loo_acc"] = round(rec["loo_correct"] / rec["loo_tested"], 3)

        # 完整 HP1/HP2 質心 (給 hp3/unphase 招募用, 非 LOO)
        C1 = {}; C2 = {}
        for cc in anc:
            v1 = [meth[x][cc] for x in hp1 if cc in meth[x]]
            v2 = [meth[x][cc] for x in hp2 if cc in meth[x]]
            if v1:
                C1[cc] = np.mean(v1)
            if v2:
                C2[cc] = np.mean(v2)

        def recruit_assign(pool):
            """回 (recruitable, recruitable_alt, assign_hp1, assign_hp2)。"""
            recn = 0; recalt = 0; a1 = 0; a2 = 0
            for rn in pool:
                pr = prof(rn)
                if len(pr) < ANCHOR_MIN:
                    continue
                recn += 1
                if hasalt.get(rn):
                    recalt += 1
                s1 = [cc for cc in pr if cc in C1]; s2 = [cc for cc in pr if cc in C2]
                if len(s1) < ANCHOR_MIN or len(s2) < ANCHOR_MIN:
                    continue
                d1 = np.mean([abs(pr[cc] - C1[cc]) for cc in s1])
                d2 = np.mean([abs(pr[cc] - C2[cc]) for cc in s2])
                if d1 <= d2:
                    a1 += 1
                else:
                    a2 += 1
            return recn, recalt, a1, a2

        rec["recruitable_hp3"], rec["recruitable_hp3_alt"], rec["hp3_assign_hp1"], rec["hp3_assign_hp2"] = recruit_assign(hp3)
        rec["recruitable_unphase"], rec["recruitable_unphase_alt"], rec["unphase_assign_hp1"], rec["unphase_assign_hp2"] = recruit_assign(unp)

        # (3) within-HP subclone: hp1,hp2 (原始) + hp3 (新)
        for hpg, hpreads in [(1, hp1), (2, hp2)]:
            sub = defaultdict(list)
            for rn in hpreads:
                if rn in geno:
                    sub[geno[rn]].append(rn)
            subs = [ids for ids in sub.values() if len(ids) >= 3]
            if len(subs) >= 2:
                rec["within_hp_subclone_tested"] += 1
                if anchor_cpgs(subs[0], subs[1], meth):
                    rec["within_hp_subclone_sep"] += 1
        sub3 = defaultdict(list)
        for rn in hp3:
            if rn in geno:
                sub3[geno[rn]].append(rn)
        subs3 = [ids for ids in sub3.values() if len(ids) >= 3]
        if len(subs3) >= 2:
            rec["within_hp3_subclone_tested"] += 1
            if anchor_cpgs(subs3[0], subs3[1], meth):
                rec["within_hp3_subclone_sep"] += 1

        records.append(rec)
        if idx % 100 == 0:
            print(f"[{chroms[0] if chroms!=['ALL'] else 'ALL'}] {idx}/{len(tgt)} "
                  f"anchored={sum(1 for x in records if x['n_anchor']>0)}", flush=True)
    tb.close()

    anchored = [r for r in records if r["n_anchor"] > 0]
    loo = [r for r in anchored if r["loo_tested"] >= 5]
    tot_t = sum(r["loo_tested"] for r in loo); tot_c = sum(r["loo_correct"] for r in loo)

    def S(key, rows=anchored):
        return sum(r[key] for r in rows)

    out = {
        "design": "V2: 拆 hp0 → hp3(過ALT) vs unphase(REF-only), 加 ALT-carriage 軸。回答 unphase/HP3 分樹決策。",
        "params": {"hp_min": HP_MIN, "anchor_min": ANCHOR_MIN, "dbeta_min": DBETA_MIN},
        "n_regions": len(records),
        "n_anchored_regions": len(anchored),
        "hp_tag_census_filtered_reads": dict(hp_census),
        "hp_recruitment_HP1HP2_LOO": {
            "regions_with_LOO(>=5 reads)": len(loo),
            "total_reads_classified": tot_t,
            "total_correct": tot_c,
            "overall_accuracy": round(tot_c / tot_t, 4) if tot_t else None,
            "per_region_acc_median": round(float(np.median([r["loo_acc"] for r in loo if r["loo_acc"] != ""])), 3) if loo else None,
        },
        "recruitable_pool_SPLIT": {
            "unphase": {
                "total_in_anchored": S("unphase_n"),
                "recruitable": S("recruitable_unphase"),
                "recruitable_carrying_ALT": S("recruitable_unphase_alt"),
                "assigned_hp1": S("unphase_assign_hp1"),
                "assigned_hp2": S("unphase_assign_hp2"),
            },
            "hp3": {
                "total_in_anchored": S("hp3_n"),
                "recruitable": S("recruitable_hp3"),
                "recruitable_carrying_ALT": S("recruitable_hp3_alt"),
                "assigned_hp1": S("hp3_assign_hp1"),
                "assigned_hp2": S("hp3_assign_hp2"),
            },
            "combined_recruitable(orig_number)": S("recruitable_hp3") + S("recruitable_unphase"),
        },
        "alt_carriage_all_anchored": {
            "hp1": {"n": S("hp1_n"), "alt": S("hp1_alt_n")},
            "hp2": {"n": S("hp2_n"), "alt": S("hp2_alt_n")},
            "hp3": {"n": S("hp3_n"), "alt": S("hp3_alt_n")},
            "unphase": {"n": S("unphase_n"), "alt": S("unphase_alt_n")},
        },
        "within_hp_subclone_HP1HP2": {
            "regions_tested": S("within_hp_subclone_tested", records),
            "regions_methyl_separable": S("within_hp_subclone_sep", records),
        },
        "within_hp3_subclone": {
            "regions_tested": S("within_hp3_subclone_tested", records),
            "regions_methyl_separable": S("within_hp3_subclone_sep", records),
        },
    }
    # rates
    up = out["recruitable_pool_SPLIT"]["unphase"]; h3 = out["recruitable_pool_SPLIT"]["hp3"]
    ac = out["alt_carriage_all_anchored"]
    for k, d in ac.items():
        d["alt_rate"] = round(d["alt"] / d["n"], 4) if d["n"] else None
    ws = out["within_hp_subclone_HP1HP2"]
    ws["separable_frac"] = round(ws["regions_methyl_separable"] / ws["regions_tested"], 4) if ws["regions_tested"] else None
    json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=1)
    tsv = out_path.replace(".json", "_perregion.tsv")
    if records:
        cols = list(records[0].keys())
        with open(tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for rc in records:
                f.write("\t".join(str(rc[c]) for c in cols) + "\n")
    print("DONE anchored=%d HP_LOO_acc=%s | UNPHASE recr=%d alt=%d | HP3 recr=%d alt=%d -> %s" % (
        len(anchored), out["hp_recruitment_HP1HP2_LOO"]["overall_accuracy"],
        up["recruitable"], up["recruitable_carrying_ALT"],
        h3["recruitable"], h3["recruitable_carrying_ALT"], out_path), flush=True)
    print("ALT-carriage rates: hp1=%s hp2=%s hp3=%s unphase=%s" % (
        ac["hp1"]["alt_rate"], ac["hp2"]["alt_rate"], ac["hp3"]["alt_rate"], ac["unphase"]["alt_rate"]), flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])

