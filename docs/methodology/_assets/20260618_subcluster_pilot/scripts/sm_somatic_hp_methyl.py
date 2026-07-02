"""[somatic HP 子相位甲基印證] 回答使用者 Q2: 先前 HP-axis 把 1-1 折進 HP1(germline 軸)。
改測 longphase-S HP tag **內含的 somatic 子相位**: 同一 germline HP 內,
  germline-only 讀(tag "1") vs somatic-sub-phase 讀(tag "1-1") 甲基是否不同?
  (HP2 同理: "2" vs "2-1")
🔑 這是**同一 germline 單倍型內**的 somatic 軸 → 非 germline-ASM 循環; 甲基若對齊 = 印證 longphase 的 somatic haplotagging。
⚠ 判讀: 1-1 由帶 somatic ALT 定義 → 甲基差異仍混 somatic-cis(突變直接效應)+subclone,需 normal 分離; 但「甲基獨立支持 somatic 相位」本身是有效觀察。
單位點隱藏 subclone(Q1): 也記錄每區 within-somatic-HP 是否有甲基結構。
argv: chroms out_path。HP 用 raw string(pysam 對 HP:Z 回 str)。
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
GRP_MIN = 3      # 每組(pure/som)至少幾條 read
CPG_MIN = 3


def bh(p):
    p = np.asarray(p); n = len(p)
    if n == 0:
        return np.array([])
    o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]; prev = min(prev, p[i] * n / (rank + 1)); q[i] = prev
    return q


def hp_norm(v):
    """raw HP → 標準字串 (兼容 longphase-to 整數)。"""
    s = str(v)
    return {"11": "1-1", "21": "2-1", "33": "3"}.get(s, s)


def read_meth_hp(a):
    try:
        mb = a.modified_bases
    except Exception:
        return None
    pairs = a.get_aligned_pairs(matches_only=True)
    qr = {q: r for q, r in pairs}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != 'm':
                continue
            for qpos, mlq in lst:
                r = qr.get(qpos)
                if r is not None:
                    meth[r] = mlq / 255.0
    hp = a.get_tag("HP") if a.has_tag("HP") else None
    return meth, hp_norm(hp) if hp is not None else "."


def per_cpg(g1, g2, meth):
    cps = set()
    for rn in g1 + g2:
        cps |= set(meth[rn].keys())
    pv = []; db = []
    for cp in cps:
        a_ = [meth[rn][cp] for rn in g1 if cp in meth[rn]]
        b_ = [meth[rn][cp] for rn in g2 if cp in meth[rn]]
        if len(a_) >= CPG_MIN and len(b_) >= CPG_MIN:
            try:
                _, p = mannwhitneyu(a_, b_, alternative="two-sided")
            except ValueError:
                continue
            pv.append(p); db.append(abs(np.mean(a_) - np.mean(b_)))
    if not pv:
        return 0, 0, 0.0
    q = bh(pv)
    nsig = sum(1 for k in range(len(pv)) if q[k] < 0.05 and db[k] >= DBETA_MIN)
    return len(pv), nsig, float(max(db)) if db else 0.0


def main(chroms, out_path):
    reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]
    tb = pysam.AlignmentFile(TBAM, "rb")
    chset = set(chroms)
    tgt = [r for r in reg if (chroms == ["ALL"] or r["chrom"] in chset)
           and r["cn"] in ("loh", "neutral")
           and r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")
           and (r.get("n_populations") or 0) >= 2]
    records = []
    for idx, r in enumerate(tgt):
        chrom = r["chrom"]
        meth = {}; hp = {}
        for a in tb.fetch(chrom, r["start"], r["end"] + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ_MIN:
                continue
            res = read_meth_hp(a)
            if res is None:
                continue
            m, h = res
            meth[a.query_name] = m; hp[a.query_name] = h
        g = defaultdict(list)
        for rn, h in hp.items():
            g[h].append(rn)
        rec = {"region": r["region"], "cn": r["cn"], "shape": r["tree_shape"],
               "n_pure1": len(g["1"]), "n_som1": len(g["1-1"]),
               "n_pure2": len(g["2"]), "n_som2": len(g["2-1"]),
               "hp1_eval": 0, "hp1_nsig": 0, "hp1_dbeta": 0.0, "hp1_aligned": 0,
               "hp2_eval": 0, "hp2_nsig": 0, "hp2_dbeta": 0.0, "hp2_aligned": 0,
               "any_aligned": 0}
        # HP1: pure "1" vs somatic "1-1"
        if len(g["1"]) >= GRP_MIN and len(g["1-1"]) >= GRP_MIN:
            rec["hp1_eval"] = 1
            _, ns, mdb = per_cpg(g["1"], g["1-1"], meth)
            rec["hp1_nsig"], rec["hp1_dbeta"] = ns, round(mdb, 3)
            rec["hp1_aligned"] = 1 if ns > 0 else 0
        # HP2: pure "2" vs somatic "2-1"
        if len(g["2"]) >= GRP_MIN and len(g["2-1"]) >= GRP_MIN:
            rec["hp2_eval"] = 1
            _, ns, mdb = per_cpg(g["2"], g["2-1"], meth)
            rec["hp2_nsig"], rec["hp2_dbeta"] = ns, round(mdb, 3)
            rec["hp2_aligned"] = 1 if ns > 0 else 0
        rec["any_aligned"] = 1 if (rec["hp1_aligned"] or rec["hp2_aligned"]) else 0
        records.append(rec)
        if idx % 100 == 0:
            print(f"[{chroms[0] if chroms!=['ALL'] else 'ALL'}] {idx}/{len(tgt)} "
                  f"eval={sum(1 for x in records if x['hp1_eval'] or x['hp2_eval'])} "
                  f"aligned={sum(x['any_aligned'] for x in records)}", flush=True)
    tb.close()

    evalr = [r for r in records if r["hp1_eval"] or r["hp2_eval"]]
    aligned = [r for r in evalr if r["any_aligned"]]
    # per-germline-HP 統計(把 HP1/HP2 各算一個測試單元)
    hp_units_eval = sum(r["hp1_eval"] + r["hp2_eval"] for r in records)
    hp_units_aligned = sum(r["hp1_aligned"] + r["hp2_aligned"] for r in records)
    dbetas = [r["hp1_dbeta"] for r in records if r["hp1_aligned"]] + [r["hp2_dbeta"] for r in records if r["hp2_aligned"]]
    out = {
        "question": "甲基是否對齊 longphase-S somatic HP 子相位(同 germline HP 內 1 vs 1-1)? = 甲基印證 somatic haplotagging",
        "params": {"grp_min": GRP_MIN, "cpg_min": CPG_MIN, "dbeta_min": DBETA_MIN},
        "note": "1-1 由帶 somatic ALT 定義 → 對齊仍混 somatic-cis + subclone, 需 normal 分離; 但非 germline-ASM 循環",
        "n_regions": len(records),
        "region_level": {
            "evaluable_regions": len(evalr),
            "aligned_regions": len(aligned),
            "aligned_frac": round(len(aligned) / len(evalr), 4) if evalr else None,
            "by_cn": {"loh": {"eval": sum(1 for r in evalr if r["cn"] == "loh"),
                              "aligned": sum(1 for r in aligned if r["cn"] == "loh")},
                      "neutral": {"eval": sum(1 for r in evalr if r["cn"] == "neutral"),
                                  "aligned": sum(1 for r in aligned if r["cn"] == "neutral")}},
        },
        "hp_unit_level": {
            "evaluable_hp_units": hp_units_eval,
            "aligned_hp_units": hp_units_aligned,
            "aligned_frac": round(hp_units_aligned / hp_units_eval, 4) if hp_units_eval else None,
        },
        "effect": {
            "dbeta_median": round(float(np.median(dbetas)), 3) if dbetas else None,
            "dbeta_max": round(float(np.max(dbetas)), 3) if dbetas else None,
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
    rl = out["region_level"]
    print(f"DONE: eval_regions={rl['evaluable_regions']} aligned={rl['aligned_regions']} "
          f"({rl['aligned_frac']}) hp_units={out['hp_unit_level']['aligned_frac']} "
          f"dbeta_med={out['effect']['dbeta_median']} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
