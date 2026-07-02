#!/usr/bin/env python3
"""[P1 候選池可測性 + confound profile] 1139 候選 join records_v6 → normal 覆蓋(可測性)、
CN/LOH 曝露(confound)、結構數、既有 Δβ。決定 cis-test 分母 + confound 暴露。純讀。"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
recs = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
idx = {f"{r['chrom']}:{r['pos']}": r for r in recs}
cands = json.load(open(f"{A}/cis_candidates_resolved.json"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


prof = []
for c in cands:
    r = idx.get(f"{c['chrom']}:{c['pos']}")
    if not r:
        continue
    prof.append({"key": f"{c['chrom']}:{c['pos']}", "set": c["set"], "cat8": c.get("cat8"), "pc_verdict": c.get("pc_verdict"),
                 "n_tumor": r.get("n_tumor"), "n_normal": r.get("n_normal"), "is_loh": r.get("is_loh"),
                 "cn_state": r.get("cn_state"), "cn_value": r.get("cn_value"), "coarse_ng": r.get("coarse_ng"),
                 "m_dbeta_tn": num(r.get("m_dbeta_tn")), "m_dbeta_hp": num(r.get("m_dbeta_hp")),
                 "m_dbeta_group": num(r.get("m_dbeta_group")), "pc_cluster_run": r.get("pc_cluster_run"),
                 "pc_overlap_max": num(r.get("pc_overlap_max")), "ls_hp_p": num(r.get("ls_hp_p"))})


def cov_buckets(vals):
    v = [x for x in vals if x is not None]
    return {"n": len(v), ">=5": sum(x >= 5 for x in v), ">=10": sum(x >= 10 for x in v),
            ">=20": sum(x >= 20 for x in v), "median": sorted(v)[len(v) // 2] if v else None}


def summ(items, name):
    nn = [p["n_normal"] for p in items]
    nt = [p["n_tumor"] for p in items]
    return {"n": len(items), "normal_cov": cov_buckets(nn), "tumor_cov": cov_buckets(nt),
            "is_loh_true": sum(1 for p in items if p["is_loh"] in (True, "true", "True", 1)),
            "cn_state": dict(Counter(str(p["cn_state"]) for p in items).most_common()),
            "coarse_ng": dict(Counter(str(p["coarse_ng"]) for p in items).most_common(8)),
            "m_dbeta_group_gt02": sum(1 for p in items if (p["m_dbeta_group"] or 0) >= 0.2)}


out = {"total": len(prof),
       "ALL": summ(prof, "ALL"),
       "TP": summ([p for p in prof if p["set"] == "TP"], "TP"),
       "FP": summ([p for p in prof if p["set"] == "FP"], "FP"),
       "by_cat8": {k: summ([p for p in prof if p["cat8"] == k], k) for k in ("A", "C-S1", "C-S2", "C-S3", "D", "B1")},
       "subclone_novel_only": summ([p for p in prof if p["pc_verdict"] == "subclone_novel"], "sn")}
json.dump({"summary": out, "per_locus": prof}, open(f"{A}/candidate_characterization.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
print("\n[-> candidate_characterization.json]")
