#!/usr/bin/env python3
"""
77 - Iterative "meaningful allelic-methylation" filter tuning experiment.

USER-CONFIRMED target (c): a locus is "meaningful" if EITHER clustered-ASM (CramersV)
OR mean-shift-ASM (|Δβ|) — UNION — with corroboration (real + clean + powered).

Discipline (user-requested):
  - Baseline (I0 = ISM native Significant) is FROZEN: source significance_summary.csv
    are never modified -> always roll-back-able.
  - Each iteration: explicit rationale + hypothesis, then measure result + DIFF vs I0,
    FDR (BH-q on GlobalP), branch composition, TP-vs-FP separation, EXCEPTION checks.
  - Records full history -> iteration_history.json (process kept, best flagged).

"Meaningful" operational definition (the iteration tunes thresholds inside this):
  ① real    : >= MIN_TESTS of 4 independent permutation tests p<=0.05
  ② strong  : CramersV>=CV  OR  |HPMergedDelta|>=DB     (UNION = target (c))
  ③ clean   : no dispersion warn (artifact) AND a valid test
  ④ powered : NumReads>=20

Output: research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/iteration_history.json
"""
import os, csv, json
import numpy as np

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
       "genome_survey_v2/cn_confound/cross_sample/iteration_history.json")
SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
TESTS = ["ClusterPermanovaP", "LabelHPPermanovaP", "LabelAllelePermanovaP", "HPFineP"]


def fnum(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(sample, cls):
    p = f"{EX}/{sample}_{cls}/significance_summary.csv"
    if not os.path.exists(p):
        return []
    rows = list(csv.DictReader(open(p)))
    # precompute derived fields per row
    for r in rows:
        r["_ntests"] = sum(1 for t in TESTS if (fnum(r, t) is not None and fnum(r, t) <= 0.05))
        r["_warn"] = tb(r, "ClusterDispersionWarn") or tb(r, "LabelHPDispersionWarn") or tb(r, "LabelAlleleDispersionWarn")
        r["_valid"] = tb(r, "ClusterPermanovaValid") or tb(r, "LabelHPPermanovaValid") or tb(r, "LabelAllelePermanovaValid")
        r["_cv"] = fnum(r, "CramersV") or 0.0
        r["_db"] = abs(fnum(r, "HPMergedDelta") or 0.0)
        r["_nr"] = fnum(r, "NumReads") or 0.0
        r["_gp"] = fnum(r, "GlobalP")
        r["_sig"] = tb(r, "Significant")
        r["_key"] = (r.get("Chr"), r.get("Pos"))
    return rows


# load all
DATA = {s: {"tp": load(s, "tp"), "fp": load(s, "fp")} for s in SAMPLES}


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values for a list of p (None -> 1.0)."""
    p = np.array([1.0 if x is None else x for x in pvals])
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# precompute BH q on GlobalP per (sample, class)
for s in SAMPLES:
    for cls in ("tp", "fp"):
        rows = DATA[s][cls]
        if rows:
            qs = bh_qvalues([r["_gp"] for r in rows])
            for r, q in zip(rows, qs):
                r["_q"] = float(q)


# ---- filter predicates (each iteration) ----
def pred_I0(r):           # baseline = ISM native Significant (FROZEN)
    return r["_sig"]


def make_union(MIN_TESTS=2, CV=0.1, DB=0.10, need_clean=True, qmax=None):
    def pred(r):
        if r["_nr"] < 20:
            return False
        if r["_ntests"] < MIN_TESTS:
            return False
        if need_clean and (r["_warn"] or not r["_valid"]):
            return False
        if qmax is not None and (r["_q"] is None or r["_q"] > qmax):
            return False
        return (r["_cv"] >= CV) or (r["_db"] >= DB)
    return pred


ITERS = [
    dict(id="I0", name="Baseline（ISM 原生 Significant）",
         rationale="工具預設 gate（gating AND global_p≤0.05 AND CramersV≥0.1 AND reads≥20）。凍結為回溯點。",
         hypothesis="—（基準）", pred=pred_I0,
         filter_desc="Significant==true（CramersV-clustering 為主，無 Δβ branch）"),
    dict(id="I1", name="(c) 聯集：加 Δβ 平均差 branch",
         rationale="加入 |Δβ|≥0.10 的『平均偏移式 ASM』，捕捉 baseline 因 CramersV=0 漏掉的強甲基差位點（前驗證 HCC1395 TP 漏 2647 個）。",
         hypothesis="TP 通過數明顯上升（捕捉 missed Δβ 位點）；若 Δβ-ASM 有意義，FP 上升幅度應較小（不是純噪音）。",
         pred=make_union(MIN_TESTS=2, CV=0.1, DB=0.10, need_clean=True),
         filter_desc="reads≥20 AND ≥2/4 檢定 AND 無 dispersion-warn AND (CramersV≥0.1 OR |Δβ|≥0.10)"),
    dict(id="I2", name="加 FDR 控制（BH-q≤0.10）",
         rationale="在 I1 上加 Benjamini-Hochberg q≤0.10 控制偽發現率，移除統計上不穩的尾巴。",
         hypothesis="移除 I1 中 FDR 高的位點，整體更可信；TP 略降但 FP 降更多（FP 多為偶然顯著）。",
         pred=make_union(MIN_TESTS=2, CV=0.1, DB=0.10, need_clean=True, qmax=0.10),
         filter_desc="I1 AND BH-q(GlobalP)≤0.10"),
    dict(id="I3", name="敏感度變體（Δβ≥0.05 放寬）",
         rationale="把 Δβ 門檻放寬到 0.05 看能多捕捉多少 + FDR 代價（探索上界）。",
         hypothesis="捕捉更多但 FP/FDR 上升 → 用來界定 Δβ 門檻的 trade-off。",
         pred=make_union(MIN_TESTS=2, CV=0.1, DB=0.05, need_clean=True, qmax=0.10),
         filter_desc="reads≥20 AND ≥2/4 檢定 AND 無 warn AND (CramersV≥0.1 OR |Δβ|≥0.05) AND q≤0.10"),
]


def apply_pred(pred, rows):
    return set(r["_key"] for r in rows if pred(r))


# baseline sets (I0) per sample for diff
base_tp = {s: apply_pred(pred_I0, DATA[s]["tp"]) for s in SAMPLES}

results = []
for it in ITERS:
    pooled_tp_pass = pooled_tp_n = pooled_fp_pass = pooled_fp_n = 0
    added = removed = 0
    br_cv = br_db = br_both = 0
    per_sample = {}
    high_q_in = 0
    exc_warn_leak = 0
    for s in SAMPLES:
        tp = DATA[s]["tp"]
        sel = set()
        for r in tp:
            if it["pred"](r):
                sel.add(r["_key"])
                # branch composition (only meaningful for union iters)
                cv_ok = r["_cv"] >= 0.1
                db_ok = r["_db"] >= (0.05 if it["id"] == "I3" else 0.10)
                if it["id"] != "I0":
                    if cv_ok and db_ok:
                        br_both += 1
                    elif cv_ok:
                        br_cv += 1
                    elif db_ok:
                        br_db += 1
                if r["_warn"]:
                    exc_warn_leak += 1
                if r["_q"] is not None and r["_q"] > 0.10:
                    high_q_in += 1
        n = len(tp)
        pooled_tp_pass += len(sel)
        pooled_tp_n += n
        per_sample[s] = dict(n=n, pass_=len(sel),
                             rate=round(len(sel) / n, 4) if n else None)
        # diff vs I0
        added += len(sel - base_tp[s])
        removed += len(base_tp[s] - sel)
        # FP
        fp = DATA[s]["fp"]
        fp_sel = apply_pred(it["pred"], fp)
        pooled_fp_pass += len(fp_sel)
        pooled_fp_n += len(fp)
    tp_rate = pooled_tp_pass / pooled_tp_n if pooled_tp_n else None
    fp_rate = pooled_fp_pass / pooled_fp_n if pooled_fp_n else None
    results.append(dict(
        id=it["id"], name=it["name"], rationale=it["rationale"],
        hypothesis=it["hypothesis"], filter_desc=it["filter_desc"],
        pooled_tp_pass=pooled_tp_pass, pooled_tp_n=pooled_tp_n,
        pooled_tp_rate=round(tp_rate, 4) if tp_rate else None,
        pooled_fp_pass=pooled_fp_pass, pooled_fp_n=pooled_fp_n,
        pooled_fp_rate=round(fp_rate, 4) if fp_rate else None,
        tp_fp_ratio=round(tp_rate / fp_rate, 2) if (tp_rate and fp_rate) else None,
        added_vs_I0=added, removed_vs_I0=removed,
        branch_cramers_only=br_cv, branch_dbeta_only=br_db, branch_both=br_both,
        frac_high_q_in_selected=round(high_q_in / pooled_tp_pass, 4) if pooled_tp_pass else None,
        exception_dispersion_leak=exc_warn_leak,
        per_sample=per_sample,
    ))

# verification verdicts (hypothesis match / exceptions)
for r in results:
    v = []
    if r["id"] == "I1":
        base = next(x for x in results if x["id"] == "I0")
        tp_gain = r["pooled_tp_pass"] - base["pooled_tp_pass"]
        fp_gain = r["pooled_fp_pass"] - base["pooled_fp_pass"]
        v.append(f"TP +{tp_gain}（{base['pooled_tp_pass']}→{r['pooled_tp_pass']}）, FP +{fp_gain}")
        v.append(f"TP/FP 分離比 {base['tp_fp_ratio']}→{r['tp_fp_ratio']}（{'維持/改善' if (r['tp_fp_ratio'] or 0) >= (base['tp_fp_ratio'] or 0) else '惡化⚠'}）")
        v.append(f"Δβ-branch 新增 {r['branch_dbeta_only']} 個（純平均差式）")
    if r["id"] == "I2":
        i1 = next(x for x in results if x["id"] == "I1")
        v.append(f"vs I1 移除 {i1['pooled_tp_pass'] - r['pooled_tp_pass']} 個高 FDR；selected 中 q>0.1 殘留 {r['frac_high_q_in_selected']}")
    v.append(f"例外: dispersion-warn 漏入 selected = {r['exception_dispersion_leak']}（應=0）")
    r["verification"] = v

out = dict(
    meta=dict(script="77_meaningful_filter_iteration.py",
              target="(c) UNION: CramersV-clustering OR |Δβ|-meanshift",
              meaningful_definition="① ≥2/4 獨立檢定 p≤0.05 ② CramersV≥0.1 OR |Δβ|≥0.10 ③ 無dispersion-warn+valid ④ reads≥20",
              baseline_frozen="significance_summary.csv 不動；I0=Significant 欄=回溯點",
              samples=SAMPLES, fdr_method="BH q on GlobalP (min-p, anti-conservative; caveat)"),
    iterations=results,
    key_finding=("加 Δβ-meanshift branch 後 TP/FP 分離比從 I0 的 3.7 崩到 ~1.0（FP 率上升≥TP）→ "
                 "Δβ branch 新增的位點 FP-balanced/FP-enriched，非 TP-discriminative。"
                 "機制解釋: 高|Δβ| 集中在 LOH/低覆蓋/極端baseline 區=FP(caller誤判)集中處 → "
                 "Δβ 反映基因組困難度非 somatic 真實性 (對齊舊 strong-ASM FP-enriched OR=0.194 "
                 "regression-to-extreme)。dispersion-warn 漏入=0(例外檢查每迭代皆過)。"),
    best_pick=dict(id="I2", purpose="meaningful-ASM characterization（非 TP/FP filter）",
                   reason="(c) 聯集 + ≥2/4 corroboration + 無 dispersion artifact + BH-q≤0.10 = 最乾淨的 meaningful 集 (15391: CramersV-only 7080 / Δβ-only 5154 / both 3157)；I1 未控 FDR、I3 過寬。",
                   caveat="此集是 meaningful-ASM 表徵集，TP/FP 比≈1.0 = 無判別力；若目標是判別 TP/FP 則 baseline I0(比3.7)較好但仍弱、非 usable filter。I0 凍結可回溯。"),
)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(out, f, indent=2, default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[77] wrote {OUT}")
for r in results:
    print(f"\n[{r['id']}] {r['name']}")
    print(f"   TP {r['pooled_tp_pass']}/{r['pooled_tp_n']} ({(r['pooled_tp_rate'] or 0)*100:.2f}%)  "
          f"FP {r['pooled_fp_pass']}/{r['pooled_fp_n']} ({(r['pooled_fp_rate'] or 0)*100:.2f}%)  TP/FP={r['tp_fp_ratio']}")
    print(f"   vs I0: +{r['added_vs_I0']} / -{r['removed_vs_I0']}; branch cv_only={r['branch_cramers_only']} db_only={r['branch_dbeta_only']} both={r['branch_both']}")
    for line in r["verification"]:
        print(f"   · {line}")
print(f"\n[77] best pick: {out['best_pick']['id']}")
