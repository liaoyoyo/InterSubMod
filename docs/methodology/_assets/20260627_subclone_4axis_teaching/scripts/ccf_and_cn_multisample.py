#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""6 樣本 CCF read 加權 + CN 套用（2026-07-10，補 item a+b；不重讀 BAM）。

- CCF：對每樣本 ambiguous 單位用 mlhp col_coverage_by_hp 算 per-mutation CCF → pigeonhole 加權（reuse ccf_full_observe 邏輯）。
- CN 套用（P0-4 修）：CN 不從 mlhp 的 neutral 預設讀，改**CN bed lookup 覆蓋**：
    HCC1395/HCC1395_DORADO → SEQC2 CN（同細胞株）；H1437/H2009/HCC1954 → SAVANA CN（purity 自洽·memory 定案可用）；
    COLO829/HCC1937 → **unavailable**（SAVANA purity mis-fit·memory 定案 cn=unknown）。
- 誠實拆 CN-clean：strictly-neutral（真 read≈CCF）vs loh vs loss（read 變少·相反 confound）vs gain（read≠CCF）vs unavailable。
- recurrence CN 裁決：CN-available 樣本可判 artifact(gain)/LOH/candidate(neutral)；unavailable 樣本標 unavailable（非 neutral 預設）。
只讀不寫源碼。輸出 ccf_and_cn_multisample.json。用法：python3 ccf_and_cn_multisample.py
"""
import json, os, math, glob, bisect
from collections import Counter, defaultdict
import tree_enumeration_solver as S

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
CNVWORK = "/big7_disk/liaoyoyo2001/cnv_sv_work"
SEQC2 = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
MARGIN = 0.05; TEMP = 0.05; GERM = ("1", "2", "3")

# (sample, mlhp_glob, layered_json, cn_bed 或 None, cn_source)
SAMPLES = [
    ("HCC1395",        f"{PILOT}/mlhp_part_*.json",             f"{PILOT}/layered_reconstruction_HCC1395.json",             SEQC2, "SEQC2"),
    ("HCC1395_DORADO", f"{MSROOT}/HCC1395_DORADO/mlhp_part_*.json", f"{MSROOT}/HCC1395_DORADO/layered_reconstruction_HCC1395_DORADO.json", SEQC2, "SEQC2(同細胞株)"),
    ("H1437",          f"{MSROOT}/H1437/mlhp_part_*.json",      f"{MSROOT}/H1437/layered_reconstruction_H1437.json",        f"{CNVWORK}/H1437/savana_wgs/H1437_smcnbed.bed", "SAVANA"),
    ("H2009",          f"{MSROOT}/H2009/mlhp_part_*.json",      f"{MSROOT}/H2009/layered_reconstruction_H2009.json",        f"{CNVWORK}/H2009/savana_wgs/H2009_smcnbed.bed", "SAVANA"),
    ("HCC1954",        f"{MSROOT}/HCC1954/mlhp_part_*.json",    f"{MSROOT}/HCC1954/layered_reconstruction_HCC1954.json",    f"{CNVWORK}/HCC1954/savana_wgs/HCC1954_smcnbed.bed", "SAVANA"),
    ("COLO829",        f"{MSROOT}/COLO829/mlhp_part_*.json",    f"{MSROOT}/COLO829/layered_reconstruction_COLO829.json",    None, "unavailable(SAVANA mis-fit)"),
    ("HCC1937",        f"{MSROOT}/HCC1937/mlhp_part_*.json",    f"{MSROOT}/HCC1937/layered_reconstruction_HCC1937.json",    None, "unavailable(SAVANA mis-fit)"),
]


def load_cn(bed):
    if not bed or not os.path.exists(bed): return None
    iv = defaultdict(lambda: ([], []))  # chrom → (starts sorted, [(end,state)])
    tmp = defaultdict(list)
    for ln in open(bed):
        p = ln.split()
        if len(p) < 4: continue
        tmp[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c, lst in tmp.items():
        lst.sort(); iv[c] = ([s for s, _, _ in lst], [(e, st) for _, e, st in lst])
    return iv


def cn_at(cn, chrom, pos):
    """CN bed lookup：命中 gain/loss/loh 回該 state；未命中=neutral（bed 只列非中性）。無 bed→unavailable。"""
    if cn is None: return "unavailable"
    if chrom not in cn: return "neutral"
    starts, ends = cn[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and pos <= ends[i][0]: return ends[i][1]
    return "neutral"


def ccf_from_colcov(colcov, positions, k):
    if not colcov or len(positions) < k: return None
    ccf = {}
    for j in range(k):
        v = colcov.get(str(positions[j])) or colcov.get(positions[j])
        if not v: return None
        nr, na = v[0], v[1]; tot = nr + na
        ccf[j] = na / tot if tot else 0.0
    return ccf


def _unlab(s):
    if s == "ROOT": return frozenset()
    s2 = s[2:] if s.startswith("H_") else s
    return frozenset(i for i, c in enumerate(s2) if c == "A")


def score_viol(edges_lab, ccf):
    sc = 0.0; vi = 0
    for (p, c) in edges_lab:
        P, C = _unlab(p), _unlab(c); acq = C - P
        if len(acq) != 1: continue
        j = next(iter(acq))
        for anc in P:
            sc += (ccf[anc] - ccf[j])
            if ccf[j] > ccf[anc] + MARGIN: vi += 1
    return sc, vi


def cn_class_of(state):
    if state == "gain": return "gain(read≠CCF)"
    if state == "neutral": return "strictly-neutral(read≈CCF·可信)"
    if state == "loh": return "loh"
    if state == "loss": return "loss(read變少)"
    return "unavailable"


def run_sample(name, mlglob, layered, cn):
    """用 layered JSON 已存 trees（免重枚舉）+ join mlhp col_coverage 算 CCF。"""
    # mlhp region 索引：region_key → (col_coverage_by_hp, positions)
    mreg = {}
    for f in sorted(glob.glob(mlglob)):
        for g in json.load(open(f)).get("groups", []):
            mreg[f"{g['chrom']}:{g['start']}-{g['end']}"] = g
    detail = json.load(open(layered))["detail"]
    amb_units = [u for u in detail if u.get("is_lineage") and "ambiguous" in u["L1_class"]
                 and u.get("trees") and len(u["trees"]) >= 2]
    n_amb = broke = win_clean = 0
    by = Counter(); by_broke = Counter()
    for u in amb_units:
        g = mreg.get(u["region"])
        if not g: continue
        positions = g.get("positions", []); fam = u["family"]
        k = len(positions)
        if k < 2 or k > 8: continue
        mid = positions[len(positions) // 2]
        cnc = cn_class_of(cn_at(cn, u["chrom"], mid))
        cbh = (g.get("col_coverage_by_hp", {}) or {}).get(fam, {})
        ccf = ccf_from_colcov(cbh, positions, k)
        if ccf is None:
            full = (g.get("populations_by_hp", {}) or {}).get(fam, {}) or {}
            if not full: continue
            tot = sum(full.values()) or 1
            ccf = {j: sum(n for gt, n in full.items() if gt[j] == "A") / tot for j in range(k)}
        n_amb += 1
        scored = [score_viol(t["edges"], ccf) + (t,) for t in u["trees"]]   # 已存 trees（≤32）
        mx = max(s for s, _, _ in scored)
        exps = [math.exp((s - mx) / TEMP) for s, _, _ in scored]; tot = sum(exps) or 1
        post = [e / tot for e in exps]; top = max(post)
        by[cnc] += 1
        if top >= 0.6:
            broke += 1; by_broke[cnc] += 1
            if scored[post.index(top)][1] == 0: win_clean += 1
    strict = by.get("strictly-neutral(read≈CCF·可信)", 0)
    strict_broke = by_broke.get("strictly-neutral(read≈CCF·可信)", 0)
    return {
        "n_ambiguous": n_amb, "tie_broken": broke, "broke_pct": round(broke / n_amb * 100, 1) if n_amb else 0,
        "winner_selfconsistent": win_clean, "winner_selfconsist_pct_of_broke": round(win_clean / broke * 100, 1) if broke else None,
        "cn_breakdown_total": dict(by), "cn_breakdown_broke": dict(by_broke),
        "strictly_neutral_trustworthy": strict, "strictly_neutral_broke": strict_broke,
        "trustworthy_reach_pct": round(strict_broke / n_amb * 100, 1) if n_amb else 0,
        "note": "winner_selfconsist=自洽率非獨立驗證；trustworthy_reach 只算 strictly-neutral（read≈CCF）",
    }


def recurrence_cn(layered, cn):
    detail = json.load(open(layered))["detail"]
    rec = [u for u in detail if u.get("is_lineage") and u["L1_class"] == "recurrence_required"]
    verdict = Counter()
    for u in rec:
        pos = (u["start"] + u["end"]) // 2
        st = cn_at(cn, u["chrom"], pos)
        v = {"gain": "artifact(CN-gain)", "loh": "LOH", "loss": "loss", "neutral": "candidate(real?)", "unavailable": "unavailable"}[st]
        verdict[v] += 1
    return {"n_recurrence": len(rec), "cn_verdict": dict(verdict)}


out = {"meta": {"date": "2026-07-10", "scope": "chr1–22",
                "cn_coverage": "HCC1395/DORADO=SEQC2; H1437/H2009/HCC1954=SAVANA(可用); COLO829/HCC1937=unavailable(mis-fit)",
                "caveats": "winner 自洽率非驗證；trustworthy reach=strictly-neutral only；CCF=未經purity校正 ALT fraction"},
       "samples": {}}
for name, mlglob, layered, bed, src in SAMPLES:
    cn = load_cn(bed)
    out["samples"][name] = {"cn_source": src, "ccf": run_sample(name, mlglob, layered, cn), "recurrence": recurrence_cn(layered, cn)}
json.dump(out, open(f"{D}/ccf_and_cn_multisample.json", "w"), ensure_ascii=False, indent=1)

print(f"{'sample':16} {'cn_src':22} {'amb':>6} {'broke':>6} {'broke%':>6} {'strictNeu':>9} {'trust%':>6} {'recCN'}")
for name, *_ , src in SAMPLES:
    c = out["samples"][name]["ccf"]; rc = out["samples"][name]["recurrence"]
    print(f"{name:16} {src[:22]:22} {c['n_ambiguous']:>6} {c['tie_broken']:>6} {c['broke_pct']:>6} "
          f"{c['strictly_neutral_trustworthy']:>9} {c['trustworthy_reach_pct']:>6}  {rc['cn_verdict']}")
