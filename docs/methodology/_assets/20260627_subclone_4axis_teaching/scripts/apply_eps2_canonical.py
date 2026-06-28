#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ε=2% 定案：產出 canonical 可逐筆驗證輸出。
規則：cell 為「present(真)」當 count > coread × 0.02（ONT 噪聲底線，隨深度縮放）。
  保留最低 1 條（低 coread 時 1 > floor 仍算）；高 coread 單讀（1 ≤ floor）判 noise。
輸出：
  pairs_eps2_annotated.tsv  每對附 floor / 強弱 / orig→new config（可手驗）
  eps2_final_band.json      orig vs ε2% config band + 重分類統計
用法：python3 apply_eps2_canonical.py
"""
import os, csv, json
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
EPS = 0.02

rows = []
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if os.path.exists(p):
            rows += list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))

def classify_flags(aaP, raP, arP):
    if not aaP: return "mutual_excl" if (raP or arP) else "sparse"
    if not raP and not arP: return "co_linked"
    if not arP: return "nested_a_in_b"
    if not raP: return "nested_b_in_a"
    return "independent"

orig = Counter(); new = Counter()
out_rows = []
n_changed = 0; n_weak_kept = 0; n_noise = 0
for d in rows:
    rr, ra, aa = int(d["RR"]), int(d["RA"]), int(d["AA"]); ar = int(d["AR"]); N = int(d["coread"])
    floor = N * EPS
    aaP, raP, arP = aa > floor, ra > floor, ar > floor
    nc = classify_flags(aaP, raP, arP)
    oc = d["config"]
    orig[oc] += 1; new[nc] += 1
    # 決定格(off-diagonal)讀數 + 強弱
    if oc == "nested_a_in_b": dec = ra
    elif oc == "nested_b_in_a": dec = ar
    elif oc == "independent": dec = min(ra, ar)
    else: dec = None
    conf = "NA"
    if dec is not None:
        if dec >= 2: conf = "strong"
        elif dec == 1 and dec > floor: conf = "weak_kept"; n_weak_kept += 1
        else: conf = "noise"; n_noise += 1
    if nc != oc: n_changed += 1
    out_rows.append({"chrom": d["chrom"], "pos_a": d["pos_a"], "pos_b": d["pos_b"], "dist": d["dist"],
                     "RR": rr, "RA": ra, "AR": ar, "AA": aa, "coread": N, "floor_2pct": round(floor, 2),
                     "orig_config": oc, "eps2_config": nc, "changed": int(nc != oc),
                     "deciding_off": "" if dec is None else dec, "confidence": conf,
                     "src_a": d.get("src_a", ""), "src_b": d.get("src_b", ""), "cn_a": d.get("cn_a", "")})

with open(os.path.join(DATA, "pairs_eps2_annotated.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out_rows)

def merge_nested(c): return {"co_linked": c["co_linked"], "nested": c["nested_a_in_b"]+c["nested_b_in_a"],
                             "independent": c["independent"], "mutual_excl": c["mutual_excl"], "sparse": c.get("sparse", 0)}
band = {"eps": EPS, "n_pairs": len(rows),
        "orig_config": dict(orig), "eps2_config": dict(new),
        "orig_merged": merge_nested(orig), "eps2_merged": merge_nested(new),
        "n_reclassified": n_changed, "pct_reclassified": round(100*n_changed/len(rows), 1),
        "weak_single_read_kept(coread<50)": n_weak_kept, "single_read_to_noise(coread>=50)": n_noise}
with open(os.path.join(DATA, "eps2_final_band.json"), "w", encoding="utf-8") as f:
    json.dump(band, f, ensure_ascii=False, indent=1)

print(f"ε=2% 套用 {len(rows):,} pairs")
print(f"  orig  merged: {merge_nested(orig)}")
print(f"  eps2  merged: {merge_nested(new)}")
print(f"  重分類 {n_changed:,} ({round(100*n_changed/len(rows),1)}%); 單讀保留(coread<50) {n_weak_kept:,}; 單讀判噪聲(coread>=50) {n_noise:,}")
print(f"  → pairs_eps2_annotated.tsv ({len(out_rows):,} 行) + eps2_final_band.json")
