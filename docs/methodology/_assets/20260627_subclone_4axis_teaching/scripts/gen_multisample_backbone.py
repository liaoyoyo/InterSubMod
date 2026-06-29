#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
逐樣本 backbone_resolution.json + candidate_scoring.json 機率回填(Part A / R6)。
對 6 multisample 樣本各 SM_DATA=<dir> 跑 backbone_resolution.py → candidate_scoring.py,
讓 parsimony_first_rank_prob 從 null 回填(純衍生自 topology_per_region.json,無 C++)。

安全機制:
 1. HCC1395(4axis) 重算到暫存目錄,與凍結 backbone/candidate_scoring 對照 MATCH(不覆寫凍結檔);MATCH 失敗即中止,不動 6 樣本。
 2. 6 樣本:回填前後比對 situation_dist 必須不變(只新增 prob,不改既有 R6 數字),變了就警示。
用法: python3 gen_multisample_backbone.py
"""
import os, json, subprocess, shutil, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA4AXIS = os.path.normpath(os.path.join(HERE, "..", "data"))
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
BACKBONE = os.path.join(HERE, "backbone_resolution.py")
SCORING = os.path.join(HERE, "candidate_scoring.py")

def run(script, smdata):
    env = dict(os.environ, SM_DATA=smdata)
    return subprocess.run(["python3", script], env=env, capture_output=True, text=True)

def jload(p):
    return json.load(open(p, encoding="utf-8"))

# ---- 1. HCC1395 MATCH 驗證(暫存,絕不覆寫凍結) ----
print("=== HCC1395(4axis) MATCH 驗證(暫存目錄) ===")
ok = True
with tempfile.TemporaryDirectory() as tmp:
    shutil.copy(os.path.join(DATA4AXIS, "topology_per_region.json"), tmp)
    r = run(BACKBONE, tmp)
    if not os.path.exists(os.path.join(tmp, "backbone_resolution.json")):
        print("backbone 跑失敗:", r.stderr[-500:]); raise SystemExit(1)
    regen_bb = jload(os.path.join(tmp, "backbone_resolution.json"))
    frozen_bb = jload(os.path.join(DATA4AXIS, "backbone_resolution.json"))
    bb_match = regen_bb["summary"] == frozen_bb["summary"]
    print("  backbone summary MATCH:", "✓" if bb_match else "✗")
    if not bb_match:
        print("   regen B1 :", regen_bb["summary"].get("B1_ambiguous_parsimony"))
        print("   frozen B1:", frozen_bb["summary"].get("B1_ambiguous_parsimony"))
        ok = False
    run(SCORING, tmp)
    regen_cs = jload(os.path.join(tmp, "candidate_scoring.json"))
    frozen_cs = jload(os.path.join(DATA4AXIS, "candidate_scoring.json"))
    cs_match = (regen_cs["situation_dist"] == frozen_cs["situation_dist"] and regen_cs["n_total"] == frozen_cs["n_total"])
    regen_nn = sum(1 for q in regen_cs["queue"] if q.get("parsimony_first_rank_prob") is not None)
    frozen_nn = sum(1 for q in frozen_cs["queue"] if q.get("parsimony_first_rank_prob") is not None)
    print("  scoring situation_dist+n_total MATCH:", "✓" if cs_match else "✗", "| prob_nonnull regen/frozen=%d/%d" % (regen_nn, frozen_nn))
    if not cs_match: ok = False
if not ok:
    print("⚠ HCC1395 MATCH 失敗 → 中止,不動 6 樣本(script 與凍結不一致)"); raise SystemExit(1)
print("  → MATCH 通過,泛化邏輯可信。\n")

# ---- 2. 6 樣本回填(situation_dist 不變保護) ----
print("=== 6 樣本 backbone + scoring 回填 ===")
print("%-16s %5s %9s %8s %12s %10s %8s" % ("sample", "B1", "median", "hi>=.7", "prob_nonnull", "順序待定", "sit不變"))
for s in sorted(os.listdir(MSROOT)):
    dr = os.path.join(MSROOT, s)
    tp = os.path.join(dr, "topology_per_region.json")
    if not os.path.exists(tp):
        continue
    csp = os.path.join(dr, "candidate_scoring.json")
    old_sit = jload(csp)["situation_dist"] if os.path.exists(csp) else None
    run(BACKBONE, dr)
    run(SCORING, dr)
    bb = jload(os.path.join(dr, "backbone_resolution.json"))
    cs = jload(csp)
    b1 = bb["summary"]["B1_ambiguous_parsimony"]
    nonnull = sum(1 for q in cs["queue"] if q.get("parsimony_first_rank_prob") is not None)
    order_pending = cs["situation_dist"].get("順序 2-3 順位待定", 0)
    sit_same = "✓" if (old_sit is None or old_sit == cs["situation_dist"]) else "✗變了!"
    print("%-16s %5d %9s %8d %12d %10d %8s" % (
        s, b1["n"], str(b1["median_first_rank_prob"]), b1["high_conf(>=0.7)"], nonnull, order_pending, sit_same))
print("DONE")
