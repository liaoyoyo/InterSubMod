#!/usr/bin/env python3
"""[P1 schema 確認 + cis-test 候選萃取] 讀 records_v6 → 印完整 schema（欄位/sample/分佈）→ 萃取 subclone
候選池（subclone_novel / 8類 A 單標籤多結構 / contingency many1×single_label_somatic）→ land cis_candidates.json
供 P1.2/P1.3 tumor-only per-CpG + normal cis-attribution。純讀，零 binary。"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
recs = json.load(open(f"{A}/phylo_cpp_wg_full_records_v6.json"))
if isinstance(recs, dict):
    # 可能是 {key: record}; 轉成 list
    recs = list(recs.values())
print(f"N records = {len(recs)}")
r0 = recs[0]
print(f"\n=== keys (first record) ===\n{sorted(r0.keys())}")
print(f"\n=== sample record (truncated) ===\n{json.dumps(r0, ensure_ascii=False)[:2000]}")

# 候選相關欄位偵測
print("\n=== candidate-defining fields ===")
for f in ("chrom", "pos", "set", "pc_verdict", "cat8", "class", "category", "contingency_type",
          "ls_hp_p", "ls_hp_disp", "ls_al_p", "best_label_axis", "coarse_ng", "coarse_label",
          "single_label", "som_marker", "sub_marker", "verdict"):
    print(f"  {f:18s} present={f in r0!s:5s} sample={str(r0.get(f))[:60]}")

# 分佈
def dist(field):
    c = Counter(str(r.get(field)) for r in recs)
    print(f"\n=== {field} dist ===\n" + "\n".join(f"  {k}: {v}" for k, v in c.most_common(15)))

for f in ("pc_verdict", "set", "contingency_type"):
    if f in r0:
        dist(f)
# 8類欄位(中文 key) 偵測
for f in ("cat8", "class", "category"):
    if f in r0:
        dist(f)

# ===== 候選萃取 =====
def keyof(r):
    return {"chrom": r.get("chrom"), "pos": r.get("pos"), "set": r.get("set"),
            "pc_verdict": r.get("pc_verdict"), "contingency_type": r.get("contingency_type"),
            "cat8": r.get("cat8") or r.get("class") or r.get("category")}

sub_novel = [keyof(r) for r in recs if r.get("pc_verdict") == "subclone_novel"]
print(f"\n=== 候選池 ===")
print(f"  subclone_novel: {len(sub_novel)}")

# many1 × single_label_somatic（若有欄位）
many1 = [keyof(r) for r in recs if r.get("contingency_type") in ("many1_結構>標籤", "many1")]
print(f"  contingency many1: {len(many1)}")

# 8類 A（中文 key 含 'A_subclone'）
cat_field = "cat8" if "cat8" in r0 else ("class" if "class" in r0 else ("category" if "category" in r0 else None))
clsA = []
if cat_field:
    clsA = [keyof(r) for r in recs if str(r.get(cat_field, "")).startswith("A")]
    print(f"  8類 A ({cat_field}): {len(clsA)}")

# union 候選（去重 by chrom:pos）
seen = {}
for grp, name in ((sub_novel, "subclone_novel"), (many1, "many1"), (clsA, "classA")):
    for r in grp:
        k = f"{r['chrom']}:{r['pos']}"
        seen.setdefault(k, {**r, "src": []})
        seen[k]["src"].append(name)
cands = list(seen.values())
print(f"  union 候選 (去重): {len(cands)}")
print(f"  union TP/FP: TP={sum(1 for r in cands if r['set']=='TP')} FP={sum(1 for r in cands if r['set']=='FP')}")
json.dump(cands, open(f"{A}/cis_candidates.json", "w"), ensure_ascii=False, indent=1)
print(f"\n[-> cis_candidates.json] {len(cands)} loci")
