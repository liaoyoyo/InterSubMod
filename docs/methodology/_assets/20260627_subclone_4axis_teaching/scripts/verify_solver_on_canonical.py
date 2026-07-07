#!/usr/bin/env python3
"""Step3-b: 新 solver 對真實 HCC1395 canonical 區的驗證 + 與現有 build_hidden_node_tree/solve_topology 單樹輸出比對。
每類挑代表區,跑 enumerate_min_trees + V1-V7,並對照舊 detail 的 edges/determinacy。"""
import json, os
from collections import Counter
import tree_enumeration_solver as S

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
d = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))
detail = {x["region"]: x for x in d["detail"]}

def pick(pred, n=3):
    return [x for x in d["detail"] if pred(x)][:n]

buckets = {
    "A_determined":       pick(lambda x: x["determinacy"].startswith("A_determined")),
    "A_ambiguous_order":  pick(lambda x: x["determinacy"].startswith("A_ambiguous")),
    "recurrence":         pick(lambda x: x["determinacy"].startswith("recurrence")),
    "E_subcube(pairwise)":pick(lambda x: "pairwise樹" in x["determinacy"]),
    "E_subcube(含衝突)":  pick(lambda x: "含衝突" in x["determinacy"]),
    "incompatible":       pick(lambda x: x["determinacy"] == "incompatible"),
    "B_pairwise":         pick(lambda x: x["determinacy"].startswith("B_pairwise")),
}

summary = Counter()
allpass = True
for cls, regs in buckets.items():
    print("\n" + "=" * 96)
    print(f"### {cls}  ({len(regs)} 代表區)")
    print("=" * 96)
    for r in regs:
        pops = r.get("populations") or {}
        subs = r.get("subread_groups") or []
        # subread_groups 可能是 dict{str:count} 或 list[str]
        if isinstance(subs, dict):
            part = list(subs.keys())
        elif isinstance(subs, list):
            part = list(subs)
        else:
            part = []
        k = r.get("genotype_len") or (len(next(iter(pops))) if pops else (len(part[0]) if part else 0))
        if k == 0:
            print(f"  {r['region']}: skip(k=0)")
            continue
        res = S.enumerate_min_trees(pops, part, k)
        newcls = S.classify(res)
        ver = S.verify_all(res, pops, part, k)
        ok = ver["overall"]
        allpass = allpass and ok
        summary[f"{cls}->{newcls}"] += 1
        vflags = " ".join(f"{a}={'✓' if b[0] else ('–' if b[0] is None else '✗')}" for a, b in ver.items() if a != "overall")
        print(f"\n  region={r['region']}  n_sSNV={r['n_sSNV']} cn={r.get('cn')}")
        print(f"    OLD determinacy={r['determinacy']}  OLD n_edges={len(r.get('edges',[]))}  OLD topology={r['topology_type']}")
        print(f"    NEW class={newcls}  n_trees={res['n_trees']}  n_hidden={res['n_hidden']}  capped={res['capped']}"
              + (f" cap_reason={res['cap_reason']}" if res['capped'] else ""))
        print(f"    verify: {vflags}  overall={'✓PASS' if ok else '✗FAIL'}")
        for i, t in enumerate(res["trees"][:3]):
            print(f"      tree{i}: {t['edges']}  rec={t['recurrence']}")
        if res["n_trees"] > 3:
            print(f"      ... (+{res['n_trees']-3} more minimal trees)")
        for a, b in ver.items():
            if a != "overall" and b[0] is False:
                print(f"    ⚠ {a} FAIL: {b[1]}")

print("\n" + "=" * 96)
print("OLD→NEW 分類遷移總結:")
for k2, v in sorted(summary.items()):
    print(f"  {v:3d}  {k2}")
print(f"\n全區 V1-V7: {'ALL PASS ✓' if allpass else 'SOME FAIL ✗'}")
print("=" * 96)
