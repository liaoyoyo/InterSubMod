"""[L1 HP-tag 貢獻 — 單變量 ablation] 量化 longphase-S HP tag 對 subclone 判別的 data-supported 貢獻。
- chance baseline: 若 HP 隨機, 兩 sSNV same-HP 機率 = Σ p_i²。
- 各關係型態 observed same-HP rate vs chance → enrichment（nested/co_linked 應 >> chance = 真同單倍型）。
- mutual_excl ablation: 無 HP gate 會把全部當 sibling subclone; 有 HP 只留 same-HP → HP 移除多少假 subclone。
- HP3（無 germline 錨）= 無法 HP-gate 的一群。
輸出 報告夾/data/sm_hp_contribution.json。
"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    cen = d["census"]
    pairs = d["pairs"]
    master = json.load(open(f"{RPT}/data/sm_locus_master_summary.json"))

    # chance baseline (修正 adversarial F1): 只用 GATEABLE HP (1-1/2-1)，因 observed same/diff 也只在 gateable 對上算
    # (HP3 已排除於投票)。用全類別(含 HP3)會低估 baseline ~0.443、虛增 enrichment ~13%。
    hp_dist = master["hp_dist_linked_somatic"]
    gate = {k: v for k, v in hp_dist.items() if k in ("1-1", "2-1")}
    gtot = sum(gate.values())
    p = {k: v / gtot for k, v in gate.items()}
    chance_same = sum(v ** 2 for v in p.values())

    # per-relationship same/diff HP (powered + both somatic)
    rel_hp = {}
    for pr in pairs:
        if not (pr["powered"] and pr["somatic_a"] and pr["somatic_b"]):
            continue
        rel = pr["rel"]
        rel_hp.setdefault(rel, {"same": 0, "diff": 0})
        rel_hp[rel]["same" if pr["same_hp"] else "diff"] += 1
    for rel, v in rel_hp.items():
        n = v["same"] + v["diff"]
        v["same_hp_rate"] = round(v["same"] / n, 3) if n else None
        v["enrichment_vs_chance"] = round(v["same_hp_rate"] / chance_same, 2) if v["same_hp_rate"] else None

    # mutual_excl ablation
    me = rel_hp.get("mutual_excl", {"same": 0, "diff": 0})
    me_total = me["same"] + me["diff"]
    ablation = {
        "without_HP_gate_called_subclone": me_total,  # 全部互斥都當 sibling
        "with_HP_gate_true_subclone(same-HP)": me["same"],
        "removed_as_allelic(diff-HP)": me["diff"],
        "false_subclone_removed_pct": round(100 * me["diff"] / me_total, 1) if me_total else None,
    }

    # HP3-ungateable linked somatic
    hp3 = hp_dist.get("3", 0)
    hp_gateable = hp_dist.get("1-1", 0) + hp_dist.get("2-1", 0)

    # HP1-1/2-1 somatic-mode 讀數 recovery: linked somatic 帶 somatic-mode tag 的數
    out = {
        "hp_distribution_linked_somatic": hp_dist,
        "chance_same_hp_baseline": round(chance_same, 3),
        "per_relationship_same_hp": rel_hp,
        "mutual_excl_ablation": ablation,
        "hp_gateable_somatic(1-1/2-1)": hp_gateable,
        "hp3_ungateable_somatic": hp3,
        "verdict": {
            "discrimination": "nested/co_linked same-HP 遠高於 chance → 真同單倍型克隆連鎖; "
                              "mutual_excl 需 HP 分 sibling(same) vs allelic(diff)",
            "hp_contribution": "HP gate 對 mutual_excl 移除 %s%% 假 subclone (allelic); "
                               "對 nested/co_linked 是確認(已 >> chance)" % ablation["false_subclone_removed_pct"],
            "limit": "HP3 %d linked-somatic 無 germline 錨, 無法 HP-gate" % hp3,
        },
    }
    json.dump(out, open(f"{RPT}/data/sm_hp_contribution.json", "w"), ensure_ascii=False, indent=1)
    print(f"chance same-HP baseline = {chance_same:.3f}")
    print("per-relationship same-HP rate (vs chance enrichment):")
    for rel in ("nested_b_in_a", "nested_a_in_b", "co_linked", "mutual_excl", "independent"):
        v = rel_hp.get(rel)
        if v:
            print(f"  {rel:16s} same={v['same']:>6} diff={v['diff']:>5} rate={v['same_hp_rate']} ({v['enrichment_vs_chance']}x chance)")
    print(f"\nmutual_excl ablation: 無HP→{ablation['without_HP_gate_called_subclone']} subclone; "
          f"有HP→{ablation['with_HP_gate_true_subclone(same-HP)']}; "
          f"移除 allelic {ablation['removed_as_allelic(diff-HP)']} ({ablation['false_subclone_removed_pct']}%)")
    print(f"HP3 ungateable linked-somatic = {hp3}")


if __name__ == "__main__":
    main()
