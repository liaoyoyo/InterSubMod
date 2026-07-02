#!/usr/bin/env python3
"""[工作站] v3.1 切割+標籤 肉眼確認 spec.json（§13-A 注入）。圖=figs_phylo_v3/ 外部 PNG。"""
import json, subprocess
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
v3render = json.load(open(f"{A}/phylo_v3_render.json"))
v2prev = {x["key"]: x["v2_ngroups"] for x in json.load(open(f"{A}/phylo_v2_final.json"))}
try:
    commit = subprocess.check_output(["git", "-C", A, "rev-parse", "--short", "HEAD"]).decode().strip()
    branch = subprocess.check_output(["git", "-C", A, "rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
except Exception:
    commit, branch = "unknown", "unknown"


def kind_for_class(cl):
    return {"CONFIRMED": "ok", "NEAR_CONFIRMED": "info", "REAL_NOVEL": "warn",
            "REAL_DIFFUSE": "warn", "NO_CLEAR_SPLIT": "no"}.get(cl, "info")


items = []
for o in sorted(v3render, key=lambda z: (-z["v2_ngroups"], -z["n"])):  # v2_ngroups 欄=render 內 v3 群數
    key = o["key"]; ng = o["v2_ngroups"]; prev = v2prev.get(key, "?")
    badges = [{"label": f"v3 {ng} 群", "kind": "ok" if ng >= 2 else "info"},
              {"label": o["fine_class"].replace("_SPLIT", ""), "kind": kind_for_class(o["fine_class"])}]
    if prev != ng:
        badges.append({"label": f"v2 {prev}→v3 {ng}", "kind": "warn"})
    metrics = [{"k": "reads(n)", "v": str(o["n"]), "src": "phylo_v3_render.json"},
               {"k": "CpG", "v": str(o["C"]), "src": "phylo_v3_render.json"},
               {"k": "v3 群數", "v": str(ng), "src": "phylo_v3_render.json"},
               {"k": "v2 群數", "v": str(prev), "src": "phylo_v2_final.json"},
               {"k": "離群", "v": str(o["n_outlier"]), "src": "phylo_v3_render.json"}]
    for L in sorted(o["align"]):
        a = o["align"][L]
        metrics.append({"k": f"群 {L}", "v": f"n{a['n']} · hp={a['hp']} · allele={a['allele']}", "src": "phylo_v3_render.json"})
    alleles = {o["align"][L]["allele"].split("(")[0] for L in o["align"]}
    hps = {o["align"][L]["hp"].split("(")[0] for L in o["align"]}
    if ng < 2:
        guide = "單群：v3 判定無清楚離散結構。看熱圖列是否大致同質 → 確認『不切』是否正確（chr20_42981498 為 borderline 範例）。"
    elif len(alleles) > 1:
        guide = f"多群且各群 allele 不同（{alleles}）= germline REF/ALT 軸（cis-ASM）。看 ALT 側欄/對角塊是否與群對齊。"
    elif len(hps) > 1:
        guide = f"多群且各群 hp 不同（{hps}）= germline 單倍型軸。看 HP 側欄是否與群對齊。"
    else:
        guide = "🔴多群但同 allele+hp = within-germline 多群（subclone 候選 vs cis-ASM 子結構需細看對角塊是否真離散）。"
    items.append({"id": key, "title": key.replace("_", ":"),
                  "badges": badges, "meta_metrics": ["reads(n)", "CpG", "v3 群數"],
                  "metrics": metrics,
                  "figures": {"mode": "png", "imgs": [{"caption": "UPGMA樹(v3群)+甲基+HP/ALT/strand+距離", "path": o["png"]}]},
                  "reading_guide": guide})

spec = {
    "meta": {"title": "phylo-v3 切割+標籤 肉眼確認工作站（FM1 修正）",
             "subtitle": "30 pilot 位點 · v3=v2+quarantine-descend（修 FM1 單離群吃群位；v2→v3 救回 6 真切割, 噪音FP 0%）· 重點=切割與標籤判定是否正確",
             "lskey": "phylo_v3_verify_v1", "tier": "⭐2-3 pilot", "scope": "30 loci",
             "build_commit": commit, "build_branch": branch},
    "changelog": {
        "data_status": {"summary": "v3.1 標籤從 phylo_v3_render.json 注入（快取矩陣重算 L1）；圖 figs_phylo_v3/ 外部 PNG。"},
        "binary_status": {"summary": "純 Python 分析；C++ binary 未改（feat/summary-nreadsvalid）。"},
        "phase": "pilot 方法確認（FM1 修正後切割+標籤正確性）",
        "corrections": [
            {"id": "C1", "what": "double-dip 修正(per-subgroup 重分群 null)", "status": "in-HEAD", "effect": "v1 假 subclone 候選消失", "src": "phylo_compare_v1v2.py"},
            {"id": "C2", "what": "FM1 復發修正(quarantine-descend 物有所值)", "status": "in-HEAD", "effect": "v2→v3 救回 6 真切割(chr22_26939195→2/30454004→4/21_40426852→2), 噪音FP 0%", "src": "phylo_v3_validate.py"}
        ],
        "audit_conclusion": {"summary": "v3.1 修 FM1 後; 用戶肉眼旗標 chr22_26939195/30454004/21_40426852 全救回; chr20_42981498 真 borderline 仍 1 群。"}
    },
    "sections": [{"id": "intro", "title": "① 怎麼看（重點=切割+標籤正確性）",
                  "html": "<p class='sub'>每位點是一筆<b>資料</b>。確認兩件事：<b>（1）切割對不對</b>—v3 切的群（樹分支色/距離對角塊）是否符合熱圖結構？<b>（2）標籤對不對</b>—各群 hp/allele 側欄是否與群一致、階層標籤合理？<br><b>v3 已修 FM1</b>（單離群不再吃整群；v2→v3 救回 6 個你眼睛看到的切割）。每項記『切割+標籤都對 / 切割存疑 / 標籤存疑 / 錯誤』。</p>"}],
    "item_config": {
        "verdict_states": [{"key": "ok", "label": "切割+標籤都對"}, {"key": "cut", "label": "切割存疑"},
                           {"key": "lab", "label": "標籤存疑"}, {"key": "wrong", "label": "錯誤"}],
        "reason_options": ["切多了(亂切)", "切少了(漏結構)", "標籤對齊錯", "離群判定錯", "邊界難判", "其他"],
        "required_metrics": ["reads(n)", "CpG", "v3 群數"], "per_page": 30
    },
    "items": items,
    "provenance": {"sources": ["phylo_v3_render.json", "phylo_v2_final.json", "phylo_v3_validate.json"],
                   "note": "v3.1 切割+標籤；圖 figs_phylo_v3/"}
}
json.dump(spec, open(f"{A}/phylo_v3_workstation_spec.json", "w"), ensure_ascii=False, indent=1)
print(f"WROTE spec: {len(items)} items · commit {commit} · branch {branch}")
