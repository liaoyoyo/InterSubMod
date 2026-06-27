#!/usr/bin/env python3
"""[工作站] 從 phylo_v2_render.json 注入 → verify-workstation spec.json（§13-A：數字不手打）。
30 位點 v2 切割+標籤，供肉眼逐一確認『切得對不對、標籤判得對不對』。圖=figs_phylo_v2/ 外部 PNG。"""
import json, os, subprocess
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
v2 = json.load(open(f"{A}/phylo_v2_render.json"))
v1 = {x["key"]: x["n_groups"] for x in json.load(open(f"{A}/phylo_groups.json"))}
try:
    commit = subprocess.check_output(["git", "-C", A, "rev-parse", "--short", "HEAD"]).decode().strip()
    branch = subprocess.check_output(["git", "-C", A, "rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
except Exception:
    commit, branch = "unknown", "unknown"


def kind_for_class(cl):
    return {"CONFIRMED": "ok", "NEAR_CONFIRMED": "info", "REAL_NOVEL": "warn",
            "REAL_DIFFUSE": "warn", "NO_CLEAR_SPLIT": "no"}.get(cl, "info")


items = []
for o in sorted(v2, key=lambda z: (-z["v2_ngroups"], -z["n"])):
    key = o["key"]; ng = o["v2_ngroups"]; v1g = v1.get(key, "?")
    badges = [{"label": f"v2 {ng} 群", "kind": "ok" if ng >= 2 else "info"},
              {"label": o["fine_class"].replace("_SPLIT", ""), "kind": kind_for_class(o["fine_class"])}]
    if v1g != ng:
        badges.append({"label": f"v1 {v1g}→v2 {ng}", "kind": "warn"})
    # metrics: n / CpG / v2群 / 各群對齊
    metrics = [{"k": "reads(n)", "v": str(o["n"]), "src": "phylo_v2_render.json"},
               {"k": "CpG", "v": str(o["C"]), "src": "phylo_v2_render.json"},
               {"k": "v2 群數", "v": str(ng), "src": "phylo_v2_render.json"},
               {"k": "v1 群數", "v": str(v1g), "src": "phylo_groups.json"},
               {"k": "離群", "v": str(o["n_outlier"]), "src": "phylo_v2_render.json"}]
    for L in sorted(o["align"]):
        a = o["align"][L]
        metrics.append({"k": f"群 {L}", "v": f"n{a['n']} · hp={a['hp']} · allele={a['allele']}", "src": "phylo_v2_render.json"})
    # 對齊型態判定（給 reading guide 用）
    alleles = {o["align"][L]["allele"].split("(")[0] for L in o["align"]}
    hps = {o["align"][L]["hp"].split("(")[0] for L in o["align"]}
    if ng < 2:
        guide = "單群：v2 判定無清楚離散結構。肉眼看熱圖列是否大致同質（無明顯分塊）→ 確認『不切』是否正確。"
    elif len(alleles) > 1:
        guide = f"多群且各群 allele 不同（{alleles}）= germline REF/ALT 軸分裂（cis-ASM）。看 ALT 側欄與群是否對齊、對角塊是否清楚。"
    elif len(hps) > 1:
        guide = f"多群且各群 hp 不同（{hps}）= germline 單倍型軸分裂。看 HP 側欄與群是否對齊。"
    else:
        guide = "多群但同 allele+hp = within-germline 多群（subclone 候選需留意）。看對角塊是否真的分離。"
    items.append({"id": key, "title": key.replace("_", ":"),
                  "badges": badges, "meta_metrics": ["reads(n)", "CpG", "v2 群數"],
                  "metrics": metrics,
                  "figures": {"mode": "png", "imgs": [{"caption": "UPGMA樹(v2群)+甲基+HP/ALT/strand+距離", "path": o["png"]}]},
                  "reading_guide": guide})

spec = {
    "meta": {"title": "phylo-v2 切割+標籤 肉眼確認工作站",
             "subtitle": f"30 pilot 位點（HCC1395 tumor-only, chr20-22）· v2=per-subgroup 重分群 null+RNULL40（已修 double-dip）· 重點=切割與標籤判定是否正確（非 TP/FP 判別）",
             "lskey": "phylo_v2_verify_v1", "tier": "⭐2-3 pilot", "scope": "30 loci",
             "build_commit": commit, "build_branch": branch},
    "changelog": {
        "data_status": {"summary": "v2 標籤從 phylo_v2_render.json 注入（快取矩陣重算，L1）；圖 figs_phylo_v2/ 外部 PNG。"},
        "binary_status": {"summary": "純 Python 分析；C++ binary 未改（feat/summary-nreadsvalid）。"},
        "phase": "pilot 方法確認（切割+標籤正確性）",
        "corrections": [
            {"id": "C1", "what": "phylo split null double-dip 修正（沿用樹分割→per-subgroup 重分群）", "status": "in-HEAD", "effect": "v1 假 subclone 候選消失（chr20_30274614/chr21_10353822）", "src": "phylo_compare_v1v2.py"},
            {"id": "C2", "what": "遞迴 null 全域→子群內；RNULL 12→40", "status": "in-HEAD", "effect": "深層標籤校正", "src": "phylo_v2_final.py"}
        ],
        "audit_conclusion": {"summary": "本工作站供肉眼確認 v2 每位點切割與標籤；判讀 localStorage 留存 + 匯出 JSON/CSV。"}
    },
    "sections": [{"id": "intro", "title": "① 怎麼看（重點=切割+標籤正確性，非 TP/FP）",
                  "html": "<p class='sub'>每個位點是一筆<b>資料</b>（不分 TP/FP）。請肉眼確認兩件事：<b>（1）切割對不對</b>—v2 切出的群（樹分支顏色 / 距離對角塊）是否符合熱圖看到的結構？該切的有切、不該切的有沒有亂切？<b>（2）標籤判定對不對</b>—各群的 hp/allele 對齊（側欄）是否與群一致？階層標籤（1/1-1/2-1）是否合理？<br>每項記『切割+標籤都對 / 切割存疑 / 標籤存疑 / 錯誤』。</p>"}],
    "item_config": {
        "verdict_states": [{"key": "ok", "label": "切割+標籤都對"}, {"key": "cut", "label": "切割存疑"},
                           {"key": "lab", "label": "標籤存疑"}, {"key": "wrong", "label": "錯誤"}],
        "reason_options": ["切多了(亂切)", "切少了(漏結構)", "標籤對齊錯", "離群判定錯", "邊界難判", "其他"],
        "required_metrics": ["reads(n)", "CpG", "v2 群數"], "per_page": 30
    },
    "items": items,
    "provenance": {"sources": ["phylo_v2_render.json", "phylo_groups.json", "phylo_compare_v1v2.json"],
                   "note": "v2 切割+標籤；圖 figs_phylo_v2/"}
}
json.dump(spec, open(f"{A}/phylo_v2_workstation_spec.json", "w"), ensure_ascii=False, indent=1)
print(f"WROTE spec: {len(items)} items · commit {commit} · branch {branch}")
