#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Part C — build the /verify-workstation spec from cases.json and render the HTML.
4-way human labels: 真結構 / germline-only / 離散度假象 / sample-confound."""
import json, os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
REPO = "/big7_disk/liaoyoyo2001/InterSubMod"
GEN = f"{REPO}/.claude/skills/verify-workstation/tools/build_workstation.py"
cases = json.load(open(os.path.join(OUTDIR, "cases.json")))

def head():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO).decode().strip()
    except Exception:
        return "?"

STRATUM = {
    "high_label_F": ("高 label-F", "warn"),
    "near_boundary": ("近邊界", "info"),
    "sample_only": ("sample-only", "gold"),
    "germline_asm": ("germline-ASM", "ok"),
    "candidate_som": ("candidate-somatic", "no"),
}

def fnum(v, nd=2):
    try:
        return f"{float(v):.{nd}f}"
    except (ValueError, TypeError):
        return str(v) if v not in ("", None) else "NA"

items = []
for c in cases["cases"]:
    m = c["metrics"]
    slab, skind = STRATUM.get(c["stratum"], (c["stratum"], "mut"))
    badges = [{"label": slab, "kind": skind}, {"label": c["cls"], "kind": "no" if c["cls"] == "FP" else "mut"},
              {"label": f"VC:{m.get('VerificationClass','?')}", "kind": "info"},
              {"label": f"強度{m.get('StrengthGrade','?')}", "kind": "mut"}]
    metrics = [
        {"k": "HP-AUC(all/tum/norm)", "v": f"{fnum(m.get('HP_AUC_All'))}/{fnum(m.get('HP_AUC_Tumor'))}/{fnum(m.get('HP_AUC_Normal'))}", "src": "phase2 significance_summary"},
        {"k": "label-F / p", "v": f"{fnum(m.get('LabelHPPermanovaF'),1)} / {fnum(m.get('LabelHPPermanovaP'),3)}", "src": "phase2"},
        {"k": "cluster-F / p", "v": f"{fnum(m.get('ClusterPermanovaF'),1)} / {fnum(m.get('ClusterPermanovaP'),3)}", "src": "phase2"},
        {"k": "CramersV", "v": fnum(m.get("CramersV")), "src": "phase2"},
        {"k": "VC (legacy→new)", "v": f"{m.get('VerificationClass_Legacy','?')} → {m.get('VerificationClass','?')}", "src": "phase2"},
        {"k": "Strength", "v": f"{fnum(m.get('StrengthScore'))} ({m.get('StrengthGrade','?')})", "src": "phase2"},
        {"k": "NumReads/Valid", "v": f"{m.get('NumReads','?')}/{m.get('NReadsValid','?')}", "src": "phase2"},
    ]
    figs = {"mode": "png", "imgs": [{"caption": "左 read×CpG 甲基 / 右 read×read NHD 距離（HP 分群）", "path": f"figs/{c['fig']}"}]} if c.get("fig") else None
    items.append({
        "id": c["id"], "title": f"{c['chrom']}:{c['pos']}",
        "subtitle": f"被 verdict 忽略的顯著 PERMANOVA（legacy {m.get('VerificationClass_Legacy','?')}）",
        "badges": badges, "metrics": metrics, "modal_stats": metrics,
        "meta_metrics": ["HP-AUC(all/tum/norm)", "label-F / p"],
        "figures": figs,
        "reading_guide": "判讀：右圖距離矩陣對角暗塊『對齊左側 HP 側欄』⇒ 真結構；甲基成帶但不對齊 HP（按甲基高低分）⇒ 可能 sample-confound；無塊/散亂高熵 ⇒ 離散度假象；若分群其實是 germline 單倍型既有差異 ⇒ germline-only。",
    })

spec = {
    "meta": {
        "title": "Part C — 被 verdict 忽略的顯著 PERMANOVA · 分層抽樣肉眼判讀工作站",
        "subtitle": f"Phase-2 corrected（±5000/SKIP/develop binary）· {cases['n_cases']} 抽樣 / 母體 {cases['n_target']} · 逐項標 真結構/germline-only/離散度假象/sample-confound",
        "lskey": "ism_phase2_inspection_v1", "tier": "L2-pilot", "scope": "PARTIAL — 分層抽樣（非全 7271）",
        "build_commit": head(), "build_branch": "research/subclonal-reconstruction-202606",
    },
    "banner": {"text": "⚠ <b>分層抽樣肉眼工作站</b>：母體 = Phase-2 corrected 中『被舊 verdict 忽略但 HP-aligned PERMANOVA 顯著』的 7,271 位點（即 reclassify 救回的假陰性候選）。此處每層抽 8 個供你逐一判真偽。"},
    "changelog": {
        "data_status": {"summary": "Phase-2 corrected ISM（develop binary、±5000、SKIP）。VerificationClass=6 類簡版、含 StrengthScore。"},
        "phase": {"current": "Part C 分層抽樣肉眼判讀", "next": "germline vs candidate-somatic 分層需 normal-BAM 重跑（develop Phase A-D）— 本 tumor-only 掃描 Δβ-sig 全 0"},
        "corrections": [
            {"id": "TARGET", "what": "目標集 = legacy VC∈{Weak,Noise} 但 LabelHP-PERMANOVA 顯著（被 verdict 忽略的結構）", "status": "authoritative", "effect": "母體 7,271", "src": "select_and_render.py"},
            {"id": "STRATA", "what": "3 分層成功：高 label-F / 近邊界(HP-AUC≈0.7) / sample-only(cluster 顯著但 HP-AUC≤0.6)", "status": "authoritative", "effect": "各抽 8", "src": "select_and_render.py"},
            {"id": "GERMLINE-STRATUM", "what": "germline vs candidate-somatic 分層：tumor-only 掃描下 Δβ-sig 全 0 → 需 normal-BAM 重跑", "status": "not-done", "effect": "本工作站暫缺此層；4-way 標記仍可標 germline-only", "src": "genome-wide GermlineAsmDbeta_Sig=0 (phase2)"},
        ],
        "audit_conclusion": {"summary": "PERMANOVA 過敏感（valid→100%顯著），故『顯著 PERMANOVA』非強判別；真偽靠肉眼看距離塊是否對齊 HP + 甲基結構。這正是『自動判不出才用 HTML』。"},
    },
    "sections": [{
        "id": "howto", "title": "① 怎麼用（4-way 判讀）",
        "html": "<div class='keybox' style='background:#07261c;border-left:4px solid #22c55e;padding:11px 14px;border-radius:8px'>"
                "<b style='color:#86efac'>母體</b>＝被舊 verdict 忽略、但 HP-aligned PERMANOVA 顯著的位點（reclassify 救回的假陰性候選）。逐項看圖後標 4 類之一：<br>"
                "<b style='color:#86efac'>真結構</b>＝距離矩陣對角暗塊對齊 HP 側欄；<b style='color:#fca5a5'>germline-only</b>＝分群其實是既有 germline 單倍型差異（非 tumor 新生）；"
                "<b style='color:#fcd34d'>離散度假象</b>＝無乾淨塊、高熵散亂（betadisper 型）；<b style='color:#7dd3fc'>sample-confound</b>＝按甲基高低/tumor-normal 分而非 HP。"
                "判讀存 localStorage，可匯出 JSON/CSV。</div>",
    }],
    "item_config": {
        "verdict_states": [
            {"key": "real", "label": "真結構"}, {"key": "germline", "label": "germline-only"},
            {"key": "dispersion", "label": "離散度假象"}, {"key": "confound", "label": "sample-confound"},
        ],
        "reason_options": [{"key": "", "label": "— 信心 —"}, {"key": "clear", "label": "明確"},
                           {"key": "ambiguous", "label": "模糊"}, {"key": "other", "label": "其他"}],
        "required_metrics": ["HP-AUC(all/tum/norm)", "label-F / p"], "per_page": 12,
    },
    "items": items,
    "provenance": {
        "sources": ["/big7_disk/liaoyoyo2001/ism_phase2_scan/HCC1395_{tp,fp}/significance_summary.csv",
                    "docs/experiments/in_progress/2026/06/20260616_phase2_inspection_workstation/{cases.json,figs/}",
                    "binary: develop snapshot /big7_disk/liaoyoyo2001/ism_phase2/bin/inter_sub_mod_develop (22e5c97)"],
        "note": "Part C 分層抽樣肉眼判讀；圖由 select_and_render.py 從 phase2 per-region 輸出渲染（read×CpG 甲基 + read×read NHD 距離）。",
    },
}
spec_path = os.path.join(OUTDIR, "spec.json")
json.dump(spec, open(spec_path, "w"), ensure_ascii=False, indent=1)
out_html = os.path.join(OUTDIR, "20260616_inspection_workstation.html")
r = subprocess.run(["python3", GEN, spec_path, "-o", out_html], capture_output=True, text=True)
print(r.stdout, r.stderr)
print("HTML:", out_html)
