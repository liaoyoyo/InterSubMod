#!/usr/bin/env bash
# palette_drift_advisor.sh — PostToolUse(Write|Edit) advisory.
#
# WHY (2026-06-23): ISM renderer 視覺規約(色盤/側欄)有 SoT (scripts/ism_heatmap_std.py + memory
# feedback_ism_case_heatmap_standard_sidebars, 自 2026-06-18) 卻反覆漂移 — 同一個「HP1=#dc2626 撞
# ALT 紅」bug 被 copy-paste 繁殖到多個 renderer (geom/phylo_cpp_render/6 個 phylo_*_render)。純紀錄
# 不足(§13 教訓: 需機械防線)。本 hook NUDGE(advisory, 永不阻擋) 當 .py renderer 本地硬編 HP/ALT/群
# 色盤而非 import ism_heatmap_std。
#
# Scope: ONLY .py 含本地語意色 dict literal (HP 子標籤 "1-1"/"2-1" → hex / ALT-REF dict / 已知撞色 hex).
#        排除 SoT 模組與稽核腳本本身。Exit ALWAYS 0 (advisory, fail-OPEN). 對應稽核 scripts/palette_drift_check.py。
# SoT: docs/references/20260623_ism_observation_panel_and_visual_convention_spec_01.md

INPUT=$(cat 2>/dev/null)
FP=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)
[ -z "$FP" ] && exit 0
case "$FP" in *.py) ;; *) exit 0 ;; esac
[ -f "$FP" ] || exit 0
# 排除 canonical 模組 + 稽核腳本(它們合法含 canonical hex)
case "$(basename "$FP")" in ism_heatmap_std.py|palette_drift_check.py) exit 0 ;; esac

# signature: HP 子標籤 → hex literal (任何本地 HP dict, 不論值對錯=漂移風險), 或 named 色盤含已知撞色 hex
if grep -qE "['\"](1-1|2-1)['\"][[:space:]]*:[[:space:]]*['\"]#[0-9a-fA-F]{6}" "$FP" 2>/dev/null \
   || grep -qE "(HPC|alc|ALC|HP_COL|ALT_COL|PALETTE)[[:space:]]*=[[:space:]]*[{[][^=]*['\"]#(dc2626|2563eb|f87171|f59e0b|7c3aed|0ea5e9)" "$FP" 2>/dev/null; then
  {
    echo "[palette-drift-advisor] ⚠️ $(basename "$FP") 似乎本地硬編 HP/ALT/群 色盤(語意撞色, 如 HP1=#dc2626 撞 ALT 紅)。"
    echo "    規約: ISM renderer 一律 import ism_heatmap_std; 用 H.HP_COL/H.ALT_COL/H.CLUSTER_COL + mpl_dual_panel/sidebar_specs, 禁本地色 dict。"
    echo "    稽核: python3 scripts/palette_drift_check.py  |  SoT: docs/references/20260623_ism_observation_panel_and_visual_convention_spec_01.md"
  } >&2
fi
exit 0
