# verify-workstation 範例

## selftest_spec.json（最小可跑範例）

```bash
python3 ../tools/build_workstation.py selftest_spec.json -o /tmp/selftest.html
```

2 項：一張 inline SVG 圖 + 一張外部 PNG 圖、changelog、自訂 section、required_metric=CV。
移除某項的 CV 再跑 → 預期 **exit 3 refuse**（§13-A）。

## 兩個真實 worked 實例（兩種尺度）

| 尺度 | 實例 | generator | 輸出 | 特徵 |
|---|---|---|---|---|
| 小（~10 項） | chr2:18M subclone | `InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/build_workstation_html.py` | `docs/explain/07_subclone-judgment-workstation-chr2-18M.standalone.html` | inline SVG、computed methyl class、conflict 校正表（authoritative-vs-draft changelog） |
| genome（30,350 項） | HCC1395 ASM 篩選 | `InterSubMod/research/tsg_promoter_asm_reviewer/scripts/102_build_workstation.py` | `…/display_v2/20260610_judgment_workstation_01.html` | compact 22-col array、60,700 外部 PNG、Tier A/B + 5 層驗證 + 即時門檻 + redesign loop + §⓪b corrected-ISM changelog |

- 小尺度 → 用本 skill 的 `tools/build_workstation.py`（`figures.mode="svg"`，單檔可攜）。
- genome-scale → 沿用 102 式 compact array + 外部 PNG（見 `../references/section_modules.md`）；§⓪b 修正過程紀錄即本 skill changelog 模組。

兩者共用：localStorage 判讀（同意/存疑/否定 + 原因）、JSON/CSV 匯出/匯入、⌖ provenance、§13-A 注入/refuse。
