# 位點判斷與顯示確認 (display_v2) — HCC1395 甲基差異篩選

互動顯示頁，讓使用者逐位點檢視「哪些位點被判為有意義甲基差異、哪些被篩掉、為什麼」，並調整篩選方法。

## 開啟

```
20260609_locus_judgment_display_01.html
```

⚠ **需與同目錄 `figs/`（4140 張 PNG）一起開**（圖用相對路徑外部載入；4140 張無法 base64 內嵌成單檔）。
`figs/` 已 gitignore（太大）→ clone 後須重生（見下）。HTML 本身已內嵌所有數值，文字/表格/門檻試算可獨立運作，只有縮圖需 `figs/`。

## 內容（四區段）

1. **篩選漏斗**：caller→longphase-S→ISM→gate 真實數字
2. **即時門檻試算**（全 30,350 位點）：拉 CramersV/|Δβ|/reads 看 TP/FP 取捨
3. **位點圖庫**：2,070 策展位點 × 2 圖（甲基 read×CpG + ISM read×read 距離矩陣）
   - PASS（通過篩選 1,308）/ FILTERED-OUT（被篩掉 762）兩分頁
   - ⚡ **latent**（487）= CramersV 卡方因表稀疏判不可靠歸零、但 PERMANOVA 仍顯著的真結構
4. **精修 gate**：PERMANOVA 納回 latent 的 minN 掃描表 + 誠實 TP/FP 判讀

## 重生（clone 後 / 數據更新）

```bash
cd InterSubMod
# 1. 策展位點分類（讀 ism_existence_scan 凍結 summary）
python3 research/tsg_promoter_asm_reviewer/scripts/83_prep_curated_loci.py
# 2. 對策展位點重跑 ISM（開 distance matrix）→ /big7_disk/.../ism_display_scan  (~286M, ~70s)
bash research/tsg_promoter_asm_reviewer/scripts/84_run_ism_curated.sh
# 3. 渲染 4140 雙圖 + manifest（多進程, ~3min）
python3 research/tsg_promoter_asm_reviewer/scripts/85_render_dual_figs.py
# 4. 修正 manifest（raw/gated CramersV 分離 + latent 旗標）
python3 research/tsg_promoter_asm_reviewer/scripts/87_fix_manifest.py
# 5. 精修 gate PERMANOVA 掃描
python3 research/tsg_promoter_asm_reviewer/scripts/88_refined_gate_permanova.py
# 6. 建 HTML
python3 research/tsg_promoter_asm_reviewer/scripts/86_build_display_html.py
```

## 數據來源（凍結 baseline，勿改）

- `ism_existence_scan/HCC1395_{tp,fp}/significance_summary.csv` — 117 欄/位點，Significant gate 真值
- gate = `PassedGating & global_p≤0.05 & CramersV≥0.1 & NumReads≥20`
- CramersV = `cramers_v_reliable ? cramers_v : 0`（可靠性 = 最小期望格≥5, Cochran；`MathUtils.cpp` / `RegionProcessor.cpp:1592`）

## 已知結論邊界

- 單樣本 HCC1395、單一 pipeline → ⭐3 上限。
- latent 納回提升 **characterization 完整性**，非 TP/FP 分辨力（甲基→filter 已 concluded NEGATIVE）。
- 不可把任何子集寫成 variant filter。
