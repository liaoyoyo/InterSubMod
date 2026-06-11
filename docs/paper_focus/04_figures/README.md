<!--
建立時間: 2026-06-09
狀態: library (圖片 — 論文需要的圖 manifest)
報告類型: paper_focus_figure_manifest
受眾: 廖子游 · PI
provenance_note: 圖對應的數字沿用已驗證集合（🟢P/🟡S）；🟡 圖須等原檔對賬才可放真值。
-->
<!-- provenance-verified: 圖清單由 02_paper_framework + 01_focus_notes/03 方向卡「預期圖表」導出。 -->

# 圖片 manifest — 論文需要的圖

> **L0 一眼結論**：論文需 ~10 張圖；**5 張今天就有 🟢 料可畫（NEGATIVE + phasing）**、4 張要等 🟡 原檔對賬、1 張（catalog）要先做 ISM-4。圖檔（PNG/SVG）放本資料夾。
>
> **L1 重點邏輯**：① 🟢 圖優先（最硬、不卡）；② 🟡 圖等 `02 文件庫` 對賬後才放真值（否則捏造）；③ 方法示意/三軸圖用 `/methods-example` 或 `/image-gen`，真實數據圖用 matplotlib + 人眼驗證。

---

## L2 — 圖清單

| 圖 | 用在哪節 | 內容 | 料(tier) | 狀態 | 工具 |
|----|---------|------|---------|------|------|
| F1 三軸結構圖 | Intro/總覽 | 主 phasing / 支撐 ASM / 底座 NEGATIVE | 🔵 framing | 可做 | /methods-example |
| F2 ISM 方法示意 | Methods | read×CpG→距離→分群→cis-test | 🔵 schematic | 可做 | /methods-example |
| F3 filter 死四道對照 | Results R1 | in-dist+0.0224 vs held-out−0.0001 / AUC=0.505 vs null / ~4%sens / ARI 0.135 vs null | 🟢 | **可做** | matplotlib |
| F4 NG=2 Inner>Outer forest | Results R2 | 7 樣本 gap + bootstrap CI[0.104,0.459] | 🟢 | **可做** | matplotlib |
| F5 ASM TP/FP/FN rate | Results R3 | subhap-matched 並列 ~3.5×；COLO829 判別消失 | 🟢 | **可做** | matplotlib |
| F6 copy-dosage falsification | Results R3 | CN-class \|Δβ\| 平坦 + signed ρ 反向 | 🟢 | **可做** | matplotlib |
| F7 AF→NGroups=phasing | Results R4 | r=0.656 + 控甲基 collider + 分布 99.5%≤2 | 🟢 | **可做** | matplotlib |
| F8 位點甲基分群 catalog | Results R6 | 每位點標籤 heatmap/表（reliable/Δβ/somatic/cis）| ⭐ 待 ISM-4 | 待做 | matplotlib |
| F9 chr17 cis-test 三路 | Results R3 | normal HP1/tumor HP1/tumor HP1-1 read×CpG | 🟡 | 待對賬 | matplotlib |
| F10 cross-sample 6/6 excess | Results R5 | +0.101–0.241，3 癌種，標 COLO829 最薄 | 🟡 | 待對賬(*_gwasm.json) | matplotlib |
| F11 V10 甲基非-copy | Discussion | depth-matched normal 0.979 > tumor 0.866 | 🟡 | 待對賬 | matplotlib |

---

## L1 — 提醒
- 🟢 圖（F3-F7）可立即畫；🟡 圖（F9-F11）**等原檔對賬才放真值**（否則撞 §13 捏造）。
- 真實數據圖用 `InterSubMod/scripts/lib/plot_setup.py`（CJK 字型）；示意圖用 /methods-example 或 /image-gen + /image-vision-check。
- 生成的 PNG 放本資料夾，命名 `F{NN}_{slug}.png`。
