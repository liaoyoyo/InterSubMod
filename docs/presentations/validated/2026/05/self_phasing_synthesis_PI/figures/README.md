# self_phasing_synthesis_PI/figures/ — 復原狀態（2026-05-11 完整版）

此目錄因 2026-05-11 git history rewrite 操作意外影響，原 13 PNG 全部消失。
復原工作於同日完成 7/13（54%），剩 6 個 F* master 圖待手動重產。

## ✅ 已復原（7/13，54%）

### IGV 截圖（3 個）— 從 v3/figures/ cp 而來

| 檔案 | 大小 | 來源 |
|------|------|------|
| `igv/D_SP1_chr19_17565944.png` | 169 KB | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/figures/` |
| `igv/D_SP2_chr19_12452332.png` | 164 KB | 同上 |
| `igv/D_SP3_chr19_12467180.png` | 151 KB | 同上 |

### Generation 概念圖（4 個）— 由 codex CLI `$imagegen` 重新生成

| 檔案 | 大小 | Prompt 來源 |
|------|------|------------|
| `G1_player_as_referee.png` | 68 KB | `../prompts/G1_player_as_referee.txt` (1694B) |
| `G2_pass_two_flow.png` | 52 KB | `../prompts/G2_pass_two_flow.txt` (1744B) |
| `G3_getVote_three_layer.png` | 130 KB | `../prompts/G3_getVote_three_layer.txt` (2097B) |
| `G4_germline_absent_three_versions.png` | 53 KB | `../prompts/G4_germline_absent_three_versions.txt` (2499B) |

**生成方式**：codex exec --full-auto + `$imagegen` (gpt-image-2 via ChatGPT OAuth)
**Token 使用**：~60k tokens × 4 = ~240k tokens (ChatGPT subscription quota)
**注意**：重生圖片內容**接近但非 byte-equal** 原圖；需人眼校對是否符合 PPT/HTML 用途

## ⚠️ 永久遺失（6/13，46%）— 待手動重產

### Master 數據圖（6 個）— 在 `master/` 子目錄

| 檔案 | 描述 | 重產建議 |
|------|------|---------|
| `master/F1_priority_bug_mechanism.png` | Priority bug 機制示意 | 找 paired_priority_bug_audit/ scripts |
| `master/F2_priority_bug_per_chr_enrichment.png` | 每染色體 enrichment | 從 phaseB/C/D run 數據重產 |
| `master/F3_binary_commit_timeline.png` | C++ binary commit 時間軸 | 手繪或 mermaid + 截圖 |
| `master/F4_chr19_752_victims_scatter.png` | chr19:752 受害位點 scatter | matplotlib + .tsv 重產 |
| `master/F5_layer15_zero_sum_4quadrant.png` | Layer 1.5 zero-sum 4 象限 | matplotlib 重產 |
| `master/F6_paired_vs_TO_HP_distribution.png` | Paired vs TO HP 分佈 | matplotlib 重產 |

### 重產候選來源

1. **`InterSubMod/research/paired_priority_bug_audit/scripts/`** — 含 phaseB/C/D 分析 .py
2. **PI Drive / Google Drive** — PI 可能有副本
3. **重新跑 .py + .tsv** — 數據還在 .gitignore 的 phaseB/C/D run 目錄

## 📦 引用方

`build_html.py` + `build_pptx.py` 引用這些 PNG 於 5/11 HTML preview slides。

開啟以下 HTML 即可看到目前可用的 7 PNG：
- `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/index.html`
- `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/slide_*.html`

## 業界規範注意

PNG 檔案**不 commit 進 git**（per `.gitignore *.png` 規則 + business standard）：
- HTML preview 開啟時 PNG 從 disk 載入 → 正常顯示
- 但其他人 clone repo 不會拿到 PNG → 需自行重生
- 預防措施：`scripts/hooks/no_binary_commit.sh` block 未來 PNG commit

如需與 PI 分享 PPT/HTML：用 Drive/email 傳完整資料夾 + PNG，不用 git。
