# 基礎概念範本貯備計畫 + 改進 backlog（保留記錄）

> 來源：workflow `wf_76c92c4c-0cf`（5 agents，Cell/Nature/Science 頂尖圖示研究）。**全期刊規格 = L3（web 通例/作者指南，非本專案真值）**。範本數值一律 schematic(synthetic) 或 data_ref 真值（§13-A renderer 缺 verified 即 refuse）。
> 目的：記錄「該建哪些基礎概念範本 + 依常用度排序 + 頂尖對齊改進」，供逐批實作。✅=已建 / ⬜=待建。

## 基礎概念 catalog（依常用度）

| Pri | 概念 | canonical 畫法（L3）| primitive | 狀態 |
|----|------|------|-----------|------|
| 高 | gene model（exon/intron/UTR）| exon `<rect>`(寬∝長)+intron `<line>`+strand 箭頭 | NEW `gene_model_track` | ✅ B1 |
| 高 | coverage/depth | 沿座標 area，y=depth；CNV 疊 ploidy 線 | NEW `coverage_track` | ✅ B2 |
| 高 | DNA 甲基化 @ CpG（單軌）| lollipop 實心=甲基/空心=未甲基 | NEW `methyl_lollipop_track` | ✅ B3 |
| 高 | germline vs somatic 二分 | 藍圓 germline / 橘圓 somatic + legend | NEW `variant_class_legendcard` | ✅ B7 |
| 高 | phasing/haplotype HP1/HP2 | 混 reads → 依 het 分 HP1/HP2 軌 | 既有 `hap_split_track`（germline-only）| ✅ B4 t1 |
| 高 | read×CpG → β → Δβ（入門）| read 狀態 → 群率 β → 兩群 Δβ 方向 | 既有 `read_cpg_matrix`+`beta_bar` | ✅ B6 t3 |
| 高 | haplotagging（read→HP，getVote）| read 跨 het allele 投票 → tag HP | 既有 hap_split_track + stage-flow | ⬜ B5/B9 |
| 高 | LOH | ideogram + interval 紫；BAF 0.5 裂雙帶 | 既有 `loh_ideogram` | ✅（U7）|
| 高 | ASM/ASE（cis vs drift）| HP1 軌 vs HP2 軌 β 差；normal triplet 拆 cis/drift | 既有 `normal_cis_triplet` | ✅（C1）|
| 高 | CNV（BAF+LogR）| 上 LogR(0=2copy)+下 BAF(0/.5/1) | NEW `baf_logr_track` | ⬜ B9 |
| 中 | SNV/indel 形態 | SNV 彩 tick / indel caret+gap | 既有 `igv_pileup`（擴 indel）| ⬜ |
| 中 | VAF/AF | alt/total 比例 → 數值 | NEW `vaf_panel` | ⬜ |
| 中 | Tumor-Normal vs TO | 兩欄 paired vs TO+PoN | NEW `tumor_normal_setup_card` | ⬜ |
| 中 | PoN filter | candidate − PoN = filtered（算式）| NEW `pon_filter_flow` | ⬜ |
| 中 | clonal evolution | fishplot：x=時間，nested polygon | NEW `fishplot_panel` | ⬜ |
| 低 | bisulfite vs nanopore modBAM | 兩欄流程對比 | NEW `methyl_platform_compare` | ⬜ |
| — | mechanism before→after（Cell GA 三段）| initial→event→consequence + numbered step | NEW `numbered_step` + stage-flow | ⬜ B8 |

## 第一批已建（2026-06-07）
B1 u8_gene_model_track · B2 u9_coverage_track · B3 u10_methyl_lollipop_track · B4 t1_phasing_101 · B6 t3_read_cpg_beta_dbeta · B7 u11_variant_class_legendcard（+ 既有 C1/U1-U7）。

## 改進 backlog（對齊頂尖；逐批做）
**renderer**：A1 panel 自動編號(a/b/c 8-9pt bold upright) · A2 export_target 多尺寸(slide_16x9 / journal_single_90mm / journal_double_180mm / cell_ga_square) · A3 字級集中 `shared.type_scale` · A4 color_scale `mode: pi_slide|journal_a11y`(Wong) + 改 hap_split_track L216-217 硬編 fallback · A5 stroke 下限 0.5pt · A6 全域 font-family sans-serif · A7 reading_flow + numbered-step/arrow connector primitive。
**lint（加 C8-C12）**：C8 CVD 安全(red-green 無形狀冗餘 WARN；甲基 filled/open enforce) · C9 grayscale 對比 · C10 panel-label(multi-primitive 無 letter WARN；italic FAIL) · C11 single take-home(每圖一 thesis) · C12 arrow legend + no-decoration(drop-shadow/gradient/裝飾格線 FAIL)。
**已對齊頂尖（保留）**：藍/橘=最安全 CVD 對 · SVG-vector-first=期刊偏好 · category-by-fill+甲基形狀冗餘=「不單靠顏色」· 4 色 legend<6 · 甲基 red/blue ramp 與 haplotype 藍橘軸正交 · generate→lint→iterate=PLOS 十規 Rule 10。
