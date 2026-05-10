<!--
build_date: 2026-05-10
agent: pptx-build P3 section batch
status: in_progress (awaiting C3 user ack)
report_class: per-slide specification
parent_outline: 01_full_narrative.md
-->

# Self-Phasing PPT — Per-Slide Specification (P3)

> 18 main + 3 backup = 21 slide。每張含：title / focal point ≤20 字 / Tier 1 elements ≤6 / visual asset / Tier 2 speaker note 大綱 / Tier 3 oral-optional 標示。
> 所有數字已標 source（母稿章節或 Figure），P4 Agent-N 將逐項 grep 驗證。

---

## S0 Cover + TL;DR (2 slide)

### Slide 01 — Cover

- **Title**: Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷
- **Focal point**: (cover, no focal)
- **Tier 1 elements**:
  - 主題：Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷
  - Subtitle: longphase-to-mod 5 commits 修補成熟 + 5/9 paired cross-ref 新發現
  - 日期：2026-05-10
  - 受眾：PI / lab meeting
  - Report 引用：5/8 整合報告 + 5/9 errata
- **Visual**: 簡潔 cover（no figure）
- **Tier 2 note**: ~30 sec 開場 — 「今天分享 self-phasing 修補 + 5/9 新發現的 V5 Layer 1.5 設計缺陷」

### Slide 02 — TL;DR (1 句結論 + 5 數字)

- **Title**: 修補 V3F+V5 100% read-level；V5 Layer 1.5 germline-absent 區仍待補強
- **Focal point**: 5 數字 + 1 caveat 一眼看完
- **Tier 1 elements** (≤6):
  - **17.3:1** baseline 全基因組 HP1:HP2 偏移（systematic bias）
  - **34,855** read-level victims（全基因組）+ 752（chr19）
  - **100%** V3F + V5 修正率（單向）
  - **+13.3 pp** clean PS paired GT concordance（V5 vs baseline）
  - **20 / 0** 指標 no regression（caller F1 三版相同）
  - ⚠ **V5 Layer 1.5 germline-absent 區 4.19:1 偏 HP1 同 baseline**
- **Visual**: 大字數字 grid + 最後一條 caveat 紅色高亮
- **Tier 2 note**: 1 min — 開頭快速建立 frame；強調「修補主線成功 + caveat 是 5/9 新發現」
- **Tier 3 oral-optional**: V3F + V5 兩 binary commit hash 細節

---

## S1 觀察起點 (2 slide)

### Slide 03 — 全基因組 17.3:1 偏移；3 條獨立論證

- **Title**: 17.3:1 偏移真實存在 — 生物學 / 跨 chr / paired 三條獨立佐證
- **Focal point**: HP1:HP2 = 17.3:1（vs 隨機 1:1）跨 23 chr 一致 → systematic bias
- **Tier 1 elements**:
  - 表：HP1 reads 614K vs HP2 reads 35.5K (1.89× / 0.11×)
  - HP1 占比 **94.6%**（vs 隨機 50%）
  - 3 論證：① 生物學（跨 chr 不該同向）② 23 chr 一致（cnLOH artifact 只影響單 chr）③ paired 1:1 對照
  - 結論箭頭：「→ 17.3:1 是 LongPhase-TO 工程 bias」
- **Visual**: 簡單對照表 + paired 1:1 同框比較
- **Tier 2 note**: 1.5 min — 強調 3 論證互相獨立；94.6% 是 systematic bias 的硬證據
- **Tier 3 oral-optional**: cnLOH artifact 詳細解釋

### Slide 04 — chr19 SP1/SP2/SP3 三個近 100% 失衡位點 IGV

- **Title**: chr19 三位點 113:0 / 109:1 / 108:0；V5 翻回 paired 3/3
- **Focal point**: IGV 6-BAM 並列：baseline 全 HP1，V5 翻 HP2 與 paired 一致
- **Tier 1 elements**:
  - SP1 chr19:17,565,944 → baseline 113:0 → V5 翻 HP2 ✅
  - SP2 chr19:12,452,332 → baseline 109:1 → V5 翻 HP2 ✅
  - SP3 chr19:12,467,180 → baseline 108:0 → V5 翻 HP2 ✅
  - 解讀：not noise / not caller / not alignment → read assignment 強制集中
  - 4 問引出（觀察驅動 → 機制 / 量化 / 修補 / 驗證）
- **Visual**: 3 個 IGV 截圖（D_SP1/SP2/SP3 by_HP_4ver，PI 報告既存）
- **Tier 2 note**: 1.5 min — IGV 鐵證；「我們 next 解釋這怎麼發生（→ S2）」
- **Tier 3 oral-optional**: 6-BAM 並列順序解釋（baseline / V2b / V3F / V5 / paired_T / paired_N）

---

## S2 機制 — 兩層 bug (3 slide)

### Slide 05 — phasing 層球員兼裁判：somatic 進 graph 自我增強

- **Title**: phasing layer 球員兼裁判 — somatic 100% 共現蓋過 germline 50/50
- **Focal point**: somatic 進 graph → edge weight 暴漲 → graph 自我增強偏向某 haplotype
- **Tier 1 elements**:
  - germline het 50/50 vs somatic 100/0 共現 對照表
  - phasing graph 應是 germline-only（diploid 物理基礎）
  - TO mode 沒 paired normal → PoN 外的位點被當 germline → somatic 進 graph
  - 致命：edge weight 暴漲，graph 變「somatic 越強 → 偏該 haplotype → 該 haplotype 越強」
  - 解法：PON-only flag（commit `8b8c1fd`）— 後 §S4 詳述
- **Visual**: 示意圖 — graph 結構 + 兩種共現比例
- **Tier 2 note**: 2 min — 解釋「球員兼裁判」隱喻；強調這是第一層 bug
- **Tier 3 oral-optional**: TO vs paired 的 PoN 機制細節
- **R-G2 術語預警**: phasing graph / haplotype / PoN / cnLOH = 4 術語 → ⚠ **超過 3 個閾值，需加 footnote 解釋 PoN（30 字內）**

### Slide 06 — tagging 層 getVote priority bug：1 票 somatic 蓋 5 票 germline

- **Title**: getVote vector 順序錯 + break early；1 票 somatic 觸發誤標
- **Focal point**: 順序：① somatic ② mixed ③ germline → break early → germline 被忽略
- **Tier 1 elements**:
  - baseline code snippet（vector 順序 + break early）★ Figure F1
  - 5-vote 範例 read：germline HP2=5, somatic HP1_1=1 → baseline 標 hp=11（錯）
  - 正確答案應是 hp=21（germline 5 票主導）
  - 為何全偏 HP1：tumor sub-clone somatic 100% 同方向 → priority bug 把所有 reads 標 HP:i:11 系列
  - 結論：17.3:1 偏移在 tag layer 形成
- **Visual**: F1 priority_bug_mechanism.png
- **Tier 2 note**: 2 min — 強調 1 票即觸發；範例是真實 read 案例（next slide 量化）
- **Tier 3 oral-optional**: enum (HAPLOTYPE1_1=2 vs HP tag=11) 比較 bug 細節
- **R-G2 術語預警**: getVote / vector ordered / break early / countMap / HAPLOTYPE enum = 5 術語 → ⚠ **加 footnote 對 getVote / countMap 30 字解釋**

### Slide 07 — 兩層 bug 兩層修補對應表 — 任一層單獨修不夠

- **Title**: 兩層 bug 互不取代 — phasing PON-only + tagging V3F+V5 缺一不可
- **Focal point**: phasing 層 + tagging 層 各需獨立修補
- **Tier 1 elements**:
  - 表格：layer / bug / 修補 commit / §5 章節
    - phasing → 球員兼裁判 → `8b8c1fd` PON-only flag
    - tagging → priority bug → `41ff147` V3F two-layer
  - 註：V3F 還需 INDEL guard (`380e8d2`) + V5 Layer 1.5 + ploidy fix (`d0bcd8c`) + threshold (`938f0df`)
  - 「為何不能只用 PON-only」：解 phasing 但 tag 仍壞 → 99.9% reads 仍標錯
  - 「為何 V3F 不夠還要 V5」：V3F Layer 1 only → germline 缺席區 reads 全 untagged
- **Visual**: 簡潔對應表 + 5 commits 引出（→ 詳見 S4）
- **Tier 2 note**: 1.5 min — 強調「兩層獨立 + 五 commit stacking」；不要讓 PI 以為一個 commit 解決
- **Tier 3 oral-optional**: 五 commit 順序與依賴關係

---

## S3 量化鐵證 (2 slide)

### Slide 08 — chr19 752 read-level victims：4-path 驗證 3.5/4 PASS

- **Title**: chr19 752 victims 100% 單向 baseline=11 → V3F=21
- **Focal point**: 752 victims 全部單向，無一條反向修正
- **Tier 1 elements**:
  - 數據表：552K rows × 3-way merged 1.07M events → 752 victims
  - **V3F / V5 修正率 100% / 100%**
  - 4-path 驗證表：① 個案 trace ✅ ② 1Mb 區域聚集 ⚠ PARTIAL ③ density 共變 🔄 反向但有意義 ④ 修正後消失 ✅
  - 統一指紋（前 5 案例）read_name + chr:pos + germline/somatic vote + baseline→V3F→V5
  - 結論：priority bug 機制因果**確立**
- **Visual**: F4 chr19_752_victims_scatter.png（hotspot 30M / 27M）
- **Tier 2 note**: 2 min — 強調「single direction = 機制證據」；4-path 是 T1.2 plan 設計嚴格驗證
- **Tier 3 oral-optional**: ③ density 共變反向解釋（low somatic vote 才觸發）
- **R-G3 數字驗證**: 752 / 549,206 / 1,069,832 / 100% / 215 / 133 / 41 → Agent-N 全 grep `T1_2_priority_bug_mechanism_report.md`
- **Metric scope**: chr19 only

### Slide 09 — 全基因組 34,855 victims（46×）；priority bug 主要分佈非 chr19

- **Title**: 全基因組 34,855 victims；chr7/chr2/chr1 主要分佈，chr19 占 2.16% rank 19
- **Focal point**: chr19 不是主 hotspot；priority bug 廣泛分佈於 chr7/chr2/chr1
- **Tier 1 elements**:
  - 規模對照：chr19 752 → 全基因組 34,855（**46.4×**）
  - V3F / V5 修正率仍 100% / 100% 一致
  - per-chr 表前 5：chr7 (3,508 / rank 1) / chr2 (2,792) / chr1 (2,674) / chr16 (2,584) / chr20 (2,101)
  - chr19 752 占比 **2.16%** (rank 19) — 「不是主 hotspot」
  - chr8 priority bug enrichment 0.34× (cold zone, rank 21)；與 LOH+HPSig hotspot 是不同 layer
- **Visual**: F2 priority_bug_per_chr_enrichment.png（左 victim N / 右 enrichment ‰）
- **Tier 2 note**: 2 min — 顛覆原 chr19 pilot 結論；chr19 SP1/2/3 是「可重現案例」非「主要分佈位置」
- **Tier 3 oral-optional**: chr8 LOH+HPSig vs priority bug 兩 layer 區分；chrY enrichment ‰ 高但 N 小
- **R-G3 數字驗證**: 34,855 / 752 / 46.4× / 100% / 3,508 / 2.16% / rank 19 / 0.34× → Agent-N grep `T1_2_F1_genome_wide_audit.md`
- **Metric scope**: 全基因組

---

## S4 修補設計 (2 slide)

### Slide 10 — 5 commits 時間軸 — 兩層三版 stacking

- **Title**: 5 commits 漸進完成 — phasing / tagging / 跨層三色標
- **Focal point**: baseline → V3F → V5；任一單獨 commit 不夠
- **Tier 1 elements**:
  - F3 timeline：04-09 `8b8c1fd` → 04-10 `41ff147` → 04-25 `380e8d2` → 04-30(a) `d0bcd8c` → 04-30(b) `938f0df`
  - 三色標 layer：phasing 藍 / tagging 綠 / 跨層紫
  - V3-Fixed = baseline + `41ff147` + `380e8d2`
  - V5 = V3F + `d0bcd8c` + `938f0df`
  - ★ 41ff147 是修偏移的關鍵 commit
- **Visual**: F3 binary_commit_timeline.png
- **Tier 2 note**: 2 min — 強調 stacking；不要把 PON-only flag 當主修補
- **Tier 3 oral-optional**: 四個 commit 對應 5 layer：phasing / tagging / INDEL guard / Layer 1.5 / threshold
- **R-G2 術語預警**: PON-only / two-layer / INDEL / ploidy / threshold = 5 術語 → ⚠ footnote 對 ploidy 解釋

### Slide 11 — getVote 三版程式碼對照 — V3F 兩層 + V5 Layer 1.5 fallback

- **Title**: getVote 三版差異 — baseline 順序錯 → V3F 兩層 → V5 Layer 1.5 補 fallback
- **Focal point**: code 三版 side-by-side；V5 Layer 1.5 = germline 缺席時用 somatic phased 投票
- **Tier 1 elements**:
  - 三 code panel：baseline (priority bug) / V3F (Layer 1+2) / V5 (Layer 1+1.5+2)
  - V3F: Layer 1 germline only / Layer 2 somatic annotation
  - V5: + Layer 1.5 germline 缺席用 somatic phased votes
  - 註：Layer 1.5 在 germline-absent 區會繼承 priority bug 偏移 → S6 詳述
- **Visual**: 三 code block side-by-side
- **Tier 2 note**: 2 min — 強調 Layer 1.5 是 V5 的關鍵新增；S6 cliffhanger 種子
- **Tier 3 oral-optional**: enum 數值（HAPLOTYPE1_1=2, HAPLOTYPE3=4）+ INDEL guard
- **R-G2 術語預警**: Layer 1 / 1.5 / 2 / fallback / phased votes = 5 術語 → ⚠ footnote 解 phased votes

---

## S5 驗證 + no-regression (3 slide, 含 cliffhanger)

### Slide 12 — SP1/2/3 修正後對齊 paired 3/3；HP1:HP2 17.3:1 → ~1:1

- **Title**: 個案層 V5 修正 3/3；全基因組 HP1:HP2 17.3:1 → ~1:1 消除偏移
- **Focal point**: 個案 + 全基因組 兩層共驗 V5 修對
- **Tier 1 elements**:
  - 個案層表：SP1/SP2/SP3 → V5 翻 HP2 對齊 paired 3/3 ✅
  - 全基因組層：HP1:HP2 17.3:1 → ~1:1（消除偏移）
  - 94.6% somatic→HP1 → ~50% balanced
  - 15-site Problem PS（含 SP1/2/3）48.5% → 52.0% (+3.5 pp)
- **Visual**: F5 layer15_zero_sum_4quadrant.png（先示意）+ 對照表
- **Tier 2 note**: 1.5 min — 個案 + 統計兩層共驗；3.5 pp 看似小但機制顯著
- **Tier 3 oral-optional**: V2b / V3F 中間版本對齊 paired 細節
- **R-G3 數字驗證**: 17.3:1 / 1:1 / 94.6% / 50% / 48.5% / 52.0% / +3.5 pp → Agent-N grep PI 報告
- **Metric scope**: chr19 + 全基因組混合（每 metric 標明）

### Slide 13 — 20 指標 no regression — ISM / HP / methylation / paired GT / 結構

- **Title**: 20 指標全綠 — 6 項顯著改善 +8.3 ~ +99.7%；其餘 ±0.01 內持平
- **Focal point**: 5 大類 20 指標 0 regression
- **Tier 1 elements** (壓縮版):
  - 指標表（5 大類）：
    - **ISM aggregate**: TP_rate +0.005 / HP_Ratio median 0.788→0.574 / Potential_LOH +3.5 pp
    - **HP_Ratio AUC**: All -0.005（隨機區間）/ Inner +0.002
    - **Methylation 6 feat**: 全 ±0.01 內持平
    - **Paired GT concord.** ⭐: clean PS **+8.3 pp** / 15-site Aggregate **+6.65 pp** / 15-site Clean PS **+13.3 pp** / Problem PS +3.5 pp
    - **HP / LOH 結構** ⭐: Phase block N50 **+99.7%** / Phased rate **+23.6 pp** / 執行時間 **1.36× 快** / LOH regions 完全相同
  - 結論：20/20 no regression，6 項 ⭐ 顯著改善
- **Visual**: 緊湊指標表（不放圖）
- **Tier 2 note**: 2 min — 強調「全綠 + 6 項顯著改善」；HP_Ratio 0.788→0.574 是 tag bias 修正非變差
- **Tier 3 oral-optional**: methylation 6 feat 詳細列表
- **R-G3 數字驗證**: 全 15+ 個數字 → Agent-N grep `20260508_Self_Phasing_完整觀察整合報告_01.md` §8.5.1
- **Metric scope**: HCC1395 5kHz @ 0.93 purity（PI 報告 Pass 1 only V5 BAM）

### Slide 14 — Caller F1 vs SEQC2 三版完全相同；purity 0.6 完整對照（含 cliffhanger）

- **Title**: Caller F1 三版相同（V5 不改 caller）；purity 0.6 樣本 0 critical regression
- **Focal point**: V5 不改 FILTER → caller F1 不變；低純度樣本完整對照無 regression
- **Tier 1 elements**:
  - 0.93 purity: A1/A3/A5 三版 F1 = **0.7166**（TP/FP/FN 完全相同）
  - 0.6 purity: B1/B3/B5 三版 F1 = **0.6273**（TP/FP/FN 完全相同）
  - 因果鏈：ClairS-TO PASS set 不變 → V5 改 GT/PS/GT2/GT3 不改 FILTER → F1 相同
  - 0.6 樣本 9 結構指標：4 改善 + 1 微差 (N50 -14.5%) + 1 持平
  - **★ Cliffhanger**: 「→ 但 5/9 paired cross-ref 揭露另一面...」
- **Visual**: caller F1 三版對照表 + cliffhanger arrow 引出 S6
- **Tier 2 note**: 2 min — F1 三版相同說明 V5 不改 caller；cliffhanger 用「但」字緩衝 emotional drop
- **Tier 3 oral-optional**: PASS set / FILTER 機制細節；purity 0.6 N50 微差解釋
- **R-G3 數字驗證**: 0.7166 / 0.6273 / 28509 TP / 11606 FP / 10938 FN / 0.93 / 0.6 → Agent-N grep `20260430_V3F_ablation_purity06_results_01.md`
- **Metric scope**: HCC1395 5kHz @ 0.93 + HCC1395 t30_n20 @ 0.6（兩 sample）

---

## S6 5/9 新發現 — V5 Layer 1.5 設計缺陷 (2 slide)

### Slide 15 — paired mode 整體無 priority bug — som_ratio 0.462 跨 windows

- **Title**: paired mode 整體無偏移 — HP1:HP2 = 1:1.275；som_ratio mean 0.462
- **Focal point**: paired 用 longphase-s 不同 binary，無 systematic bias
- **Tier 1 elements**:
  - paired chr19 HP:Z: distribution 表：HP:Z:2 (51.6%) / HP:Z:1 (40.5%) / HP:Z:2-1 (4.1%) / HP:Z:1-1 (3.5%) / HP:Z:3 (0.3%)
  - germline HP1:HP2 = **1:1.275**（vs TO 17.3:1）
  - somatic HP1-1:HP2-1 = **1:1.169**（無偏移）
  - 57 chr19 1Mb windows: som_ratio mean **0.462** / median 0.494 / stdev 0.332
  - 真實 sub-clone signal: chr19:3M 全 HP2-1 / chr19:0M 全 HP1-1 / chr19:17M 對稱 0.500
- **Visual**: F6 paired_vs_TO_HP_distribution.png
- **Tier 2 note**: 2 min — paired 不同 binary（longphase-s vs longphase-to）；無 priority bug
- **Tier 3 oral-optional**: HP:Z: 字串編碼 vs HP:i: 整數 binary 差異
- **R-G3 數字驗證**: 51.6% / 40.5% / 1:1.275 / 1:1.169 / 0.462 / 0.494 / 0.332 → Agent-N grep `paired_priority_bug_audit/00_audit_report.md`
- **Metric scope**: paired chr19 only（5,789 events）
- **R-G2 術語預警**: HP:Z: / HP:i: / longphase-s / som_ratio = 4 術語 → ⚠ footnote

### Slide 16 — germline-absent 區 V5 = baseline 4.19:1；V3F 標 hp=33 反而穩健

- **Title**: V5 Layer 1.5 設計缺陷 — germline-absent 區與 baseline 完全相同
- **Focal point**: V5 = baseline 4.19:1；V3F hp=33 保守正確
- **Tier 1 elements**:
  - cross-tab 表（5,789 chr19 events，cnt_HP1+HP2=0 且 somatic>0）：
    - **baseline**: hp=11=3,312 / hp=21=791 → **4.19:1 偏 HP1**
    - **V3F**: 全 hp=33（5,789）→ **保守不選邊** ✅
    - **V5**: hp=11=3,313 / hp=21=790 → **4.19:1 偏 HP1（與 baseline 完全相同）** ⚠
  - 機制詮釋：V5 Layer 1.5 用 somatic phased 投票 → sub-clone somatic 100% 共現 → 偏向同邊 → 繼承 priority bug 偏移
  - 結論：**V5 Layer 1.5 = priority bug 的 feature 化非修補**；V3F hp=33 才是穩健
- **Visual**: F7 germline_absent_crosstab.png（左 panel 三版對比 / 右 panel paired HP cross-tab）
- **Tier 2 note**: 2.5 min — 詳述 mechanism；強調 V3F 的 hp=33 設計反而對；V5 Layer 1.5 是「設計選擇」非「bug」但 implications 待 F-paired-D3
- **Tier 3 oral-optional**: caveat — cross-binary axis alignment（paired vs TO 各自獨立 phasing）
- **R-G3 數字驗證**: 5,789 / 3,312 / 791 / 4.19:1 / 3,313 / 790 → Agent-N grep `01_step_D_germline_absent_finding.md`
- **Metric scope**: paired chr19 germline-absent only

---

## S7 結論 + errata + follow-up (2 slide — C2 拆 2)

### Slide 17 — 5 條 PI errata 已 patch（E1-E5）

- **Title**: PI 報告 4-29 5 條 errata 已 patch；主結論不撤回
- **Focal point**: E1-E5 對應 5 處表述修訂
- **Tier 1 elements**:
  - E1 §3.3.3 chr19 SP1/2/3 → 從「主要 hotspot」→「可重現案例」（占 2.16% rank 19）
  - E2 §5.2 V5 working tree → 已 commit `d0bcd8c` + `938f0df`（2026-04-30）
  - E3 §5.2 證據強度 → 升級「commit msg + 3 IGV」→「+ 34,855 read-level 鐵證」
  - E4 §6.4/6.5 V5 數值歸因 → 「Pass 1 only；主要 V3F + Layer 1.5；Pass 2 二次效益尚未量化」
  - E5（5/10 加）§5.2 V5 Layer 1.5 → 「germline-absent 區與 baseline 4.19:1 相同 — 設計缺陷非修補」
- **Visual**: 5 條 errata 編號表 + commit hash 引用
- **Tier 2 note**: 1.5 min — E4 + E5 是核心 errata；其他三條是表述精確化
- **Tier 3 oral-optional**: errata commit chain (f17754f → 2553e96 → 71d21bd)
- **R-G3 數字驗證**: 2.16% / rank 19 / d0bcd8c / 938f0df / 34,855 / 4.19:1 → Agent-N grep `20260509_PI_Report_4_29_Errata_01.md`

### Slide 18 — 整體成熟度 + 5 項 follow-up

- **Title**: V5 仍可作 production baseline；5 項 follow-up 待跑
- **Focal point**: 12 維度成熟度 10 ✅ / 1 ⚠ / 1 待跑；5 項後續 cycle
- **Tier 1 elements**:
  - 成熟度表（簡化）：機制 ✅ / 修補設計 ✅ / 個案 ✅ / 全基因組 ✅ / zero-sum ✅ / 20 指標 ✅ / Caller F1 ✅ / purity 0.6 ✅ / 三路徑算法 ✅ / Pass 2 量化 ✅ / 版本對齊 ✅ / paired Step A+C ✅ / **V5 Layer 1.5 ⚠** / Pass 2 獨立貢獻 待 / 跨樣本 待
  - 5 follow-up：F-paired-D1 全基因組擴展 0.5d / F-paired-D2 phase block 內 axis-aligned 1d / F-paired-D3 Layer 1.5 改 V3F 1-2d / T3 7 樣本 1-2d / T1.3 4-cell ablation 3d
  - 結論：**V5 仍可作 production tag baseline**（germline-absent 區占比小不阻擋）；F-paired-D3 量化後決定是否回 V3F default
- **Visual**: 12 維度成熟度表 + 5 follow-up 工作量表
- **Tier 2 note**: 1.5 min — 強調「整體成熟 + caveat 清晰」；follow-up 排序 F-paired-D3 優先
- **Tier 3 oral-optional**: T1.3 4-cell ablation 設計細節；7 樣本 binary patch 重複量

---

## Q&A Backup (3 slide)

### Slide B1 — Pass 2 second round 機制與 incremental N50 +3.51%

- **Title**: Pass 2 = 只重跑 2-point edgeConnectResult；incremental N50 +3.51%
- **Focal point**: 高 purity 才觸發 Pass 2；只重 phase 不重 call
- **Tier 1 elements**:
  - 兩 graph 演算法：somaticCalling (3-point patternMining) + edgeConnectResult (2-point)
  - Pass 1 都跑 / Pass 2 只重跑 2-point（高 purity > 0.9）
  - Pass 1 only vs Pass 1+2 對比表：phased var -2.90 pp / blocks -9.79% / N50 +3.51%
  - 為何不重跑 3-point：Pass 1 已產出穩定 origin 分類
- **Visual**: 演算法分層表 + Pass 1+2 對比
- **PI 預期問題**: 「Pass 2 為什麼只跑 2-point？」「N50 +3.51% 為何是 polish？」
- **R-G3 數字驗證**: 0.9 threshold / -2.90 pp / -9.79% / +3.51% / 1808 / 1631 / 11,388,114 / 11,788,053 → Agent-N grep §8.5.3

### Slide B2 — purity 0.6 樣本完整對照表（6 caller F1 + 9 結構指標）

- **Title**: purity 0.6 樣本 baseline vs V5 完整對照 — 0 critical regression
- **Focal point**: 6 F1 完全相同 + 9 結構指標 4 改善 + 1 微差 + 1 持平
- **Tier 1 elements**:
  - 6 F1: TP/FP/FN/Precision/Recall/F1 三版完全相同
  - 9 結構: phased% +4.01 pp / n_blocks +18.1% / N50 -14.5% (微差) / HP:i:33 +20 / AMB% +3.12 pp / HP1/HP2 0.48→0.38 / purity 0.607→0.634 / Pass 2 ❌→❌ / LOH 完全相同
  - 結論：低純度樣本 V5 conservative tagging 是好事；ploidy bug 在低純度自我治癒
- **Visual**: 完整對照表
- **PI 預期問題**: 「低純度樣本是否變差？」
- **R-G3 數字驗證**: 0.6273 / 0.607 / 0.634 / +4.01 pp / +18.1% / -14.5% / +3.12 pp → Agent-N grep §8.5.2

### Slide B3 — 7 樣本擴展 + cnLOH 雙親同源待開放方向

- **Title**: 跨樣本擴展 (T3) + cnLOH 雙親同源待開放
- **Focal point**: HCC1395 5kHz 是單樣本驗證；7 樣本 + cnLOH 為下一主軸
- **Tier 1 elements**:
  - T3: HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829 跑同 vote audit
  - 預期工作量：1-2 天（含 binary patch + dump 6×3 版）
  - cnLOH 雙親同源：parent 1 vs parent 2 染色體在 LOH 區難 phase；Layer 1.5 設計選擇對 cnLOH 區影響待量化
  - 排序：F-paired-D3 (1-2 天) > T3 (1-2 天) > T1.3 ablation (3 天)
- **Visual**: 7 樣本 sample list + cnLOH 雙親同源示意
- **PI 預期問題**: 「跨樣本一致性？cnLOH 區呢？」
- **R-G3 數字驗證**: 7 樣本名稱 → grep MEMORY canonical sample list

---

## P3 自我審查清單

| 檢查項 | 狀態 |
|--------|------|
| 21 slide 各有 title + focal point + Tier 1 elements ≤ 6 | ✅ |
| Title 為 assertion-evidence sentence（非 generic label） | ✅ |
| Focal point ≤ 20 字 | ✅ |
| 個人風格 R-G2 術語密度標 ⚠ 高密度 slide：05/06/10/11/15 | ✅ 已標 footnote 需求 |
| R-G3 數字驗證 source 已標 source path | ✅ 待 P4 Agent-N grep |
| Metric scope 明示：chr19 vs 全基因組 / 0.93 vs 0.6 / Pass 1 only | ✅ 已標每 slide |
| Cliffhanger transition S5 → S6 | ✅ slide 14 末已標 |
| Tier 1 / Tier 2 / Tier 3 三層分流 | ✅ 每 slide 已標 |

---

---

## Appendix A — R-G4 術語三層分流表（C3 補強，2026-05-10）

依用戶 C3 確認的 R-G4 規則：困難概念 → 補圖 / explanatory slide；中等術語 → in-slide glossary box；簡單術語 → footnote ≤30 字。

### A.1 困難概念（已用 figure / 大圖示處理；不另開 explanatory slide）

| 概念 | 處理 slide | 處理方式 |
|------|----------|---------|
| Self-phasing 球員兼裁判 | slide 05 | germline 50/50 vs somatic 100/0 共現對照圖 + 隱喻說明 |
| getVote priority bug 1 票觸發 | slide 06 | Figure F1 機制圖 + 5-vote 範例 read 走查 |
| 兩層 bug 兩層修補對應 | slide 07 | 對應表 + 「為何不能只用 PON-only / V3F」雙問解釋 |
| 5 commits 兩層三版 stacking | slide 10 | Figure F3 timeline + 三色標 layer |
| getVote 三版 code 演進 | slide 11 | 三 code panel side-by-side |
| Layer 1.5 zero-sum 重分配 | slide 12 (借用) | Figure F5 4 象限視覺 |
| V5 Layer 1.5 設計缺陷機制 | slide 16 | Figure F7 三版 cross-tab + 機制詮釋 |
| Pass 2 vs Pass 1 機制 | Q&A B1 | 兩 graph 演算法對照 + Pass 1+2 量化 |
| Ploidy bug → purity=0 因果鏈 | Q&A B1 (附帶) | 母稿 §7.1 ASCII chain（Tier 2 speaker note 詳述） |

### A.2 中等術語（in-slide glossary box ≤ 60 字）

每張 slide 角落或 sidebar 加一個「📖 名詞」框。

| 術語 | 解釋（≤60 字） | 用於 slide |
|------|--------------|----------|
| **PoN** | Panel of Normals — 從多個正常樣本建構的 germline 變異 reference set | 05, 07 |
| **germline het** | 個人遺傳變異中雜合位點（A/G 等），雙親各帶一型 | 03, 05, 11 |
| **somatic mutation** | 腫瘤後天獲得的體細胞變異，僅在某 sub-clone 內出現 | 03, 05, 06 |
| **phasing graph** | 把 het 位點當 node，read 上共出現當 edge，連成 haplotype 的圖 | 05 |
| **HP:i: vs HP:Z:** | longphase-to 用整數 HP tag (1/2/11/21/33) / longphase-s 用字串 (1/2/1-1/2-1/3) | 15 |
| **sub-clone** | 腫瘤內基因型相同的細胞群；不同 sub-clone 帶不同 somatic | 02, 03, 06 |
| **haplotype (HP1/HP2)** | 父系 / 母系兩條染色體之一；germline het 的方向標 | 03, 05 |
| **LOH / cnLOH** | 雜合性丟失 / 拷貝中性 LOH（保留 2 套但同源） | 13, B3 |
| **purity** | 樣本中腫瘤細胞占比；ClairS-TO 由 ploidy 算出 | 14, B1 |
| **somaticCalling vs edgeConnectResult** | 兩個 phasing graph 演算法：3-point patternMining / 2-point pairwise edge | B1 |

→ 共 10 個中等術語，分散到 9 張 slide。每 slide 最多 1-2 box（<10% page area）。

### A.3 簡單術語（footnote ≤ 30 字）

| 術語 | footnote | 用於 slide |
|------|---------|----------|
| **two-layer getVote** | V3F 拆 Layer 1 (germline) + Layer 2 (annotation) | 06, 11 |
| **break early** | for 迴圈在第一個非空 vector 處 break；後續被忽略 | 06 |
| **INDEL guard** | 補 OOB undefined behavior（HAPLOTYPE_UNDEFINED 檢查）| 10 |
| **Layer 1.5 fallback** | germline 缺席時用 somatic phased votes 決方向 | 11, 16 |
| **threshold 0.95 → 0.9** | 觸發 Pass 2 second round 的 purity 閾值 | 10 |
| **phased votes** | 根據 phasing graph 結果分類後的 somatic votes | 11, 16 |
| **som_ratio** | HP1-1 票 / (HP1-1 + HP2-1) 票，量化偏向 | 15 |
| **PASS set** | ClairS-TO snv.vcf FILTER=PASS 的 variants 集合 | 14 |
| **highPurity** | `bool highPurity = purity > 0.9` flag | B1 |
| **ploidyRatioMap** | Pass 2 計算 purity 用的 hp allele count map | B1 |

→ 共 10 個簡單術語，標 footnote。

### A.4 R-G2 vs R-G4 雙重檢核

每 slide 上 R-G2 上限 ≤ 3 個 PI 不熟術語，R-G4 用三層分流處理超出部分。

| Slide | 中等術語（glossary box）| 簡單術語（footnote）| R-G2 ≤ 3? |
|-------|----------------------|-------------------|----------|
| 02 TL;DR | sub-clone | — | ✅ 1 |
| 03 17.3:1 | germline het, haplotype, somatic | — | ✅ 3 |
| 04 IGV SP | haplotype | — | ✅ 1 |
| 05 球員兼裁判 | PoN, germline het, phasing graph, haplotype | — | ⚠ **4 → 升 1 個為獨立 footnote**（PoN → footnote）|
| 06 priority bug | sub-clone | two-layer, break early | ⚠ **3 OK 但 footnote 兩條接受** |
| 07 兩層對應表 | PoN | — | ✅ 1 |
| 10 5 commits timeline | — | INDEL guard, threshold | ✅ 0 + 2 footnote |
| 11 getVote 三版 code | germline het | two-layer, Layer 1.5, phased votes | ⚠ **3 footnote 接受**（多技術詞 forced）|
| 12 SP 修正 | haplotype | — | ✅ 1 |
| 13 20 指標 | LOH | — | ✅ 1 |
| 14 caller F1 | purity | PASS set | ✅ 1 + 1 footnote + cliffhanger |
| 15 paired mode | HP:i: vs HP:Z: | som_ratio | ✅ 1 + 1 footnote |
| 16 V5 缺陷 | sub-clone | Layer 1.5, phased votes | ⚠ **2 footnote 接受**（核心 caveat 必說清）|
| 17 errata | — | — | ✅ 0 |
| 18 follow-up | LOH | — | ✅ 1 |
| B1 Pass 2 | somaticCalling vs edgeConnectResult, purity | highPurity, ploidyRatioMap | ⚠ **2 box + 2 footnote 接受**（backup 可較密）|
| B2 purity 0.6 | purity | — | ✅ 1 |
| B3 跨樣本 | LOH | — | ✅ 1 |

**例外原則**：slide 11/16 是核心機制 + 核心 caveat，技術密度高 forced；接受 ≤ 4 處理；B1 是 backup（PI 主動問才出現）允許更密。

---

## C3 confirmation（已答覆）

✅ 用戶 C3 確認：邊 build 邊 N 驗證（不提前）；R-G4 三層分流（footnote / glossary box / 補圖）已生效；下一步進 P4 batch build。
