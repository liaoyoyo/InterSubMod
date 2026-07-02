<!--
建立時間: 2026-06-25
類型: methodology — chr17:48360161 subclone worked example 完整解釋敘述（論文可用 canonical 例）
狀態: in_progress（單樣本 HCC1395 ⭐3）
build_branch: docs/method-comparison-ism-external-202606
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/chr17_complete_data.json, docs/methodology/_assets/20260618_subcluster_pilot/chr17_tree_data.json, docs/methodology/_assets/20260618_subcluster_pilot/p2_linkage.json
-->

# chr17:48360161 — 一個 subclone 從「候選」到「完整重建」的 worked example

> **敘述框架**：Discovery-Narrative（發現歷程）— 用一個位點把整套 subclonal reconstruction 的方法、陷阱、與洞見串成可重述的故事。**這是論文方法章的 canonical 範例。** HCC1395 單樣本 ⭐3；數字全來自 data_sources JSON。

## §0 為什麼選這個位點

chr17:48360161 把整套方法的**每一環**都示範了一遍：① 甲基如何「提出候選」、② 為何甲基單獨不能定 subclone、③ sSNV 連鎖如何「非循環確認」、④ 一個矛盾如何揭露**漏掉的 somatic**、⑤ 完整克隆樹如何重建、⑥ 甲基最終扮演的正確角色。讀懂這一個，就讀懂整個 pipeline。

## §1 起點：甲基 cis-test 篩出的「候選」

全基因組 34,736 位點 → 甲基結構篩出 1,139 候選 → normal-anchored cis-test（NACT）逐層洗掉循環假象 → 收斂到 9 個「tumor-specific、非循環」候選。**chr17:48360161 是其中之一**。但此時我們只知道「這裡的甲基結構不是 germline cis、像是 tumor 特有」——**還不知道它是不是真 subclone**（甲基單獨永遠不能回答，因為甲基切群有 double-dip 循環風險）。

## §2 第一層遺傳錨：sSNV read-level 連鎖

要非循環地確認，需要**遺傳證據**：同一條 long read 上**多個 somatic SNV 的共現**。此窗內 filtered-TP 集有 3 個 somatic：48362515(β1)、48365089(α)、48365161(β2)。逐 read 取等位 → 共現 2×2：

- `α × β2`：ALT-REF 20、ALT-ALT 19、**REF-ALT 0** → β2 的 ALT **從不在沒有 α 的 read 上出現** → **β2 嵌套在 α 內**（α 祖先、β2 衍生）。
- `β1 × β2`：只有 REF-REF 與 ALT-ALT → β1、β2 **完美共連**（同屬一個衍生事件）。

→ 第一層結論：α 是較早的 somatic（廣），其上累積了 β（β1+β2）形成一個**後代 subclone**。這是 nested（線性）演化，**有方向**。

## §3 一個矛盾：HP1-1 但三個 sSNV 全 REF

把 read 按 α/β 基因型分群後，有一群 read **在我們 3 個 sSNV 全是 REF**（看似「無 somatic 的祖先」）。但檢查 longphase-S 的 HP tag——**這群裡有 5 條被標 HP1-1**。

🔑 **關鍵推理**：LongPhase-S（v1.7.3 `--somaticMode`）的 HP1-1 定義 = 「在 germline hap1 上、被 phase-block 內**任一 somatic SNV** 標記的 read」。所以**「HP1-1 但我們已知 sSNV 全 REF」⟹ 這條 read 必然帶著一個我們沒看到的 somatic SNV**。HP tag 比我們手上的 VCF「知道更多」。

## §4 找出漏掉的 somatic γ → 完整克隆樹

順著矛盾追：這 5 條 read 的 phase-set PS=**46496608**（phase-block ~1.9Mb），read 自身跨度 48326-48381kb。在其 span 內、filtered-TP 集**沒有**第 4 個 somatic——但 longphase 的輸入集 `HCC1395_methyl_PASS` 與 filtered-FP 集裡有 **48357368(C→T)，叫它 γ**。逐 read 定型：

- **那 5 條 HP1-1 read 在 γ 是 ALT**；γ 的 normal = **35/0 REF**、tumor VAF 0.18 → **ONT 偵測候選**（somatic-like，但 SEQC2 未收 HighConf/superSet ＝缺金標準確認，**非「真 somatic」定論**）。
- `γ × α`：**ALT-ALT 0** → γ 與 α **互斥** → γ 與 α 是**兩個 sibling 分支**（不是祖先-後代）。

🔴 **γ 不在 SEQC2 HighConf/superSet**（落在 callable HC region 內＝缺確認，**非 SEQC2 主動判 FP**；它落到我們 pipeline 的 filtered-FP 集），所以只用 filtered-TP 集時我們**漏掉它**，把 γ subclone 誤標成「ancestral root」。補回 γ 後，**完整局部克隆樹**浮現：

```
            germline hap1（ancestral，6 reads，無 somatic）
            ├── γ subclone（+48357368）      5 reads   ┐ sibling
            └── α subclone（+48365089）     39 reads   ┘ 兩分支互斥
                ├── L1 α-only               20 reads
                └── L2 α+β（+48362515,+48365161 共連）  19 reads  ← nested 後代
```

**4 個細胞群**，全 4 位點 normal=REF（α/β1/β2 = SEQC2 HighConf 確認；γ = ONT 候選，SEQC2 未收）。這是 read-level 連鎖**直接重建**的乾淨局部 subclonal phylogeny（既有 sibling 分支、也有 nested 深度）。

## §5 甲基的正確角色：characterize 已驗 lineage（非循環）

lineage 一旦由 **sSNV（遺傳）定義好**，再看甲基就**不循環**了。實測：
- **L1 祖先 vs L2 後代：16 個 CpG 甲基顯著不同**（|Δβ|≥0.2）。
- **per-CpG 可歸因到驅動者**：23 個 CpG 關連 α（祖先 somatic 的 ASM）、6 個關連 L1→L2 轉變（後代特有甲基）、少數關連 β。

🔴 **對照——甲基單獨做不到**：ISM 的無監督甲基 clustering（coarse_label）只把 reads 分成 1-2(8 條，L0/γ 為主) 與 1-1(42 條，**L1+L2 混**)；BERNOULLI 距離 UPGMA 樹也一樣——**只切得出「有沒有 somatic」，切不出 L1↔L2 subclone**。→ 證實：甲基**提出候選 + 刻畫已驗 lineage**，但**細分 subclone 必須靠遺傳錨**。

## §6 三個方法學教訓（可直接寫進論文 Methods/Discussion）

1. **甲基 ≠ subclone caller**：genotype-同質群內的甲基結構永遠 double-dip；甲基的價值是 characterize/corroborate，subclone 的非循環確認來自**同 germline-HP 上 ≥2 somatic SNV 的 read 共現**。
2. **「HP-tag 但已知 sSNV 全 REF」= 偵測漏掉 somatic 的訊號**；完整克隆重建應用 **longphase 輸入集（methyl_PASS）+ normal 確認**，而非只用 SEQC2-filtered TP（TP filter 會移除 clonally-informative 的候選，如 γ——ONT 偵測 somatic-like 但 SEQC2 未收 HighConf）。
3. **兩個尺度的克隆結構**：HP tag 在 **phase-block（~Mb）** 尺度標 lineage（靠 phasing）；read 共現只在 **≤read-長（~10-50kb）** 尺度直接連。兩者互補，但 phase-block 跨區無法用 read 連（reconstruction GAP）。

## §7 對論文的意義

chr17:48360161 是論文標題「**Subclonal reconstruction using somatic haplotagging and methylation profiles**」的具體體現：
- **somatic haplotagging（sSNV 連鎖 + LongPhase-S HP tag）= 重建骨幹**（非循環定出 γ/α/L1/L2 克隆群）；
- **methylation profiles = characterize**（刻畫各克隆群的表觀差異，非偵測驅動）。

🔴 **誠實邊界**：⭐3 單樣本；這是**局部**克隆階層（此 phase-block 內），非 genome-wide clone tree；γ 等部分 somatic 在 single-cell/multi-region 才是 confirmation 黃金標準。但作為 read-level 分子證據，這是一個**乾淨、完整、可重述**的 subclone 範例。

> 紅線：⭐3 單樣本 HCC1395；subclone = 局部 ≥2 lineage 分子證據非 genome-wide tree；甲基 corroborate 非偵測；§13 數字全 JSON 溯源。相關 HTML：`InterSubMod/docs/methodology/20260625_chr17_complete_subclone_example_01.standalone.html`（完整樹）+ `..._chr17_48360161_subclone_workstation_01.standalone.html`（距離/UPGMA/歸因）。
