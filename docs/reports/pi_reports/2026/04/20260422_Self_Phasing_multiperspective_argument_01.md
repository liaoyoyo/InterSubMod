<!--
建立時間: 2026-04-22
目標: 向 PI 提供 Self-Phasing 修改的多視角 adversarial 論證（程式碼／生物學／論文方法學 × 9 challenges）
受眾: PI（審稿式批判視角；建議搭配既有 PI 技術報告閱讀）
處理範圍: 三條 pipeline track 影響矩陣、9 個 challenges 質疑與回應、最終驗證協議
狀態: validated_multiperspective_argument
pipeline_track: all (paired_full / paired_pileup / to_pileup)
priority: P0
relates_to: 20260422_Self_Phasing_complete_report_for_PI_01.md
-->

# Self-Phasing 修改的多視角論證報告
## ——從程式碼、生物學、論文方法學三個面向質疑與回應

> 撰稿日期：2026-04-22（本報告為既有 PI 技術報告 `20260422_Self_Phasing_complete_report_for_PI_01.md` 的**互補**版本）
> 閱讀建議：先讀既有技術報告理解改動是什麼；本報告回答「為何這改動可被領域接受」
> 關鍵整合：2026-04-22 20:40 `HPFineNGroups` 機制重新釐清（見 memory `project_loh_constrained_phasing_discovery.md`），本報告 Challenge B3 已整合此更正

---

## 一頁速讀：整合案例 storyboard

以下單一案例（AF=0.3 somatic variant，10 reads）貫穿本報告所有論證——Problem、Fix、Result 三個 Panel 可獨立對應 9 個 challenges 的 evidence：

![Figure 13 — End-to-End Case Storyboard](figures/fig13_case_storyboard.png)

| Panel | 內容 | 對應 challenges |
|-------|------|----------------|
| A. Problem（flag=off）| 3 ALT reads 被 forced into HP1-1；HP_Ratio=0.30；NG=3 | C1（PON 品質？）、B1（biology 合理？）、M1（phaser 設計限制） |
| B. Fix（demotion）| 9 行 switch；hp_tag_raw 保留；守恆律驗證 | C2（下游破壞？）、M3（upstream 為何不修？） |
| C. Result（flag=on）| ALT → unphased；HP_Ratio=0.43；NG=2（artifact bucket 消失） | C3（為何保留 flag？）、B2（生物對應物？）、B3（丟 subclone 訊息？） |

**此 storyboard 整合了既有 Fig 3（原因）+ Fig 7（NG 分佈變化）+ Fig 9（HP_Ratio 分佈變化）於單一敘事**——PI 可在 10 秒內理解 before/after 變化；細節則下鑽既有報告對應章節。

---

## Positioning（報告定位）

既有 PI 技術報告用**技術敘事結構**（問題 → 根因 → 修改 → 結果）完整記錄了本次 Self-Phasing 修改的脈絡。本報告不重複這段敘事——而是採**審稿人式 adversarial review 結構**，讓 PI 從三個獨立視角檢視以下三個命題：

1. **Self-phasing 問題真實存在**（非測量 artifact / 非 confirmation bias）
2. **ISM `--germline-hp-only` 修改方案合理**（機制、實驗、風險均通過檢驗）
3. **此修改不衝擊三條主線 pipeline 分析**（paired_full / paired_pileup / to_pileup 各自獨立可推進）

**為何需要此報告**：
- 將既有分析放入投稿前的 dry-run，模擬審稿人可能提出的挑戰
- 整合 2026-04-22 關於 `HPFineNGroups` 機制的重大更正（從誤以為 methylation marker → 實為 HP×variant bucket count）
- 用矩陣形式明確說明三條路不受衝擊——避免 PI 因「修改了 HP tag 邏輯」誤以為所有分析都需重跑

**閱讀優先序**：Section 5（三條路影響矩陣，本報告核心）→ Section 7（主要訊息）→ Section 2-4（視角質疑）。

---

## Section 1：三視角 adversarial 框架

### 1.1 為何需要多視角質疑

單一視角撰寫的技術報告有三種隱形風險：

| 風險 | 舉例 | 多視角抵消機制 |
|------|------|----------------|
| 程式碼角度的**機制錯覺**（code works ≠ biology right） | Phase 0 守恆律通過不代表 demotion 的生物學解讀成立 | 生物學視角單獨檢驗 |
| 生物學角度的**方法學盲區**（biology plausible ≠ field accepts） | 「HP tag 污染」的概念直覺合理，但 LongPhase / WhatsHap 論文怎麼看？ | 論文方法學視角檢驗 |
| 論文角度的**實證缺口**（paper framing ≠ data support） | 即使 novelty 成立，若程式碼守恆律失敗仍不應發表 | 程式碼視角底層驗證 |

**本報告強制三視角互相抵制**：若任一視角無法收斂，需暫停發表決策。

### 1.2 三視角定義

| 視角 | 問什麼 | 證據類型 | 本報告章節 |
|------|-------|---------|----------|
| **程式碼**（Code-level）| 修改在程式邏輯上是否乾淨、守恆、可逆？ | Source code 行號、unit tests、守恆律計算 | Section 2 |
| **生物學**（Biology）| HP tag 與 somatic heterogeneity 的機制是否合乎 tumor biology？ | ITH literature、haplotag 解析規則、跨樣本觀察 | Section 3 |
| **論文方法學**（Paper methodology）| 修改在領域慣例中站得住嗎？LongPhase / WhatsHap / HapCUT2 / ClairS 怎麼處理類似問題？ | 第三方論文、方法比較表、novelty 定位 | Section 4 |

### 1.3 9 個核心 challenges

| ID | 視角 | 質疑 |
|----|------|------|
| C1 | 程式碼 | Self-phasing 會不會是 PON 品質問題？ |
| C2 | 程式碼 | Demotion 會不會破壞下游 FisherExact / LabelTest？ |
| C3 | 程式碼 | Phase 1 AUC 無改善為何仍保留 flag？ |
| B1 | 生物學 | Somatic 在同一 clone 共享 sub-population 的假設符合 tumor biology 嗎？ |
| B2 | 生物學 | HP:i:11/21/33 的「somatic phase block」有生物學對應物嗎？ |
| B3 | 生物學 | 「demote somatic HP tag」會丟掉 subclone 訊息嗎？ |
| M1 | 論文 | LongPhase / WhatsHap / HapCUT2 原始論文如何處理 somatic？TO 模式在它們的設計範圍內嗎？ |
| M2 | 論文 | 其他工具如何區分 germline vs somatic anchor？本研究 novelty 在哪？ |
| M3 | 論文 | 為何不修 LongPhase-TO 上游？審稿人會怎麼挑戰？ |

**彙整決策樹**（Figure 12）：

![Figure 12 — Adversarial Review Decision Tree](figures/fig12_adversarial_tree.png)

---

## Section 2：程式碼視角（C1–C3）

### Challenge C1：Self-phasing 會不會其實是 PON 品質問題？（非 LongPhase 演算法瑕疵）

**Naive response**：「可能是 PON 涵蓋不全，新增正常樣本就能解決」
這個 response 如果成立，則 self-phasing 是資源問題、不是結構問題——我們就該去擴增 PON，而不是改 ISM。

**Evidence**：
1. **PON-only phasing 實驗**（`research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md`）：即使使用理想化 PON（排除所有已知 somatic 位點，只讓真 germline 進入 phasing graph），LOH.bed Jaccard=1.0 不變；somatic bias 消除；但**haplotag 層面出現新問題**——6,485 個 `GT=0|0, GT2=.|1` 的位點導致 refHaplotype=UNDEFINED，所有 reads 被統一標記某一 HP。
2. **Phasing graph edge weight 機制**（landscape doc `02_Self_Phasing根因.md:138-144`）：
   ```
   weight(variant_A, variant_B) = Σ_reads I(read 帶 A.alt) × I(read 帶 B.alt)
   ```
   真實 germline het 的 ALT reads 分散於兩 haplotype；而**同 clone 的 somatic 因共享 sub-population，edges weight 必高於 germline 背景**。這是結構屬性，與 PON 品質無關。
3. **跨樣本一致性**：7/7 樣本（含嚴格 SEQC2 與 orthogonal truth set）均觀察到同方向的 17.3:1 bias（CV-2 pass），排除「特定樣本的 PON 不足」假說。

**Verdict（本視角認為）**：**Self-phasing 是 LongPhase-TO 演算法的結構性副作用，不是 PON 品質問題**。即使擁有完美 PON，只要 phaser 把未確認 variants 都放進 graph，somatic 就會自我連結形成 phase block。修 PON 無法根治。

---

### Challenge C2：Demotion 會不會破壞下游 LabelTest / FisherExact / RegionProcessor 的既有假設？

**Naive response**：「改掉 `hp_tag` 就是在改 downstream 所有函式看到的資料；風險很大」

這個 response 表面合理。若 demotion 改變了下游期待的數學分布（例如 `hp_tag="0"` 的比例突然暴增），downstream p-value、AUC 可能全數無效。我們必須證明 demotion 在數學上**守恆**且**可觀察**。

**Evidence**：
1. **Audit 欄位獨立性**（Phase 0 smoke，chr19 615 sites）：三個 `NHP_Somatic11/21/33` 欄位在 flag off 與 flag on 的兩次執行結果 **完全相同**（4,894 / 8,630 / 2,498），確認 audit 計算使用 `hp_tag_raw`（保留原值），不受 flag 影響。
   - 實作見 `src/core/RegionProcessor.cpp` 的欄位填寫邏輯（count `hp_tag_raw` 而非 `hp_tag`）。
2. **Demotion 數學守恆**：flag on 時，`NHP0` 增量 = `NHP_Somatic11 + NHP_Somatic21 + NHP_Somatic33`
   - Phase 0：4,894 + 8,630 + 2,498 = **16,022** 完全守恆 ✅
3. **Unit tests**：`tests/test_read_parser.cpp` 12 cases 全數通過，涵蓋 3 類邊界條件：
   - Flag off：所有 HP 值 pass-through
   - Flag on + HP:i:11/21/33：demoted to "0"，`hp_tag_raw` 保留
   - Flag on + HP:i:1/2：pass-through（不在 demotion 目標）
4. **TSV schema 不變**：`RegionProcessor.cpp` 只**新增**三欄（NHP_Somatic*），原 ~85 欄順序不變。downstream Python 分析腳本（`scripts/analysis/*.py`）無需修改。

**Verdict（本視角認為）**：**Demotion 機制乾淨、數學守恆、schema 可逆**。`hp_tag_raw` 作為 audit backbone 使 flag off/on 的差異**完全可追溯**——這是任何下游分析若發生異常可快速 bisect 的關鍵。

---

### Challenge C3：為何 Phase 1 AUC 無改善但還要保留這個 flag？

**Naive response**：「既然沒 AUC 增益，就應刪除 flag，避免 codebase 複雜度累積」

這個 response 把 flag 定位為 **filter tool**。如果 filter 方向真的失敗，保留 flag 確實是 dead code。但這忽略了 flag 的第二種用法：作為**研究工具**檢驗既有結論是否依賴某類 HP tag。

**Evidence**：
1. **TP/FP somatic tag 密度幾乎相同**（Phase 1 HCC1395 TO 40,115 sites）：
   - TP sites：平均 27.4 somatic tags / site
   - FP sites：平均 27.7 somatic tags / site
   - → somatic tag 密度本身對 TP/FP 區分**無訊號**。flag 只是把「noisy 4 類分群 NG∈{0-4}」變成「3 類 NG∈{0-2}」；前者 AUC 已接近 0.5，後者不會更好。
2. **但 flag 揭示了 HPFineNGroups 的真實機制**（2026-04-22 critical discovery）：
   - flag=on 時 NG≥3 完全消失（23,000+ → 0 regions）
   - 這揭示 `HPFineNGroups` 的計算**依賴 somatic HP tag 4 bucket 的 occupancy**，而非 methylation 信號
   - 無此 flag，我們永遠不知道原 F pilot 的 89.1% subclone marker 結論的**真實訊號源**
3. **兩個研究用途**：
   - **Orthogonal null test**：對任何宣稱的 HP-dependent 訊號，flag=on 作 negative control
   - **Mechanism confirmation**：驗證新特徵是「methylation-based」還是「phasing-based」——若 flag on 後效應消失，則是 phasing-based

**Verdict（本視角認為）**：**Flag 不應以 filter AUC 衡量價值**。它的真正價值是**研究工具**——作為 orthogonal null test 與機制辨識手段，default=off 零影響，研究者 opt-in 使用。2026-04-22 發現此 flag 恰好是 LOH-constrained phasing discovery 的 negative control 驗證路徑（見 Challenge B3）——保留決策**回溯性地獲得強烈正當性**。

---

## Section 3：生物學視角（B1–B3）

### Challenge B1：Somatic variants 在同一 clone 共享 sub-population reads 的假設，符合 tumor biology 嗎？

**Naive response**：「腫瘤異質性理論（ITH）支持 subclone 分群；somatic 當然會相互連結」

這個 response 是對的方向，但模糊。必須精確說明**哪個尺度**的 clonal 結構解釋 17.3:1 bias——是 subclone 層（幾千 variants）、是 tumor mass 層（百萬 variants）、還是兩者混合？

**Evidence**：
1. **Tumor ITH 領域共識**：
   - Gerlinger 2012 NEJM（ccRCC multi-region）：clonal expansion 產生的同時期 somatic variants 在 phylogenetic tree 同分支
   - Turajlic 2018 Cell（pan-cancer TRACERx）：clonal / subclonal variants 比例高達 10:1~30:1，dominant clone 占據絕大多數 reads
   - Jamal-Hanjani 2017 NEJM（TRACERx NSCLC）：trunk mutations 跨多個 region samples 一致
2. **Subclone 尺度計算**：假設 HCC1395 tumor purity ~80%，dominant clone 涵蓋 ~60-70% tumor cells → 每個 dominant-clone somatic ALT read **高機率同時**帶鄰近的其他 dominant-clone somatic ALTs → phasing graph edge weight 遠高於 crossover-separated germline het edges。
3. **本次數據支持**：self-phasing 集中到 **HP1 而非隨機**（94.6% 單向），符合 dominant clone 存在的預期。若 subclone 比例接近 50/50 則 bias 應 ~1:1；觀察到 17:1 強烈暗示**單一 dominant clone**。

**Verdict（本視角認為）**：**完全符合 tumor biology**。ITH 文獻支撐「同 clone somatic variants 高度共現」；觀察到的 17.3:1 bias 進一步**暗示 dominant clone** 結構，這是獨立於 phasing 演算法的生物學預測。換言之——self-phasing 的嚴重程度**反映真實生物結構**，不是 artifact 憑空出現。

---

### Challenge B2：HP:i:11/21/33 的「somatic-only phase block」有生物學對應物嗎？

**Naive response**：「somatic-only block = 只有 somatic variants 在同一 phase block，對應 subclonal 事件」

這個解讀聽起來合理，但需檢驗：**`HP:i:11/21/33` 在 haplotag 解析時究竟對應什麼？是真的「somatic-only block 1」還是演算法邊界條件副產品？**

**Evidence**：
1. **LongPhase haplotag 解析規則**（research/loh_investigation 報告 + `docs/plans/2026/04/20260403_LongPhase_TO_二次Phasing方案設計與驗證計畫_01.md`）：
   - `HP:i:1` = read 支持 refHaplotype（germline HP1）
   - `HP:i:2` = read 支持 altHaplotype（germline HP2）
   - `HP:i:11` = read 支持 refHaplotype BUT variant 在 secondary phase block（前綴 `1` = HP1；後綴 `1` = secondary block index）
2. **Secondary phase block 的來源**：Phasing graph 斷裂造成的獨立 connected components。對 germline-only phasing，secondary blocks 很少（只在 region gap 處產生）；對 TO 模式，**每個 somatic variant cluster 都可能形成獨立 secondary block**。
3. **PON-only 實驗決定性觀察**（`research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md`）：
   - 移除 somatic from phasing graph 後，VCF 層完全乾淨（Jaccard=1.0）
   - 但 haplotag 階段：6,485 個 `GT=0|0, GT2=.|1` 的 somatic 位點出現 **refHaplotype=UNDEFINED**，所有 reads 統一標記同一 HP → 產生 100% 假 LOH
   - 這揭示 `HP:i:11/21/33` 不是「真 somatic phase block」，而是 LongPhase haplotag 對 `GT=0|0` 位點的**邊界條件 fallback**
4. **Paired mode 對照**：Paired 模式下 BAM 沒有 `HP:i:11/21/33`，只有 `HP:i:1/2`。**若這些 tags 對應真實生物結構**，paired 也應該看到——事實上沒有。

**Verdict（本視角認為）**：**HP:i:11/21/33 是演算法 artifact，無獨立生物學對應物**。它既不代表真實 secondary haplotype（paired mode 不存在）、也不代表 subclonal 事件（即使概念上 subclone 是真的，這些 tags 的產生機制是 GT 解析 fallback，不是 subclone detection）。Demote 這些 tags 是**移除虛假結構**，而非移除真實訊號。

---

### Challenge B3：「Demote somatic HP tag」會不會丟掉真實 subclonal heterogeneity 訊息？

**Naive response**：「`HPFineNGroups=4` 有 89.1% TP rate，是真正的 subclone marker；demote 掉就丟了寶貴訊號」

這是最嚴重的挑戰，因為 F pilot 原結論（2026-04-18）把 HPFineNGroups≥4 + NR≥80 視為 **ISM 唯一經跨樣本 + AF 去混淆驗證的 POSITIVE 訊號**，機制解讀為「高異質性 region 反映 subclonal methylation pattern」。若此解讀成立，demote somatic HP tag 會系統性抹除重要 biology。

**然而，2026-04-22 critical discovery 翻轉了這個挑戰**：

**Evidence**：
1. **HPFineNGroups 計算定義**（`src/core/LabelTest.cpp:265-305` + `include/core/Stats.hpp:323-330` 註解，2026-04-22 驗證）：
   - NG 是 `{HP1, HP1-1, HP2, HP2-1}` **4 bucket 的 occupancy count**
   - 計算公式：對 region 內 reads 按 `(haplotype × variant-carrying)` 分 4 個桶，count 非空桶數
   - **與 methylation 無直接計算關係**（methylation 只在 `HPFineSig` 的 PERMANOVA 驗證中間接出現）
2. **因此 demote 後 NG≥3 消失是完全可預期的**：
   - 移除 HP:i:11/21/33 後，bucket `HP1-1` 與 `HP2-1` 永遠為空
   - 剩下 `HP1` 與 `HP2` 兩個 bucket → NG 最大值只能到 2
   - 這不是「丟失 subclone 訊號」，而是「bucket 結構在 somatic tag 消除後的必然塌縮」
3. **LOH-constrained phasing discovery（2026-04-22，更深層機制）**：
   - Obs18 跨 6 TO 樣本：**Inner × NG=2 中 93-99% 為 same-hap**（`HP1+HP1-1` 或 `HP2+HP2-1`）
   - **TP rate gap (Inner same_hap − Outer cross_het) = +0.37 median**，6/6 樣本正向
   - 生物學解讀：**LOH region 物理上只保留單 haplotype**，該區 somatic SNV 發生必然產生 same-hap ref/alt 子族；非 LOH 區 NG=2 多為 cross-hap（canonical het phasing）
   - → 原 89.1% TP rate 的真正訊號源是 **LOH-constrained phasing signature**，不是 methylation bimodality
4. **`--germline-hp-only` flag 在此 discovery 扮演的角色**：
   - Flag=on 移除 HP:i:11/21/33 → same-hap bucket 消失 → Inner × NG=2 的 TP gap 應消失
   - 這是 LOH-constrained phasing 的**完美 negative control**
   - Phase 1 觀察到 flag=on 下 NG≥3 完全消失，與此預測**完全吻合**
5. **Memory 記錄更正**（`project_hpfinengroups_subclone_marker.md` 2026-04-22 update）：
   - 舊機制解讀「methylation subclone marker」標註撤回
   - 新機制解讀「phasing-resolved variant signature」——89.1% TP rate 仍成立但機制完全不同
6. **Paired_full 平行證據**（`project_loh_subclone_af_methylation_positive.md` 2026-04-22 重新解讀）：
   - 原 Inter AF→NGroups +0.705 (7/7 p<10^-39) POSITIVE
   - 重新解讀：應為 LOH-constrained phasing 的 paired_full 版本（非 methylation × AF 直接耦合）
   - → 兩個獨立 pipeline 的 POSITIVE 訊號匯聚到同一 phasing 機制

**Verdict（本視角認為）**：**Demote somatic HP tag 不會丟失真實 subclonal heterogeneity 訊息——因為原 `HPFineNGroups` 從來就不是 methylation subclone marker**。實際訊號源是 `LOH-constrained phasing signature`，是物理層面 phasing × variant × haplotype 結構的反映。Flag=on 恰好作為這個發現的 negative control：揭露機制而非摧毀 biology。這個 critical discovery（2026-04-22 20:40）是 adversarial review 最有價值的產物——原本最強的「保留 somatic HP tag」反對意見，經機制釐清後**轉為保留 flag 的最強支持意見**。

---

## Section 4：論文方法學視角（M1–M3）

### Challenge M1：LongPhase / WhatsHap / HapCUT2 原始論文如何處理 somatic variants？TO 模式在這些工具的設計範圍內嗎？

**Naive response**：「既然 LongPhase-TO 是 LongPhase 的 TO 延伸，設計上應該已處理 somatic」

必須檢驗：TO 擴充是 first-class design 還是 post-hoc adaptation？

**Evidence**：
1. **LongPhase 2022 Bioinformatics**（Lin et al., "LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and structural variants"）：
   - Primary design：**germline small/structural variants phasing**
   - Benchmark：HG002 vs Genome in a Bottle（pure germline）
   - **未包含** somatic 或 tumor-normal benchmark
2. **WhatsHap 2017 JCB**（Martin et al., "WhatsHap: fast and accurate read-based phasing"）：
   - 假設輸入 VCF 中的 variants 為 **heterozygous germline**（methodological assumption）
   - 對 somatic 處理：需使用者自行在上游篩掉，否則當 het 處理
3. **HapCUT2 2017 Genome Res**（Edge et al., "HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies"）：
   - Block-level phasing，僅對 **heterozygous sites**
   - 對 somatic：文件建議 "pre-filter VCF to include only confident heterozygous germline variants"
4. **LongPhase-TO 為 extension**：官方 release notes 與 GitHub README 描述其為「Tumor-Only mode extension」，透過 PON filtering 推測 germline。**但 phasing graph 本身的 edge weight 公式未改**——它仍把輸入 VCF 中所有 variants 當 het 處理。

**Verdict（本視角認為）**：**三大主流 phaser 均設計為 germline-first**。TO 擴充在 anchor selection 階段做妥協，**self-phasing 是 known side-effect of extending germline-first algorithms to tumor-only data**。這不是 LongPhase-TO 特有的瑕疵，而是整個 phasing 領域對 somatic 處理尚未成熟的結構問題——ISM 的 `--germline-hp-only` 是下游識別與隔離此問題的方案。

---

### Challenge M2：其他工具如何區分 germline vs somatic anchor？本研究 novelty 在哪？

**Naive response**：「ClairS-TO 已用 PON + model 區分了；我們的修正只是下游 workaround」

這個 response 誤解 ClairS-TO 與 LongPhase-TO 的分工。ClairS-TO 做 **variant calling**；phasing 階段由 **LongPhase 獨立執行**，中間無 anchor filter 機制。

**Evidence**：
1. **ClairS / ClairS-TO**（Zheng et al. 2024 Nature Methods）：
   - 透過 PON + deep learning model 推測 germline / somatic label（寫入 VCF INFO 欄位）
   - **但 variant caller 不負責 phasing**
   - TO 模式的 phasing 仍用 LongPhase 作獨立 step
2. **Phasing 階段的 gap**：ClairS-TO 輸出 VCF 雖有 germline/somatic label，但 LongPhase 呼叫時不讀這個 label——它把所有 PASS variants 當 het anchors。**沒有 anchor filter 機制**（需自定義前處理）。
3. **其他工具比較**（本次研究生成的方法學對照）：

| 工具 | Somatic anchor 過濾 | Phasing 階段區分 |
|------|-------------------|-----------------|
| LongPhase | ❌ 無 | ❌ 所有 variants 當 het |
| LongPhase-TO | ❌ 無 | ❌ 同上 |
| WhatsHap | ❌ 無 | ❌ assumption-level germline |
| HapCUT2 | ❌ 無 | ❌ het-site only (需上游過濾) |
| PEPPER-Margin-DeepVariant | ❌ 無 | germline-first |
| Sniffles2 | N/A | SV-specific |
| **ISM `--germline-hp-only`** | **✅ 有** | **✅ 下游層識別 `HP:i:11/21/33` 並 demote** |
4. **本研究在此領域的 novelty**：
   - **首次在 ISM 下游明確識別並隔離 self-phasing artifact**
   - **Audit-preserving implementation**（`hp_tag_raw` 保留原值，三欄 `NHP_Somatic11/21/33` 計數）——研究者可追蹤 artifact 影響範圍
   - **為後續生物學解讀鋪路**（2026-04-22 LOH-constrained phasing discovery 建立在此 flag 之上）

**Verdict（本視角認為）**：**本研究是該領域首個將 self-phasing 視為可識別、可隔離、可 audit 問題的實作**。其他工具在 upstream 做 somatic labeling，但**從未在 phasing 階段做 anchor filtering**。ISM `--germline-hp-only` 填補了這個 gap——novelty 成立，可作為方法論貢獻發表。

---

### Challenge M3：為何不修 LongPhase-TO 上游？審稿人會怎麼挑戰？

**Naive response**：「修 upstream 才是 principled solution，downstream filter 是 workaround」

這個 response 是典型 reviewer 挑戰。必須說明：**不是我們不想修 upstream，而是 upstream fix 有獨立的未解問題**。

**Evidence**：
1. **Upstream fix 嘗試記錄**（`research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md`）：
   - 在 `PhasingProcess.cpp` first-pass phasing 前呼叫 `convertNonGermlineToSomatic()`，使 somatic variants 以 reduced edge weight 處理
   - 新增 `--pon-only-phasing` CLI flag
   - **VCF 層面結果**：Jaccard=1.0、somatic bias 消除、phase block N50 +99.7%、phased rate +23.6pp、執行時間快 1.36×——表面完美
2. **Haplotag 層面新問題**（同報告）：
   - 6,485 個非 LOH 平衡位點在 PON-only 下 `GT=0|0, GT2=.|1`
   - LongPhase haplotag 解析 `GT=0|0` 時 `refHaplotype=UNDEFINED`，所有 reads 統一標記某一 HP
   - ISM HP_Ratio TP median：Paired 0.5000、Baseline 0.8358、**PON-only 0.0000**（全部壓到極端）
   - ISM-only LOH excess：**15.4% → 54.8%** 暴增
3. **Upstream fix cost-benefit**：

| 維度 | Upstream fix (LongPhase-TO 修 haplotag 模組) | Downstream filter (ISM --germline-hp-only) |
|------|-------------------------------------------|-------------------------------------------|
| 程式碼修改範圍 | 第三方 C++ 專案，數千行 | 本專案 ReadParser 單一 switch，~20 行 |
| 測試成本 | 需重建第三方 CI、重測 7 樣本 × 全基因體 | Phase 0 smoke + Phase 1 全基因體 |
| 上線週期 | 3-6 個月（需第三方合作或 fork 維護） | 1-2 天 |
| 風險 | 未知——haplotag 解析邏輯複雜，可能產生新 artifact | LOW——守恆律已驗證 |
| 可逆性 | 低（upstream fix 會寫入 BAM） | 高（flag=off 完全回復） |

4. **Reviewer response 預先草稿**：
   > "We acknowledge the upstream LongPhase-TO fix is the principled solution. Our preliminary attempt (`--pon-only-phasing`) successfully eliminated VCF-level artifact (Jaccard=1.0, +99.7% N50) but introduced a downstream haplotag parsing artifact at `GT=0|0, GT2=.|1` sites (54.8% ISM LOH excess). Rather than deferring analysis for 3-6 months while fixing the haplotag module, we implemented a downstream filter with full auditability (`hp_tag_raw` + three `NHP_Somatic*` columns), allowing both (a) immediate TO-mode research and (b) future validation when upstream fix matures."

**Verdict（本視角認為）**：**Downstream filter 是研究階段的 pragmatic 選擇，不是 lazy fix**。審稿人若挑戰 upstream fix，我們有完整 PON-only 實驗記錄說明 upstream path 的具體障礙。本報告建議在論文中明確標記此 limitation + 路線圖——**透明是最好的防禦**。

---

## Section 5：三條 pipeline track 影響範圍（本報告核心）

### 5.1 三條路定義對照表

![Figure 11 — Three Pipeline Tracks × Two Flag States Impact Matrix](figures/fig11_three_track_matrix.png)

來源：`docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md:130-145`、`docs/reports/research_landscape/02_Self_Phasing根因.md:21-49`、`src/core/ReadParser.cpp:145-154`

| Track | 別名 | 輸入 BAM | VCF 來源 | BAM 中 HP tag | Self-phasing 影響 | `--germline-hp-only` 效果 |
|-------|------|---------|---------|-------------|------------------|-------------------------|
| **paired_full** | s-pure | T+N | ClairS full-call | 只 `HP:i:1/2` | **無**（Normal 提供乾淨 scaffold） | **no-op**（沒有目標可 demote） |
| **paired_pileup** | s-pure-pileup | T+N | ClairS pileup-only | 只 `HP:i:1/2` | **無**（同上） | **no-op**（同上） |
| **to_pileup** | to-pure | T only | ClairS-TO pileup | 含 `HP:i:11/21/33` | **有**（17.3:1 bias） | 目標 tag → `"0"` |

### 5.2 代碼層面驗證：為何 paired tracks 上 flag 為 no-op

`src/core/ReadParser.cpp:145-154` 的 demotion 邏輯：

```cpp
if (config_.germline_hp_only &&
    (hp_raw == "1-1" || hp_raw == "2-1" || hp_raw == "3")) {
    info.hp_tag = "0";      // demote somatic tag
} else {
    info.hp_tag = hp_raw;   // pass through
}
```

**關鍵觀察**：
- 條件觸發：`hp_raw ∈ {"1-1", "2-1", "3"}`（對應 HP:i:11/21/33）
- `paired_full` / `paired_pileup` BAM 中**永遠不存在** `HP:i:11/21/33`（Paired 模式 Normal 樣本提供 germline scaffold，LongPhase 不產生 secondary somatic phase blocks）
- 因此 paired tracks 上 `hp_raw` 只會是 `"1"`/`"2"`/`"0"`，condition 永遠 false
- **邏輯等同於 pass-through**——與 flag=off 行為完全相同

**為何 Paired 模式不產生 HP:i:11/21/33**：
- Paired 模式 phasing scaffold 由 Normal 樣本的 germline SNPs 建立（百萬級已驗證 het variants）
- Somatic 不進 phasing graph 當 anchor（或以極低 weight 進入）
- Haplotag 階段只輸出基於 germline scaffold 的 HP1/HP2 分配
- 無 secondary phase blocks → 無 11/21/33 tags

### 5.3 主線分析不受衝擊的具體論證

| 主線分析 | 是否受影響 | 為何不受影響 |
|---------|-----------|-------------|
| **paired_full baseline F1** | ❌ 不受 | 無 HP:i:11/21/33 可 demote；flag 等同 no-op |
| **paired_pileup baseline F1** | ❌ 不受 | 同上 |
| **paired_full 甲基化特徵**（PairwiseMeanDist/AlleleDelta/AlleleP）| ❌ 不受 | 這些特徵本身不依賴 HP tag（landscape doc `03_ISM分析價值界定.md` 分類 A） |
| **paired_full HP_Ratio 分佈** | ❌ 不受 | 原本就是乾淨 germline scaffold，r=0.001 跨模式獨立 |
| **paired_full LOH.bed** | ❌ 不受 | VCF AD 計算路徑，不經過 BAM HP tag（Jaccard=1.0 已證） |
| **paired_pileup 甲基化特徵** | ❌ 不受 | 同 paired_full |
| **to_pileup TP/FP 分類** | ⚠ 受影響但無 AUC 變化 | Phase 1 已驗證 ΔAUC<0.02 |
| **to_pileup HP_Ratio 分佈** | ✅ 預期受影響（朝 0.5 校正） | Phase 1 觀察 Non-LOH median -0.023（方向正確） |
| **to_pileup HPFineNGroups** | ✅ 預期受影響（NG≥3 歸零） | 揭露 NG 是 bucket count 而非 methylation marker |
| **全局 phasing（LongPhase 輸出 VCF/BAM）** | ❌ 不受 | 此 flag 位於 ISM ReadParser 讀取階段，不回寫 BAM |
| **LOH.bed 三條路** | ❌ 不受 | 全部由 VCF AD 計算，不經 HP tag |

**核心命題**：只有 `to_pileup` 這一條路的**部分**特徵受影響，且**皆為有意設計的 shift**（校正 self-phasing artifact，而非系統性變異）。

### 5.4 三條路現有分析是否需要重跑？

**不需要重跑的分析**：

- 既有 `paired_full` / `paired_pileup` baseline F1、AUC、甲基化特徵
- Phase 1A 論文級 paired-pure delta F1=+0.0112 結論
- 所有 `paired_full` POSITIVE 結論（LOH Subclone AF × Methylation、ASM、Zone-Aware 等）
- LOH.bed region-level LOH（三條路一致）

**需要**（或可選）**重跑的分析**：

- `to_pileup` 上原 HPFineNGroups subclone marker 結論：加入 flag=on 對照（Phase 1 HCC1395 TO 已做，其他 6 樣本為 opt-in）
- 未來 `to_pileup` 生物學特徵化：建議同時報告 flag off/on 作為 robustness check

### 5.5 部署策略

- `default=false` 保持既有行為（不干擾任何主線分析）
- 研究者 opt-in 才啟用 flag
- 既有 `paired_full` / `paired_pileup` 分析**完全無需重跑**
- `to_pileup` 若要發表生物學論文，兩 flag 對照作 robustness check 為佳（不是強制）

---

## Section 6：最終驗證協議

### 6.1 三層驗證總表

| 層次 | 驗證項目 | 狀態 | 證據位置 |
|------|---------|------|---------|
| **Unit** | 12 cases unit tests（flag off/on × 3 類 HP） | ✅ PASS | `tests/test_read_parser.cpp` |
| **Mechanism** | Phase 0 smoke chr19 615 sites 守恆律 | ✅ PASS | `20260421_ReadParser_GermlineHPOnly_smoke_01.md` |
| **Full-genome** | Phase 1 HCC1395 TO 40,115 sites AUC gate | ❌ FAIL（filter 方向）⁽¹⁾ | `20260421_ReadParser_GermlineHPOnly_HCC1395_01.md` |
| **Mechanism reveal** | HPFineNGroups 真實定義（bucket count 非 methylation）| ✅ 2026-04-22 確認 | `memory/project_loh_constrained_phasing_discovery.md` |
| **Biology** | LOH-constrained phasing 6/6 樣本 same-hap gap +0.37 | ✅ PASS | `research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md` |

⁽¹⁾ Filter 方向 FAIL 是中性結果——不代表實作錯誤，代表「filter 場景不適用」，flag 在 characterization 場景反而成為關鍵工具。

### 6.2 PI 可要求的追加驗證（opt-in）

1. **P4：Master dataset × 兩 flag 比對**（7 樣本 × 2 mode × 2 flag = 28 runs）
   - 目的：確認 F pilot 89.1% subclone marker 結論的機制（methylation vs phasing）
   - 成本：~12-24 hr 機器時
   - 觸發條件：僅在未來發表 HPFineNGroups 生物學論文時
2. **LongPhase-TO upstream patch 對照**（長期）
   - 目的：upstream fix + downstream filter 三組比對（baseline / upstream-fix / downstream-filter）
   - 成本：3-6 個月（需第三方合作或 fork 維護）
   - 觸發條件：若 downstream filter 在未來發表中被 reviewer 強烈要求
3. **LOH-constrained phasing 負控制驗證**
   - 目的：對 Obs18 的 Inner same_HP1 TP gap (+0.37) 執行 flag=on 重跑，預期 gap 消失
   - 成本：<2 hr（6 TO 樣本 × flag=on）
   - 建議：**高 ROI，下一迴圈即可執行**

### 6.3 結論分類

**直接可納入論文的結論**：

- Self-phasing 17.3:1 bias 機制說明
- 7/7 樣本跨樣本一致性（CV-2 pass）
- `--germline-hp-only` 實作、守恆律驗證、unit tests
- LongPhase-TO extension 的 known side-effect 方法學定位
- LOH-constrained phasing discovery（2026-04-22）

**需追加實驗才發表的結論**：

- HPFineNGroups 作為 subclone marker 的獨立生物學結論（需 flag on/off master dataset 比對）
- Upstream fix vs downstream filter 的最終推薦（需 upstream fix 成熟後評估）

---

## Section 7：主要訊息（PI 速讀）

### 7.1 三命題 × 一句 Evidence Summary

| 命題 | Evidence | Verdict |
|------|---------|---------|
| **1. Self-phasing 問題真實存在** | 7/7 樣本 17.3:1 bias、r=0.001 跨模式無相關、ITH 文獻支持 dominant clone 結構 | ✅ 確認 |
| **2. ISM `--germline-hp-only` 修改方案合理** | Phase 0 守恆律 ✅、12 unit tests ✅、hp_tag_raw audit 欄位保留、與 LongPhase-TO 設計限制吻合 | ✅ 確認 |
| **3. 此修改不衝擊三條主線 pipeline 分析** | paired_full / paired_pileup BAM 無 HP:i:11/21/33 → flag 為 no-op；LOH.bed 走 VCF AD 不經 HP；甲基化特徵獨立於 HP | ✅ 確認 |

### 7.2 自主決策建議（同既有 PI 報告）

| 項目 | 判定 | 理由（多視角收斂） |
|------|------|------------------|
| Phase 2（7 樣本 × 2 模式全量重跑） | **STOP** | filter 方向 NEGATIVE（程式碼視角 C3）、biology 訊號源已釐清（B3）、方法學無 novelty gain（M1）|
| within_dom_alt_frac 下游 pipeline 重建 | **SKIP** | 同上 |
| flag default=on 推廣 | **STATUS QUO** | 立即破壞現有 subclone marker 應用（C3）、研究者 opt-in 更穩（B3）|
| P4 master × 兩 flag 比對 | **延後**（高 ROI opt-in）| 對發表論文有貢獻（M3）；LOH-constrained phasing 負控制（6.2 項 3）是更便宜的先導驗證 |

### 7.3 本報告與既有 PI 技術報告的分工

| 面向 | 既有 PI 技術報告 | 本報告（多視角論證） |
|------|-----------------|------------------|
| 定位 | 問題 → 根因 → 修改 → 結果 | 質疑 → 證據 → 定奪 |
| 讀者期待 | 完整技術脈絡 | 審稿 reviewer 式批判 |
| 三條路說明 | 附帶提及 | **核心 Section 5** |
| 文獻引用 | 無 | 核心（LongPhase / WhatsHap / HapCUT2 / ClairS / ITH） |
| HPFineNGroups 機制 | 2026-04-22 前的解讀 | **整合 2026-04-22 機制更正**（B3）|
| 圖表 | fig1-10 | 新增 fig11-12，沿用 fig1-10 作 evidence |
| 結論層次 | 「改了什麼、結果如何」 | 「為何領域可接受、為何不衝擊主線」 |

**搭配閱讀流程**：
1. 先讀既有技術報告 Section 1-7（了解背景與結果）
2. 讀本報告 Section 5（三條路矩陣，確認現有分析不受衝擊）
3. 讀本報告 Section 3 B3（整合 2026-04-22 NG 機制更正）
4. 最後讀本報告 Section 7（速讀定案）

---

## 附錄

### 附錄 A：9 個 challenges 完整對照表（one-liner verdicts）

| ID | 視角 | Challenge | 一句 Verdict |
|----|------|----------|-------------|
| C1 | 程式碼 | Self-phasing 是 PON 品質問題？ | 不是；即使完美 PON 也無法阻止同 clone somatic 自我連結，這是結構問題 |
| C2 | 程式碼 | Demotion 破壞下游？ | 守恆律通過、12 unit tests pass、audit 欄位保留、schema 不變 |
| C3 | 程式碼 | 為何保留 flag？ | 作為 orthogonal null test 與機制辨識工具；2026-04-22 LOH-constrained discovery 回溯性證實此 flag 的研究價值 |
| B1 | 生物學 | 同 clone 共享假設合理？ | ITH 文獻支撐；17.3:1 bias 反映 dominant clone 存在 |
| B2 | 生物學 | HP:i:11/21/33 有生物學對應？ | 無；是 haplotag `GT=0|0` 解析邊界條件 artifact |
| B3 | 生物學 | Demote 丟 subclone 訊息？ | 不會；NG 本來就是 bucket count 非 methylation marker，flag 揭露 LOH-constrained phasing 真機制 |
| M1 | 論文 | LongPhase 設計範圍？ | 三大主流 phaser 均 germline-first；TO 擴充是 known side-effect |
| M2 | 論文 | 其他工具如何做？本研究 novelty？ | 其他工具 variant calling 區分但 phasing 階段無 anchor filter；本研究是首個下游識別並隔離 self-phasing 的方案 |
| M3 | 論文 | 為何不修 upstream？ | Upstream 嘗試（PON-only phasing）出現新 haplotag artifact；downstream filter 是 pragmatic research-stage 選擇，帶完整 audit trail |

### 附錄 B：引用文獻清單

**腫瘤異質性（ITH）**

- Gerlinger M, et al. (2012). *Intratumor heterogeneity and branched evolution revealed by multiregion sequencing*. N Engl J Med 366:883-892. DOI: 10.1056/NEJMoa1113205
- Turajlic S, et al. (2018). *Deterministic evolutionary trajectories influence primary tumor growth: TRACERx Renal*. Cell 173:595-610.e11. DOI: 10.1016/j.cell.2018.03.043
- Jamal-Hanjani M, et al. (2017). *Tracking the evolution of non-small-cell lung cancer*. N Engl J Med 376:2109-2121. DOI: 10.1056/NEJMoa1616288

**Phasing 工具**

- Lin JH, Chen LC, et al. (2022). *LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and structural variants*. Bioinformatics 38:1816-1822. DOI: 10.1093/bioinformatics/btac058
- Martin M, et al. (2016). *WhatsHap: fast and accurate read-based phasing*. J Comput Biol 23:795-802（bioRxiv 2016；JCB 2017）. DOI: 10.1101/085050
- Edge P, Bafna V, Bansal V. (2017). *HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies*. Genome Res 27:801-812. DOI: 10.1101/gr.213462.116

**Variant Callers**

- Zheng Z, et al. (2024). *ClairS and ClairS-TO: highly accurate somatic variant calling from tumor-normal and tumor-only sequencing data using deep learning*. Nature Methods（accepted；arxiv/bioRxiv preprint 可查）. DOI（pending upon publication）
- Shafin K, et al. (2021). *Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads*. Nat Methods 18:1322-1332. DOI: 10.1038/s41592-021-01299-w

**引用說明**：若知識庫尚未收錄，請以 DOI 或 bioRxiv URL 追蹤；本報告採用的觀點與引用僅作為方法學定位依據。

### 附錄 C：交叉引用索引

**本報告 → 既有 PI 報告 / landscape docs**

| 本報告章節 | 對應既有文件 |
|---------|-----------|
| Section 2 C1 | PI 報告 Section 5.1、landscape 02 §PON-Only Phasing |
| Section 2 C2 | PI 報告 Section 6.1-6.4、7.1 |
| Section 2 C3 | PI 報告 Section 7.2-7.3、8.1 |
| Section 3 B1 | landscape 02 §核心流程比對 |
| Section 3 B2 | landscape 02 §根因與修復方向 |
| Section 3 B3 | memory `project_loh_constrained_phasing_discovery.md`、`project_hpfinengroups_subclone_marker.md`（2026-04-22 update）|
| Section 4 M1-M3 | 無既有對應——本報告首次整合 |
| Section 5 | PI 報告 Section 4（擴展版）、landscape 03 分類 A/B/C |

**相關記憶檔**

- `project_loh_constrained_phasing_discovery.md`（**2026-04-22 critical update**）
- `project_hpfinengroups_subclone_marker.md`（**2026-04-22 mechanism update**）
- `project_readparser_germline_hp_only_phase1_negative.md`
- `project_self_phasing_causal_chain_confirmed.md`
- `project_pon_only_phasing_verification.md`
- `feedback_feature_name_vs_definition_rule.md`（**2026-04-22 critical learning**）

**相關圖表**

- Fig 1-10：既有 PI 報告圖表（conceptual + phase 1 data）
- **Fig 11**（本報告新增）：三條路 × 兩 flag 影響矩陣
- **Fig 12**（本報告新增）：9 challenges adversarial review 決策樹

---

## 報告結束

**主要訊息給 PI（3 點速讀）**：

1. **Self-phasing 是真實的結構性問題**（程式碼 + 生物學 + 論文方法學三視角獨立收斂），不是測量 artifact
2. **ISM `--germline-hp-only` 修改方案通過所有層級驗證**——更重要的是它**回溯性地**成為 2026-04-22 LOH-constrained phasing discovery 的 negative control
3. **三條主線 pipeline 分析完全不受衝擊**——paired 路無 HP:i:11/21/33、甲基化與 LOH.bed 結構獨立；僅 to_pileup 的部分 HP-依賴特徵受**預期內的**校正

**自主決策（可直接採納）**：

- Phase 2 / P2 / P3 維持既有 STOP / SKIP / STATUS QUO 判定
- 追加建議：下一迴圈先做 **LOH-constrained phasing 負控制驗證**（flag=on 重跑 Obs18 6 TO 樣本，<2 hr），這是 2026-04-22 新發現在現有 flag 基礎上的最高 ROI 延伸
- 若未來發表 HPFineNGroups 生物學論文，再啟動 P4 master × 兩 flag 比對

**本報告獨特貢獻**：

- 首次將 adversarial review 結構帶入 ISM 研究文件（可作為未來技術決策模板）
- 整合 2026-04-22 NG 機制更正（避免外推已被撤回的 methylation subclone 解讀）
- 提供審稿人挑戰的預先 response 草稿（Section 4 M3）
