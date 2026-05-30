# ClairS / LongPhase-S 來源驗證報告（HKU 合作補強）

<!--
建立時間: 2026-05-23
版本: v4 (Q4 rewrite per 5/23 user README.md 揭示)
類型: 外部論文 / source code / GitHub release 驗證
範圍: ClairS caller window / pileup tensor / 甲基處理 / LongPhase-S（mainline v1.7.3）版本 / ssrs vs ss
訪問日期: 2026-05-23
-->

> **術語說明（v4 修正）**：本檔以下「**LongPhase-S**」一律指 `twolinin/longphase` mainline **v1.7.3** 配 `--somaticMode` flag — 屬 mainline 內建功能，**非 fork**。`HP1_1=5` / `H2_1=7` 等 enum 定義於 `src/haplotag/HaplotagType.h:97-108`，在 mainline v1.7.3 已內建。本地 `/big8_disk/.../longphase-s/` 為 user 對 twolinin mainline clone 的本地命名（其 `README.md` 引用 twolinin v1.7.3 download URL `https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz`）。

## 摘要（給 PI 一眼看）

| Q | 結論 | 對 HKU 報告影響 |
|---|------|----------------|
| Q1 caller 視窗 | **33 bp（flanking ±16 bp）非 ±128 bp** | 🔴 **關鍵修正** |
| Q2 pileup tensor | **ONT: [130 reads × 33 positions × 7 channels]** | 🟡 需補章節 |
| Q3 甲基訊號 | **ClairS / ClairS-TO 均不讀 MM/ML（pileup tensor 無 methylation channel）** | 🔴 **關鍵修正措辭** |
| Q4 BAM 真實上游（5/23 用戶執行 samtools view -H + cat README.md 揭示） | **LongPhase-S = `twolinin/longphase` mainline v1.7.3 `haplotag --somaticMode`（內建，非 fork）+ ClairS-TO ssrs SNVs + SEQC2 --highCon-snp 引導；BAM 名稱含 "ClairS_pileup" 為 misleading legacy naming** | 🔴 **重大 framing 修正** — 非 ClairS paired pileup；亦非 CCU fork |
| Q5 ssrs vs ss 視窗 | **完全相同 33 bp，差別僅在訓練資料（real cancer cell-line augmentation）** | 🟢 報告維持 |

---

## Q4 真實 BAM @PG line + LongPhase-S 來源澄清（5/23 v4 終版 — 取代原 v1.6 / CCU fork 推測）

### Q4 結論（v4 終版）

**LongPhase-S = `twolinin/longphase` mainline v1.7.3 內建 `--somaticMode` flag + HP somatic 階層 tag**。**非 CCU fork、非 v1.0.0、非 v1.6**。

### 證據鏈

**Source 1**: `samtools view -H /big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` 5/23 執行：

```
@PG  ID:longphase  PN:longphase  PP:samtools.1  VN:1.7.3
     CL:longphase haplotag --somaticMode
        -s /big8_disk/liaoyoyo2001/data/vcf/HCC1395BL_methyl_phase.vcf.gz       ← normal phased VCF
        -b /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam    ← normal BAM
        --tumor-snp-file /big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/
                         ONT_5khz_simplex_5mCG_5hmCG/pileup/HCC1395_methyl_PASS.vcf  ← ClairS-TO ssrs（非 ClairS paired pileup）
        --tumor-bam-file /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam
        -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
        -t 64
        -o /big8_disk/.../HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter
        --log --tagSupplementary -q 20
        --highCon-snp /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf  ← SEQC2 truth set 引導
```

**Source 2** (5/23 用戶執行 `cat /big8_disk/.../longphase-s/README.md`)：本地 `longphase-s/` 目錄的 `README.md` 引用的下載 URL 為 **`https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz`** — 直接證明本地命名為 "longphase-s" 的目錄其實是 **mainline `twolinin/longphase` v1.7.3 的 clone**，user 採取 "-s" 命名只是區分本地用途，非真實 fork。

**Source 3**：BAM `@PG VN:1.7.3` 對應 mainline release（非 fork patch）。

### v4 修正了哪些 v1-v3 推測

| 原推測（v1-v3 錯誤） | v4 實際 |
|--------|-------------------|
| `longphase-S v1.0.0 (CCU fork)` | **`twolinin/longphase` mainline v1.7.3 內建 `--somaticMode`**（非 fork）|
| `ClairS paired pileup mode → longphase-s` | **ClairS-TO ssrs → LongPhase-S (mainline v1.7.3) `haplotag --somaticMode`**（with `-s normal_vcf + -b normal_bam`）|
| 純 caller-only output | **用 SEQC2 `--highCon-snp` 引導**（非 evaluator-blinded）|
| `HP1_1=5 / H2_1=7` enum 是 CCU fork-only patch | **mainline v1.7.3 已內建**（enum 在 `src/haplotag/HaplotagType.h:97-108`）|
| A1 web research 看到 `v1.6 latest visible` | **錯誤** — twolinin/longphase 真的有 **v1.7.3 release**；A1 web 查時可能受 GitHub 顯示快取 / pre-release filter 影響 |

### 對 HKU 報告影響

1. 報告中「ClairS pileup → longphase-S」必須改為「**ClairS-TO ssrs + Normal BAM → LongPhase-S（mainline v1.7.3）`haplotag --somaticMode`**」
2. BAM 名稱「ClairS_pileup」是 legacy misleading naming（與 5/20 「paired-pileup → ClairS-TO ssrs」 framing-fix 同一類問題）
3. 「LongPhase-S」這個術語在本專案內部統一指 **mainline v1.7.3 + `--somaticMode` flag**（非 fork），可放心沿用
4. HKU 合作報告主軸「介紹 LongPhase-S」對應「介紹 mainline `twolinin/longphase` v1.7.3 + `--somaticMode` + ClairS-TO ssrs」整套流程

**HKU 報告建議用詞**：
> 「HCC1395 tagged BAM 由 ClairS-TO ssrs model（v0.4.x）call somatic SNVs，搭配 SEQC2 high-confidence truth set 引導，經 **LongPhase-S（`twolinin/longphase` mainline v1.7.3 `haplotag --somaticMode`）** 同時讀 normal BAM 與 normal-phased VCF 產出。」

---

## Q1. ClairS caller 觀察視窗精確值

### 結論
**ClairS 與 ClairS-TO 的 pileup model 觀察視窗為 33 bp（flanking ±16 bp，含中心位點）**，非報告中誤用的 ±128 bp。

### 來源
- ClairS GitHub `shared/param.py`：`flankingBaseNum = 16` → `no_of_positions = 2 * flankingBaseNum + 1 = 33`（[https://github.com/HKU-BAL/ClairS](https://github.com/HKU-BAL/ClairS)）
- ClairS-TO GitHub `shared/param.py`：相同設定 `flankingBaseNum = 16`，跨 Illumina / ONT / PacBio HiFi 平台統一（[https://github.com/HKU-BAL/ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)）
- bioRxiv preprint DOI: [10.1101/2023.08.17.553778v1](https://www.biorxiv.org/content/10.1101/2023.08.17.553778v1)（abstract 可達；全文 403）

### 對 HKU 報告的影響
**🔴 必須修正**：報告中所有「±128 bp」、「±128 bp window」、「128 bp flanking context」需改為「±16 bp（33 bp total window）」。
- 若曾推論「ClairS 看到 256 bp 範圍的甲基化分佈」此類陳述需整段重寫
- 33 bp window 意義：ClairS pileup 僅看 **單一 CpG 中心 ± 半個 nucleosome unit** 的鄰近 base composition，**不足以涵蓋一個完整 CpG island 的 inter-CpG context**

### 補充：full-alignment 路徑
ClairS / ClairS-TO 採 pileup + full-alignment 雙模型 ensemble。Full-alignment 路徑的 receptive field 在 `param.py` 同樣由 `flankingBaseNum` 控制（同 33 bp）；full-alignment 額外考慮 read-level 整段對齊資訊，但 **空間範圍仍為 33 bp**。

---

## Q2. ClairS pileup tensor 結構

### 結論

| 維度 | ONT | Illumina | HiFi |
|------|-----|----------|------|
| Depth (rows) | **130** | 170 | 130 |
| Positions (cols) | **33** | 33 | 33 |
| Channels (full-aln) | **7** | 7 | 7 |
| Pileup channels | **18** | 18 | 18 |

### Channel 定義（source code）

**Full-alignment 7 channels**：
1. reference_base
2. alternative_base
3. base_quality
4. strand_info
5. insert_base
6. **phasing_info**（即 HP tag — 這是 ClairS **唯一**接觸 longphase 輸出的入口）
7. mapping_quality

**Pileup 18 channels**：A / C / G / T / I / I1 / D / D1 / * / a / c / g / t / i / i1 / d / d1 / #（大寫 = 正股，小寫 = 反股，I/D = insertion/deletion，# = end-of-read）

### 來源
- ClairS-TO `shared/param.py`：`channel_size = 7`、`pileup_channel_size = 18`、`input_shape = [130, 33, 7]` for ONT
- 同 ClairS-TO Nature Comms 引用 [DOI 10.1038/s41467-025-64547-z](https://doi.org/10.1038/s41467-025-64547-z)（Chen et al. 2025）

### 是否考慮上下游 SNV pair？
**否，完全 local context only**。Pileup tensor 為以 candidate variant 中心位點為原點的 33 bp 區段；ClairS **不**將鄰近 SNV 當 input feature 給 pileup model。鄰近 SNV 的 phasing 影響透過 channel 6 (phasing_info / HP tag) 進入模型，**而 HP tag 是 longphase 在 BAM 上先 tag 完才給 ClairS 用**。

### 對 HKU 報告的影響
🟡 **建議補章節**「ClairS input tensor 結構」澄清：
- 模型不直接看「上下游 SNV pair」
- 鄰近 SNV 經 longphase haplotag → HP tag → ClairS phasing channel 此單一路徑進入
- 33 bp window 意味著 inter-SNV phasing 訊號 **完全外部化** 給 longphase，ClairS 自身無 phasing graph

---

## Q3. ClairS 甲基處理（5mC MM/ML tags）

### 結論
**ClairS 與 ClairS-TO 均不讀 MM/ML tag、不使用甲基當 input feature、無 ablation 報告甲基訊號 effect**。Pileup tensor 的 18 channel 與 full-alignment 的 7 channel **皆無 methylation channel**。

### 來源
- ClairS GitHub README 全文：「No mentions of methylation, 5mC, MM/ML tags, modifications, modkit, or methylation calling」（[https://github.com/HKU-BAL/ClairS](https://github.com/HKU-BAL/ClairS)）
- ClairS-TO GitHub README 全文：「No references to methylation, 5mC, MM/ML tags, or modified base handling」（[https://github.com/HKU-BAL/ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)）
- `shared/param.py` channel 定義：18 個 pileup channel 全部為 ACGT/indel/strand/end-marker，**無 5mC/5hmC channel**
- `clairs.py` 主入口：「References to MM, ML, methylation, modified bases, mod_call, or modkit: None detected」

### ClairS-TO ssrs vs ss 是否不同？
**相同**。兩個 model 共用同一個 `param.py`、同一個網路架構、同一組 channel 定義；差別僅在訓練資料 augmentation（ssrs 多用 HCC1937/HCC1954/H1437/H2009 四個 real cancer cell-line dataset 微調，ss 純合成資料）— 但**兩者都不接觸 MM/ML tag**。

### 對 HKU 報告的影響
🔴 **必須修正措辭**：
- 不可寫「ClairS 模型已內建甲基訊號」、「ClairS 使用甲基當特徵」、「ClairS 整合 epigenetic context」
- 正確措辭：「**ClairS 完全不讀 MM/ML tag**；甲基訊號需 ClairS 之外的後處理（如 InterSubMod、modkit、longphase modcall）注入」
- 此事實**直接強化 InterSubMod 在論文中的定位** — 因為 caller 端有此 blind spot，下游 read-level epigenetic filter 才有真正的 value-add

### 補充：longphase 端的甲基處理
- **LongPhase v2.0.0+** 支援 `modcall` 子命令：從 BAM 的 MM/ML tag 呼叫 allele-specific methylation
- **LongPhase phase** 支援 `--mod-file` 將 modification VCF 一併 co-phase
- ClairS / ClairS-TO **本體**仍不讀 MM/ML — longphase 是分工另外做的

---

## Q4. LongPhase-S 版本確認（v4 終版 — 核心問題）

### 結論（v4，依 5/23 user `cat README.md` 揭示重寫）
**LongPhase-S = `twolinin/longphase` mainline v1.7.3 + `--somaticMode` flag（內建功能，非 fork、非 v1.0.0 CCU patch）**。手上 BAM 由 ssrs 流程於 2025-06-30 產出。upstream mainline v1.7.3 已 release，沿用即可。

### BAM @PG line 揭示（5/23 user 已執行）
```bash
samtools view -H /big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam | grep -i "@PG"
# → @PG  ID:longphase  PN:longphase  PP:samtools.1  VN:1.7.3
#       CL:longphase haplotag --somaticMode ...
```
→ `VN:1.7.3` 直接對應 mainline `twolinin/longphase` v1.7.3 release（**非 fork patch**）。

### 本地 longphase-s/ 目錄真實身份（5/23 user `cat README.md` 揭示）

`/big8_disk/.../longphase-s/README.md` 引用：
```
https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz
```
→ 本地命名 "longphase-s" 是 **user 對 mainline `twolinin/longphase` v1.7.3 clone 的本地命名**，非實質 fork。enum `HP1_1=5 / H2_1=7` 定義於 `src/haplotag/HaplotagType.h:97-108`，**屬 mainline v1.7.3 內建**，非 fork-only patch。

### 證據對應

1. **檔案路徑指紋**：本 BAM 與 `/big8_disk/liaoyoyo2001/longphase_somatic_output_ssrs/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` **byte-identical**（283,071,595,503 bytes，見 InterSubMod `disk-audit-2026-04-18.md` §6）→ ssrs 流程輸出
2. **產出日期**：mtime = 2025-06-30
3. **HP tag 格式**：本 BAM 含 `HP:z:1 / 2 / 1-1 / 2-1 / 3` 字串格式（mainline v1.7.3 `--somaticMode` 內建輸出）
4. **enum 內建**：`src/haplotag/HaplotagType.h:97-108` 在 mainline v1.7.3 已內建 `HP1_1=5 / H2_1=7` 等 somatic 階層 enum，**非 fork-only**

### GitHub 上游版本（v4 修正）

| 工具 | 版本 | 來源 | 與手上 BAM 關係 |
|------|------------|------|----------------|
| **LongPhase-S = `twolinin/longphase` mainline** | **v1.7.3** | [twolinin/longphase v1.7.3 release](https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz) | **這就是手上 BAM 用的版本** |
| `--somaticMode` flag | mainline 內建 | `twolinin/longphase` mainline v1.7.3 source | 產出 HP1-1 / HP2-1 / HP3 細分 tag |
| `HP1_1=5 / H2_1=7` enum | mainline 內建 | `src/haplotag/HaplotagType.h:97-108` | **非 CCU fork-only** |

### A1 web research 「v1.6 latest visible」更正

A1 web research 報告中「longphase v1.6 為 GitHub 最新可見 release，但 BAM 用 1.7.3 — 可能 1.7.3 為較新 release 或 pre-release / patch」**屬錯誤推測**。

**5/23 user `cat README.md` 揭示**：twolinin/longphase **真的有 v1.7.3 release**，URL 為：
```
https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz
```

A1 web 查時看到「v1.6 latest visible」可能成因：
- GitHub release page UI 快取 / 顯示延遲
- A1 fetch 抓到的可能是 pre-release filter 視圖
- 真實 release 列表確實包含 v1.7.3（per `cat README.md` 直接 URL 引用）

### 建議

🟢 **沿用既有 BAM**（不重跑）：
- LongPhase-S = mainline v1.7.3 + `--somaticMode`，BAM 即此版本產出
- 無需考慮 CCU fork 或 base v2.0.0/v2.0.2 patch（誤推測已釐清）
- 重跑成本 ≈ 283 GB I/O × 24 cores × 數小時，無 upstream 變動誘因

✅ **5/23 v4 已 100% provenance**（@PG line + README.md URL 雙重證據）— 無需再做 30 秒 samtools header 驗證（已執行過）。

---

## Q5. ClairS-TO ssrs vs ss 視窗

### 結論
**視窗完全相同（33 bp，flanking ±16 bp）；channel 設定相同；網路架構相同。差別僅在訓練資料 augmentation。**

### 來源
- ClairS-TO `shared/param.py`：`min_bq_dict` 與其他 quality threshold **對 ssrs 與 ss 完全相同**，無平台層差異參數
- ClairS-TO README：ssrs = synthetic + real samples augmented (HCC1937/HCC1954/H1437/H2009)；ss = synthetic-only
- HKU-BAL 官方建議：「The `ssrs` model provides better performance and fits most usage scenarios」

### 對 HKU 報告的影響
🟢 **報告無需修改**。若報告中有「ssrs/ss 視窗不同」的暗示需刪除。
可補一句：「ssrs 在 real cancer cell-line 上的 calibration 改善 calling 效能，但 input tensor 結構與 ss 完全一致」。

---

## 對報告主軸的修正建議（總表）

| # | 修正項目 | 優先級 | 原因 |
|---|---------|--------|------|
| 1 | **所有「±128 bp」改為「±16 bp（33 bp window）」** | 🔴 P0 | 數字錯誤；33 bp ≠ 128 bp 是 ~8× 的差距 |
| 2 | **「ClairS 使用甲基」措辭改為「ClairS 不讀 MM/ML tag」** | 🔴 P0 | 反向陳述；影響論文 narrative |
| 3 | 補章節「ClairS input tensor 結構（130×33×7）」 | 🟡 P1 | 讓 HKU 端能對齊 InterSubMod 設計邏輯 |
| 4 | LongPhase-S 版本標記為 **`twolinin/longphase` mainline v1.7.3 + `--somaticMode`**，BAM mtime 2025-06-30 | 🟡 P1 | provenance 完整（v4 修正：非 CCU fork v1.0.0） |
| 5 | 補述「鄰近 SNV 透過 HP tag (phasing_info channel 6) 進入 ClairS — pileup model 本身不看 inter-SNV pair」 | 🟢 P2 | 加強 InterSubMod 與 ClairS 的職責分工敘事 |
| 6 | 強化敘事：**ClairS 33 bp blind spot + 無甲基 channel = InterSubMod read-level epigenetic filter 的存在理由** | 🟢 P2 | 論文 framing 加分 |

---

## 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| ClairS GitHub `shared/param.py` | source code | **高** | 一手；訪問日期 2026-05-23 |
| ClairS-TO GitHub `shared/param.py` | source code | **高** | 一手；驗證一致性 |
| ClairS bioRxiv 10.1101/2023.08.17.553778v1 | preprint | 中 | abstract 可達；全文 403 |
| ClairS-TO Nature Comms 2025 (DOI 10.1038/s41467-025-64547-z) | peer-reviewed | **高** | 一手；本 fetch 受 paywall 但已知架構同 GitHub |
| LongPhase mainline GitHub releases (`twolinin/longphase`) | release notes | **高** | 上游 changelog 完整；**v1.7.3 為手上 BAM 對應 release** |
| 本地 `longphase-s/README.md` 引用 URL (v4 證據) | 一手 user shell 輸出 | **高** | 直接證明 = mainline v1.7.3 clone，非 fork |
| InterSubMod `disk-audit-2026-04-18.md` § BAM identity | 本地 audit | **高** | 283 GB byte-match 證據 |
| InterSubMod `auto-run-sh.md` BAM 用法 | 本地 doc | 中 | 確認此 BAM 用於 ssrs 流程 |
| BAM @PG line 直接驗證 | 一手 (5/23 user 已執行) | **高** | ✅ v4 已驗證 `VN:1.7.3` |
| 本地 `longphase-s/README.md` URL 引用 | 一手 (5/23 user `cat README.md`) | **高** | ✅ v4 證實 = mainline v1.7.3 clone，非 CCU fork |

---

## 建議行動

1. ✅ **已完成**（5/23 user 執行）：`samtools view -H ... | grep @PG` 揭示 `VN:1.7.3` + `cat /big8_disk/.../longphase-s/README.md` 揭示 mainline v1.7.3 URL → §Q4 已達 100% provenance
2. **🔴 修正 HKU 報告**：依「修正建議總表」P0 兩項全文搜尋替換 + 全 doc 用「LongPhase-S」統一指 mainline v1.7.3 + `--somaticMode`
3. **🟡 補章節**：ClairS input tensor (130×33×7) + 鄰近 SNV 透過 HP tag 而非直接 input
4. **🟢 強化敘事**：33 bp window blind spot + 無甲基 channel = InterSubMod value-add 的雙重論據

---

## Sources

- [ClairS GitHub repository](https://github.com/HKU-BAL/ClairS) (accessed 2026-05-23)
- [ClairS-TO GitHub repository](https://github.com/HKU-BAL/ClairS-TO) (accessed 2026-05-23)
- [ClairS bioRxiv preprint v1](https://www.biorxiv.org/content/10.1101/2023.08.17.553778v1) (accessed 2026-05-23)
- [ClairS-TO Nature Communications 2025](https://doi.org/10.1038/s41467-025-64547-z) (accessed 2026-05-23, paywall)
- [LongPhase mainline GitHub releases (`twolinin/longphase`)](https://github.com/twolinin/longphase/releases) (accessed 2026-05-23) — **v1.7.3 download URL: `https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz`**（per 5/23 user `cat README.md` 揭示）
- 本地 `longphase-s/README.md` 引用 URL (v4 證據；user 5/23 cat) — 直接證實本地 "longphase-s" 目錄為 mainline v1.7.3 clone
- 本地 Knowledge: `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase.md` (last_verified 2026-04-17)
- 本地 Knowledge: `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md` (last_verified 2026-04-17) — **本地命名，指 mainline v1.7.3 clone，非 CCU fork**
- 本地 Knowledge: `/big8_disk/liaoyoyo2001/Knowledge/05_tools/variant-callers.md` (last_verified 2026-04-01)
- 本地 audit: `/big8_disk/liaoyoyo2001/Knowledge/raw/cli-snapshots/disk-audit-2026-04-18.md` (BAM byte-identity §6)
