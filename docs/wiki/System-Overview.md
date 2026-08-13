# System Overview 系統全景

[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> **這套系統由哪些部件組成，各自做什麼？**
> 本頁回答四個問題：**目標是什麼**、**有哪些部件**、**每個部件吃什麼吐什麼**、**現在哪些真的能跑**。
> 每個部件的細節在各自分冊頁，本頁只給地圖與狀態。
>
> 受眾：指導教授 + 實驗室同學（第一次接觸這套系統的人）

---

## 一分鐘版本

這套系統研究的問題是：**一顆腫瘤裡混著好幾群帶不同突變的細胞，我們想約束可能的局部突變狀態關係**。

傳統短讀長資料主要提供位點邊際頻率，單靠邊際頻率無法唯一識別聯合結構。本系統改用 **ONT 長讀長**：一條 DNA 分子可以**同時跨過好幾個突變位點**，所以同一物理分子上的共現是直接觀測；細胞共屬、突變先後與譜系仍是 model-dependent inference，不能直接觀測。

> **宣稱邊界**：exact-PS 主線產出的是 **local recurrence-allowed minimum mutation-state candidate arborescence／topology**。read-AF 是未經 CN／LOH 校正的條件式排序；甲基化只做 pattern-conditioned association。它不確認 cellular clone、subclone 或 biological ancestry。

系統由 **兩支 C++ 主程式**（InterSubMod、LongLineage）＋ **一層上游工具鏈** ＋
**commit-pinned Python research solvers** ＋ **validator／HTML 呈現層**組成。其中只有一部分
現在真的跑得動；Python 是實作語言，不是 evidence tier，producer 與 presentation 必須分開記錄。

> **✅ 本系統最值得注意的設計**
> 它在**機器欄位的層級**把「技術跑通」和「科學可用」分開：canonical `cohort_receipt.json` 承載 `technical_all_pass = true`，而 `summary/all7_summary.json` 承載 `validation_evidence_eligible = false`；兩者不是 `authority_manifest.json` 的 top-level 欄位。
> Frozen cohort receipts 的具名 checks 對其引用的 frozen artifacts 通過；這不是 current-source certification，也不是 production／release gate PASS，因此系統仍宣告這批結果**還不能當作驗證證據**。
> 這種「不確定就明講不確定」的紀律，貫穿整份設計。

---

## 01 · 全景圖 — 資料從哪裡進來，經過誰，變成什麼

由下而上共五層。方框顏色 = 現在的可用狀態：✅ 可跑 · ⚠️ 有限制 · 🔴 被鎖住。

![architecture-overview](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/0e1822b5957cb4bf8720f0f6ce84f5d0297ef1df/docs/images/architecture-overview.png)

> **圖 1 · 五層全景與資料流**
>
> **怎麼讀這張圖：** 由下往上是資料流方向。HP/PS 有兩種 commit-scoped representation：
> InterSubMod 從 BAM aux tags 讀 HP/PS；exact-PS／LongLineage 使用 sidecar。兩者是平行
> provenance contracts，不是兩支引擎直接串接的執行期介面；feature/preview 的 tagged-BAM
> writer 能力另依明示 branch/commit 判讀。
>
> **數字出處：** 292 GB tumor BAM 與 40,859,727 列 sidecar 是具名歷史檔案量測。2026-08-13 盤點識別 7 個 PRE-FIX `paired_full` BAM，exact paths/bytes/mtimes 合計 1,840,983,466,353 bytes（1.67436 TiB），另有 7 個 `paired_pileup`；14 個總計 3,709,322,840,333 bytes。Sampled-content set identity 已凍結，但未讀全檔 SHA-256、generation-level correspondence 未閉合，全為 `PARTIAL`/`NON_FINAL`。`paired_full` 除以現行 7 sidecar exact stat bytes 為 **294.2669×**，只是跨世代 storage-footprint quotient；不是因果壓縮效果或 frozen authority，舊 287× 錯誤。

### 五層逐層說明

#### ① 原始資料

| 部件 | 說明 |
|---|---|
| `tumor BAM`（ONT） | HCC1395 = 292 GB |
| `normal BAM` | 配對正常組織 |
| `reference FASTA` | GRCh38 + `.fai` 必備 |

← 這三樣是所有分析的共同起點；BAM 需自帶 `MM`/`ML` 甲基化標記。

#### ② 上游工具（都是外部既有軟體，不是本專案寫的）

| 工具 | 狀態 | 做什麼 | 產出 |
|---|---|---|---|
| **Dorado basecaller** | ✅ | 電訊號 → 鹼基 + 甲基化 | `MM` / `ML` tag |
| **ClairS v0.4.0** | ✅ | 腫瘤/正常配對比對 | somatic SNV VCF |
| **LongPhase-S** | ✅ | 判斷每條 read 來自哪一套染色體 | `HP` / `PS` tag |
| **SAVANA（cna）** | ⚠️ | 拷貝數變異 | 目前 `NOT_INTEGRATED` |

前三個是必要輸入的生產者。⚠️ SAVANA 的 CN 結果尚未接進主線。

#### ③ Read-label representations — BAM aux tags 與 sidecar 是平行契約

**HP/PS sidecar TSV**（`.tsv.gz` + `.tbi` 索引）

```text
#CHROM  START0  END0  QNAME  FLAG  MAPQ  CIGAR_B2  HP  PS
```

每條 read 一列，記錄它落在哪、屬於哪個單倍型。HCC1395 這一份有 40,859,727 列、1.43 GB。

InterSubMod 的 supported core input 是含 HP/PS aux tags 的 BAM，**不直接讀這份 sidecar**；
sidecar 是 exact-PS／LongLineage 的輸入契約。LongLineage preview 的 tagged-BAM writer 則是
branch/commit-specific capability，不能由檔名推成 public supported workflow。

為什麼用 sidecar：FIFO mode 可避免為該次 run 新增落地 BAM，sidecar只保留 frozen 九欄 contract。這不代表歷史 tagged BAM 從未落地；7 `paired_full` + 7 `paired_pileup` PRE-FIX inventory 仍是 `PARTIAL/NON_FINAL`，不可將 294.2669× footprint quotient 寫成壓縮效果。

#### ④ 兩支 C++ 主程式（本專案核心）

| | **InterSubMod（`inter_sub_mod`）** ✅ 可跑 | **LongLineage（`longlineage`）** ⚠️ 部分可跑 |
|---|---|---|
| 定位 | 研究原型。對每個 somatic SNV 開一個視窗，把跨過它的 read 的甲基化狀態排成矩陣 → 算 read-read 距離 → 分群 → 統計檢定。 | **PRIVATE research-preview candidate**；source-origin／license／dependency audit pending，`NOT_READY`/non-production，P3/P4/P5/P7/P8 BLOCKED。Schema/receipt 與 fail-closed 是工程契約，不證明生物正確性。 |
| 輸出粒度 | **per-region**（每個位點一包結果） | **per-read**（每條 read 一列） |
| 現況 | frozen release baseline `ddd8909a` 的釘選 source census 為 29 個 CLI 選項；此數量 branch-sensitive，不是穩定 API。build output 不入版控，fresh clone 需 exact checkout + clean build。2.9 秒只是一份具名 internal fixture receipt，不是 runtime promise。 | 🔴 candidate `b9aaa12` 的 `run` / `probe` 被鎖；frozen HCC1395 receipt 中只有 `dataset-gate` 觀察到 research artifacts |
| 吃 | tumor BAM + reference + somatic VCF（必要） | 8 角色 manifest（BAM + VCF + sidecar + reference） |
| 吐 | `methylation.csv` · `distance_matrix` · `tree.nwk` · `significance_summary.csv` | `site_reads` · `methyl_calls` · `cooccurrence_*` · `topology_unit` |

#### ⑤ Python science producers 與 validator／HTML 呈現層

| 部件 | 說明 |
|---|---|
| **Commit-pinned research solver／分析腳本** | exact-PS 共現骨幹建構 · local recurrence-allowed minimum mutation-state candidate family 枚舉；只有 producer commit、input、command receipt、schema 與 hash 齊全時才可作 science evidence |
| **Validator／publication builder／standalone HTML** | 驗證既有 JSON/TSV、明示分母、產生圖表資料與離線人工判讀頁；不得暗中重跑或改寫 science |

Python research solver 是獨立 science producer；validator／publication builder 與 HTML 才是呈現路徑。
HTML 只能讓人檢查 validated artifact，人工判讀可由 `localStorage` 匯出，但不能反過來取代
producer receipt 或 authority data。

⚠️ 這層目前最碎片化，見「分析與呈現層」分冊。

---

## 02 · 兩支程式是什麼關係？（和直覺不同）

最容易誤會的一點：LongLineage **不是** InterSubMod 的新版本，也不是它的分支。

| | InterSubMod（ISM） | LongLineage（LL） |
|---|---|---|
| **定位** | 研究原型：C++ 核心 ＋ 大量 Python 探索腳本 | PRIVATE research-preview candidate：truth-isolated、fail-closed；source-origin／license／dependency audit pending |
| **git 起源** | `71ce312e` | `7a8f75e8` — **完全不同的根 commit** |
| **血緣** | — | `ORIGIN.md` 記錄 21 條 source mapping 與 17 條 authority；在 origin/license gate 完成前，不把 clean-room／可公開來源當成已驗證結論 |
| **輸出粒度** | **per-region**（每個 SNV 位點一包） | **per-read**（每條 read 一列） |
| **單倍型來源** | 直接從 BAM 的 aux tag 讀 | **明文禁止讀 BAM tag**，只認 sidecar TSV |
| **對「甲基化能不能參與重建」的立場** | 🔴 **絕不參與** — 只做事後關聯描述 | ⚠️ **是必要關卡** — 甲基化分不出群的位點不會有拓撲解 |

> **🔴 兩系統的科學立場直接衝突**
> ISM 的既有結論是「甲基化**絕不**進入重建的 likelihood」，因為會產生循環論證（見下方 §04）。
> 但 LL 現行程式碼把「甲基化有沒有分出群」當成拓撲單元的**前置關卡**。
> 後果：在 LL 裡，**甲基化沒分出群的位點永遠不會有拓撲解** —— 這是一個選擇性偏誤來源。
> 任何跨系統的比較都必須明講這點，否則會得出與 ISM 既有結論矛盾的敘述。

---

## 03 · 現在到底哪些能跑？—— 誠實狀態表

這張表是本頁最實用的部分。所有狀態都是本輪**實際執行過**或**讀原始碼確認**的結果，不是文件宣稱。

| 部件 | 狀態 | 實測依據與限制 |
|---|---|---|
| **`inter_sub_mod`**<br>ISM 主程式 | ✅ 可跑 | Build output **不入版控**；fresh clone 必須從 exact checkout 在 clean build directory 建置。具名 internal fixture 曾有單次 **2.9 秒、exit 0** receipt，僅限該硬體／commit／輸入，不是讀者 runtime promise。 |
| **`longlineage preflight`** | ✅ 可跑 | 驗證 8 個角色的 manifest 完整性。 |
| **`longlineage dataset-gate`** | ⚠️ 限制 | 在 private candidate `b9aaa12` 的 frozen HCC1395 receipt 中，唯一觀察到能產生 research artifacts 的介面；該次產出 8 個 artifact，不等於 production science validation。<br>🔴 樣本名稱硬寫死在該 revision 原始碼裡。 |
| **`longlineage run` / `probe`** | 🔴 被鎖 | `KernelBlocked`（exit 6）是預期 fail-closed；`release_attestation.json` 固定 P3/P4/P5，phase ledger／public-safety receipt 固定 P7/P8，五者仍 **BLOCKED**；整體 `NOT_READY`／non-production。 |
| **`longlineage` 的 topology 產出** | 🔴 一份 receipt 為 0 candidate units | Candidate `b9aaa12` 的 frozen HCC1395 dataset-gate receipt：79,687 個位點 → 通過甲基化關卡 5 個 → mutation-state candidate topology units **0 個**。134,278 個配對中 134,276 個不合格；不可外推其他 revision、dataset 或 run。 |
| **`longlineage-query` / `export-legacy`** | 🔴 未實作 | 仍是 fail-closed 的空殼。 |
| **Research exact-PS topology-AF solver**<br>Python runners + separate C++ solver | ✅ 可跑 | 這是目前實際產出 2026-07-24 funnel 數字的路徑，7 個資料集全跑完；不是 `inter_sub_mod` binary（見下方 §05）。 |
| **7 個 technical datasets／6 個 biological IDs 的 sidecar 資料** | ✅ 齊全 | 實測 `run_status.tsv` 顯示 **7/7 technical datasets PASS**；HCC1395 與 HCC1395_DORADO 是同一 biological ID 的技術資料集。 |
| **輸出帶標籤的 BAM** | ⚠️ 依版本而異 | `inter_sub_mod` 不寫 BAM；LongLineage private baseline/main snapshot **`5daf50f`** 也沒有 writer；private public-preview candidate **`b9aaa12`** 含 `longlineage-tag-bam`。這是 revision-scoped、尚未發布的 `NOT_READY` capability；現存檔案來源仍須依 provenance 判定。 |
| **兩支程式串起來跑** | 🔴 無 | **沒有任何腳本同時執行兩個 binary**（雙向搜尋皆無命中）。目前是兩條各自獨立的線。 |

> **🔴 給教授看的三句重點**
> 1. **exact-PS funnel 來自獨立的 research `exact_ps_topology_af` C++ solver ＋ Python runners**；`inter_sub_mod` 產生 per-region 甲基化／統計輸出。兩者不可寫成同一支程式。
> 2. **LongLineage 的 `P3/P4/P5/P7/P8` 具名 gates 尚未通過**——包含 parity／validation、source-origin、license、dependency 與 release-safety；單一 dataset-gate 執行部分 kernels 不等於 production 或公開資格完成。
> 3. **沒有一條打通兩支程式的路**。零件很齊、契約很嚴，但接縫還沒建。

---

## 04 · 為什麼現行流程不用甲基化建立 candidate graph？

這是整套方法論最核心、也最常被外部讀者質疑的設計決定，值得單獨解釋。

直覺上會想：既然甲基化模式在不同細胞群之間有差異，那用甲基化把 read 分群，不就能找出亞群了嗎？
**現行形式不行 —— 若用甲基化定義群，再用同一甲基化訊號驗證群，會產生循環論證。這不是對所有可能甲基訊號或新設計的總檢定。**

![methylation-circularity](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/0e1822b5957cb4bf8720f0f6ce84f5d0297ef1df/docs/images/methylation-circularity.png)

> **圖 2 · cis-ASM 循環：為什麼「用甲基化確認亞群」證明不了東西**
>
> **這個設計的代價與收穫：** 代價是甲基化這個資訊量很大的軸被降級成輔助；收穫是避免「用甲基定義群、再用甲基驗證群」這一項可預見的循環論證批評，不代表已通過外部同儕審查。
> 實際數字上，在 frozen formal family 中，甲基化輔助層共 1,045 個正式單元、811 個可評估、最終只有 **3 個**（0.37%）達到該定義下的 robust pattern-conditioned association。這是現行設計不把該 formal family 當骨幹的一個理由，不是所有 methylation signal／backbone design 的總檢定。

### 觀察到「兩群甲基化模式不同的 read」時，至少有四種成因，單一 bulk 分不開

**觀察**：同一位點的 read 分成甲基化高 / 低兩群。

| # | 成因 | 說明 |
|---|---|---|
| ① | **germline 等位特異性甲基化** | 天生兩套染色體的甲基化就不同（如印記基因） |
| ② | **LOH 解遮蔽** | 腫瘤丟掉一套染色體，原本被平均掉的差異浮現 |
| ③ | **拷貝數增益的劑量效應** | 某段被複製多份，比例被拉偏 |
| ④ | **真正的亞群譜系差異 ← 我們想要的** | 不同細胞群確實帶不同甲基化 |

**🔴 循環論證**：要用甲基化「確認」某群 read 是亞群，得先排除 ①②③。但要排除 ①②③，就得先知道哪些 read 屬於哪個亞群 —— **這正是待證明的結論。** 在目前 single-bulk measurements 下，若沒有 orthogonal data 或額外假設，四種成因不可識別。

> **✔ 本系統的解法：改用「同一條分子上的突變共現」當骨幹**
> 一條 read 上同時看到兩個候選 allele = 物理上的分子連鎖。它只在不使用甲基化衍生標籤的限定下非循環；仍依賴上游 variant calling、alignment、basecalling 與 haplotag 假設。
> 甲基化被結構性地隔離在 exact-PS 主線之外：local recurrence-allowed minimum mutation-state candidate topology 先由遺傳證據定好，甲基化只在事後做 association-only 描述，改不動任何一條邊。

---

## 05 · 實際跑出來長什麼樣？—— 7 technical datasets／6 biological IDs funnel

這是目前的 canonical 結果（2026-08-01 凍結，7 個 technical datasets／6 個 biological IDs、chr1–22）。每一層的數字都能在 `denominator_registry.tsv` 或 `all7_summary.json` 裡查到。

| 指標 | 數值 |
|---|---|
| sSNV dataset records（7 technical datasets／6 biological IDs） | **469,849** |
| k=1 strict read-linkage components | ⚠️ **170,131 / 255,752（66.52%）** |
| frozen read-AF model 可排序的 candidate units | ✅ **71,955** |
| confirmed 細胞亞群 | 🔴 **0** |
| confirmed 線性祖先關係 | 🔴 **0** |

![funnel-7samples](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/0e1822b5957cb4bf8720f0f6ce84f5d0297ef1df/docs/images/funnel-7samples.png)

> **圖 3 · 從 469,849 筆 dataset records 到 63,506 個 model-assigned candidate-shape records：每一層如何改變分析 grain**
>
> **算術自洽性已驗證：** 39,648 + 23,858 = 63,506；+ 8,449 = 71,955；+ 3,224 + 45 = 75,224；+ 10,717 = 85,941；+ 13,014 = 98,955。每一層加總都對得起來。
>
> **🔴 最常被誤引的數字就是 88.2579%。** 正確讀法是：frozen model 對 **63,506 / 71,955 個 rankable candidate units** 指派單一 rooted-unlabeled candidate-shape signature。這是 **local、recurrence-allowed、model-conditional** 的圖形統計，不是 biological topology、accuracy、prevalence 或「腫瘤演化關係已解出 88%」。另有 170,131 / 255,752（66.5219%）個 strict components 為 k=1；該比例不能套用到 469,849 筆 dataset records。

### funnel 逐層數字

| 層 | 內容 | 數值 | 該層流失 |
|---|---|---|---|
| L1 | sSNV dataset records（7 technical datasets／6 biological IDs · chr1–22） | **469,849** | — |
| L2 | 嚴格 read-linked 連通成分 | **255,752** | ⚠️ k=1 strict components：**170,131 / 255,752（66.52%）**；這是 component grain，不是 dataset-record grain |
| L3 | k≥2 有共現 → 有界切分後的分析單元 | **85,621 → 98,955** | k>12 的大成分被切開，所以單元數會變多 |
| L4 | 帶突變的單元 | **85,941（86.85%）** | ⚠️ 無活躍 ALT，不推論：**13,014（13.15%）** |
| L5 | local recurrence-allowed minimum mutation-state candidate family **完整枚舉完成** | **75,224（87.53%）** | 🔴 算力上限，主動放棄：**10,717（12.47%）**<br>非隨機缺失，集中在複雜區 |
| L6 | 可用 read-AF 排序 | **71,955（95.65%）** | ⚠️ 分母為 0 / 遞迴篩選：**3,224 + 45** |

### L6 之後的三分支

| 分支 | 數值 | 說明 |
|---|---|---|
| ✅ 唯一 AF-best candidate arborescence | **39,648（55.10%）** | 未經 CN／LOH 校正的 read-AF 算術上唯一勝出；不是確認的生物譜系 |
| 並列，但 candidate shape signature 相同 | **23,858（33.16%）** | 並列 candidates 在 frozen model 下共享同一 rooted-unlabeled shape signature；非 biological topology／accuracy／prevalence |
| ⚠️ 並列且跨不同樹形 | **8,449（11.74%）** | 「定不出來」← 這也是答案 |

**frozen model 指派單一 rooted-unlabeled candidate-shape signature = 39,648 + 23,858 → 63,506 / 71,955 rankable candidate units = 88.2579%**

> **🔴 排序器同時也是篩選器 — 一個容易被忽略的細節**
> read-AF 不只用來「排序」候選樹，它在同分的候選集合裡**取最大值並丟掉其餘**，因此它實際上也**改變了候選集合本身**（實測丟棄 30.03%）。
> 後果：**純 parsimony（不靠 AF）的單一拓撲率只能說是一個區間，不能報單一數字**。
> ⚠️ 本輪查證發現：該區間的上界 88.26% 有 registry 佐證，但下界 64.89% **在本 repo 內未 grep 到原始出處**，僅存在於工作記憶索引 —— 對外引用前必須先補上出處。

---

## 06 · 能回答什麼、不能回答什麼

這節直接取自 `authority_manifest.json` 的 `claim_boundary` 欄位 —— 是機器可讀的宣稱邊界，不是散文。

### ✔ 可以主張（permitted）

- 嚴格 read-linked 的**局部區域**結構
- 家族完整時，完整的最小 read-相容突變狀態候選家族
- 允許遞迴的 Hamming-1 父子樹候選
- **未經 CN/LOH 校正**的 deterministic read-AF 排序
- 全域 AF-best 數學樹的 exact rooted-unlabeled 拓撲普查
- pattern-conditioned 的區域甲基化**關聯**（非因果）

### ✘ 明文禁止主張（forbidden）

- **confirmed 細胞層級的 clone / subclone**
- CN/LOH 校正過的 CCF
- 校準過的 likelihood 或 posterior
- 把 read cluster 當成 cell clone
- 用甲基化差異重排拓撲
- 把 HCC1395 與其 DORADO 重跑當成生物學重複

### 已被判定為死路的方向（不要重開）

| 方向 | 為什麼不行 |
|---|---|
| **對整棵樹做機率排序** | 候選樹在 likelihood 下**等機率**，posterior 只會吐 P=1/N。換成連續 posterior「只是換個記法，變不出 reads 裡不存在的資訊」。 |
| **用甲基化排序拓撲** | 三個理由全否：cis 構造對稱、四配子最常見成因是拷貝數多重性而非二次突變、等機率區沒有遺傳真值可校驗（等於拿甲基驗甲基）。 |
| **甲基化非監督分群確認亞群** | 單一位點的四種成因不可拆（見 §04），屬已判 NEGATIVE 的 tumor-only double-dip。 |
| **甲基化當變異 TP/FP 過濾器** | 已測 formulation 為 concluded negative；只有提出 materially new hypothesis 並先通過 pre-decision audit 才可重開。 |
| **用舊的 50 kb 鄰近取代嚴格 read linkage** | 幾何鄰近不等於同分子連鎖。已被 read-linked 方法取代。 |

> **證據等級的天花板**
> 在目前單一 bulk observation／model 且未整合 CN／LOH 的條件下，不能確認 cellular
> subclone 或 linear ancestry，因此目前封頂 ★3。single-cell、multi-region、orthogonal
> CN／purity 等額外獨立證據可能提高識別性；不宣稱只有某一種 assay 可行。Steiner
> 節點會被誠實標成 `inferred`。

---

## 07 · 接下來看哪一頁

依你的目的挑一頁進去。

| 分冊 | 狀態 | 一句話 | 吃 → 吐 |
|---|---|---|---|
| [ISM 部件與輸入輸出](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) | ✅ 可跑 | `ddd8909a` 的釘選 29-option census、內部 8 個處理階段、每個輸出檔的真實欄位與範例行；option count 是 branch-scoped。 | tumor BAM + reference + somatic VCF → `methylation.csv` · `distance_matrix` · `tree.nwk` |
| [LongLineage 部件與契約](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) | ⚠️ 部分可跑 | 4 個子命令的真實狀態、M1→M2→topology 科學鏈、artifact schema 與 `read_id` 的算法。 | 8 角色 manifest → `site_reads` · `cooccurrence_*` · `topology_unit` |
| [上游工具鏈與資料](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) | ✅ 齊全 | sSNV / HP-PS / 甲基化 / CN 各自從哪個工具來、7 technical datasets／6 biological IDs 的實體資料在哪。 | Dorado · ClairS · LongPhase-S · SAVANA · sidecar TSV |
| [Python 分析層與 HTML 呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) | ⚠️ | 主力腳本各自做什麼、工作站 HTML 怎麼生成、這層目前的碎片化問題。 | 分析腳本 → workstation spec → standalone HTML |
| [怎麼跑 — 操作手冊](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run) | 可複製指令 | 從編譯到跑出第一個結果的完整步驟，每條指令都已實際執行驗證過。 | build → 最小跑法 → 完整跑法 → 看結果 |
| [科學方法細節](https://github.com/liaoyoyo/InterSubMod/wiki) | 既有頁面 | 名詞定義、ISM 方法本體、甲基化讀取與篩選、chr2:18M 真實案例逐步拆解。 | 解釋中心第 01–10 頁 |

---

## 本頁的資料來源與驗證方式

本頁所有數字、狀態、檔案格式，來自 2026-08-06 對兩個 repo 的一次系統性實測收集（9 個維度並行，每個維度獨立收集後再做對抗式驗證）。

- **歷史自述、目前 `UNVERIFIED`**：2026-08-06 文件曾記錄「334 個 source refs＋111 個
  paths，0 missing／0 out-of-bounds」。公開 repo 目前缺 commit-pinned inventory、command log、
  hashes 與 replay receipt，因此不能把該自述當成 current blanket guarantee；ALG-023 維持
  `UNVERIFIED`，待可重播 receipt 補齊才可升格。
- **可跑狀態**：所有標示「可跑」的部件，均為本輪**實際執行並檢查 exit code** 的結果。
- **科學數字**：funnel 各層取自 `denominator_registry.tsv` 與 `all7_summary.json`（2026-08-01 凍結的 canonical 輸出），並已驗證各層加總自洽。

### 已知的資料缺口（誠實標註）

- 「純 parsimony 單一拓撲率」的**下界 64.89%** 在 repo 內未找到原始出處，對外引用前需補齊；上界 88.26% 有 registry 佐證。
- LongLineage 對 7 technical datasets／6 biological IDs cohort 的總執行時間與記憶體上界**無實測證據**，且其文件明文禁止由單一 dataset 外推。
- 拷貝數（CN／LOH）目前狀態為 `NOT_INTEGRATED`，尚未接入主線。

**來源**：`InterSubMod/docs/explain/11_system-map-overview.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06
