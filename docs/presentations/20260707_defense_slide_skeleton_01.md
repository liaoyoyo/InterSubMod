<!--
建立時間: 2026-07-07
類型: 口試 PPT 逐頁骨架（規劃階段；講稿待逐頁填）
狀態: draft / planning（結構已與使用者對齊；講稿 + 個案數字待下一階段填 + 核實）
build_branch: research/subclonal-reconstruction-202606
資料來源:
  - docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json  (聚合真值，本 session 讀回)
  - docs/methodology/_assets/20260706_layered_reconstruction_workstation.standalone.html    (P15 工作站，2026-07-07 建)
  - docs/methodology/_assets/layered_workstation/{7 樣本}.html                              (P16 跨樣本可攜版)
  - docs/references/20260701_thesis_research_handoff_onboarding_01.md                        (主軸/framing SoT)
-->

# 口試 PPT 逐頁骨架 — Subclonal reconstruction using somatic haplotagging and methylation profiles

> **本檔用途**：口試簡報的**逐頁規劃工作面**。結構（頁序 / 必要 vs 補充 / 圖）已與使用者對齊；`講稿{{待填}}` 與 `預期提問{{待填}}` 兩格為下一階段逐頁填寫。
> **狀態**：⭐3 exploratory / proof-of-concept / 單樣本 HCC1395 主 + 7 ONT 樣本比例佐證 / single-pipeline。
> **⚠ 數字紀律（§13.0）**：標「✅已驗證」的聚合數字＝本 session 從 census json 讀回；標 `⚠核實` 的個案/跨樣本數字＝來自 memory（較舊），**投影片動筆前必重新 grep**，不先填預期值。

---

## 0. 全場設定

| 項目 | 值 |
|---|---|
| 時長 | **50 分鐘**（內容 ~46 min + 緩衝 4 min；Q&A 另計） |
| 聽眾 | 生物資訊背景、懂基因體學的教授 |
| 頁數 | **20 頁核心**（方法/結果可各拆 1–2 頁 → 上限 ~24） |
| 敘事框架 | macro = ABT + SCQA；Results = Rule(脈絡)+Exception(exemplar climax)；標題 = Assertion-Evidence；結論 = Verdict-Pyramid |
| 深度分配 | **背景講快**（聽眾已懂）→ **深度壓在方法段 P6–P11** → 結果兩個高潮各歸各訊號源 |

### ✅ 已驗證聚合真值（本 session 從 `layered_reconstruction_HCC1395.json` 讀回，可直接引用）

| 指標 | 值 |
|---|---|
| regions_total | **7100** |
| regions_mixed_germline（multi-HP） | **3992（56.2%）** |
| family units / lineage units / germline units | **17073 / 12475 / 10988** |
| somatic3 units / none units | 1487 / 4598 |
| recurrence → 送 CN | **43** |
| V1–V7 驗證 | **全 pass、0 verify-fail** |
| L3 甲基層狀態（工作站內） | **pending（bounded-auxiliary；非循環；不 rank 樹集）** |
| 工作站可攜版 | 7 樣本齊，2026-07-07 02:40 建 |

> `⚠核實` 清單（memory 07-07，投影片前重 grep）：all-det 32.6% · avg 樹 3.61 · 7 樣本跨表（COLO829 78.7% multi-HP 等）· chr2:18M 個案 V/F/reads · 甲基 Δβ 0.25–0.33 · Q5 7 候選 · HP contribution 9187/5238。

---

## 1. 🔴 全場防守校準（Defense-critical — 每頁都受此約束）

### 校準 A — 甲基「排序」用詞（P8/P9 核心）
- 遺傳**已定**的樹 → 甲基做**一致性佐證（corroborate）**＝強、可辯護（cn-clean Δβ `⚠核實` 0.25–0.33）。
- 遺傳**真無法定**的 attachment → 甲基只給**候選提示（candidate / exploratory）**，需正交確認。
- 依據：four-gamete demo（07-03）乾淨 same-HP 案例從不顯著（min p=0.066, n=6, underpowered），唯二顯著全 CN-gain confound。
- **講詞**：說「甲基**佐證 / 一致性支持**」，**不說**「甲基**排序 / 確認**樹結構」。

### 校準 B — 「HP 排除 cis-ASM」講到位
- HP 條件化移除 **germline HP1-vs-HP2 ASM**（62% neutral、Δβ 0.97 大訊號）✅ 方法強項，必凸顯。
- 殘餘 = **somatic-cis**（甲基是該 somatic 突變的 cis 下游）→ 與 subclone 綁定、**非完全獨立** → 用「cis-clone / cis-subclone」誠實用詞；完全分離需 **matched-normal**（COLO829 normal 缺）。
- **講詞**：甲基是 **corroborate（同管線 cross-check）**，非 independent confirmation。

### 🔴 全場禁區（越線即被打臉）
1. 不稱「genome-wide tree」→ 只講 **regional partition（≤ read-span）**。
2. 不稱「甲基偵測 subclone」→ 甲基 **corroborate 非 detect**。
3. 不稱「分子證據＝single-cell confirmation」。
4. 不稱「對手只用二代定序 / 對手缺顯著性檢定」（cvlr/ASMS/MethylBERT 都 ONT + 有 randomization 檢定）。
5. HP1-1/HP2-1 ≠ 已確認 subclone（haplotag 類別）。

---

## 2. 逐頁骨架

> 圖標籤：**【示意】**概念圖(image-gen,無真數字) · **【實據】**真資料圖(需驗數字) · **【既有】**已存在可重用 · **【需製】**要新做 · `⚠核實`=數字用前重 grep。
> 🎯必要 = 一定要講/放 · ➕補充 = 有時間/深度才加，可砍。

### 一、背景段 P1–P5（~10 min）— **生物 block(P2–P3) → 技術/資料 block(P4–P5)**；講快，只留能導出方法的點
> 2026-07-08 對調 P3/P4：先堆滿生物 need（二次打擊→演化→subclone），ONT 才作為 answer 登場；技術+資料+gap 三頁連續。

#### P1 標題 ⏱1'
- 🎯 問好 + 題目 + 姓名 / 指導教授 / 日期
- 🖼【既有】標題頁（無圖亦可）
- 講稿{{待填}} | 預期提問{{待填}}

#### P2 癌症與二次打擊 ⏱2'
- 🎯 germline first-hit + somatic/LOH/**甲基** second-hit → 抑癌基因雙等位失活 → **所以要看「每個單倍型上的突變＋修飾」**（通往方法的橋）
- ➕ SNP vs sSNV、致癌 vs 抑癌基因細分
- 🖼【既有】第 2 頁圖，**修「倍體」→「單倍型」**；定序表移去 P3
- 講稿{{待填}} | 預期提問{{待填}}

#### P3 癌症演化 + clone/subclone ⏱2.5'
- 🎯 (a) main clone vs subclone 定義；(b) 為何 subclone 重要（抗藥/復發/轉移）；(c) **HP1-1/HP2-1 是 haplotag 類別，≠已確認 subclone**（埋防雷）
- 🎯 收尾一句 need：「要**在每個單倍型上、跨演化**看突變＋甲基＋它們的**共現**」→ 為 P4 ONT 鋪路
- ➕ CCF 概念、bulk 混合細胞挑戰
- 🖼【既有】第 3 頁演化圖
- 講稿{{待填}} | 預期提問{{待填}}

#### P4 為何 ONT ⏱2'
- 🎯 承上 need → 需要能**一次讀全且保留單分子共現**的工具 → **ONT 單一長讀同時帶 SNP + sSNV + 5mC**，不需拆多次實驗
- ➕ 定序比較表（長度/準確度/價格）當佐證，別逐格唸
- 🖼【既有】定序比較表 `⚠核實`（讀長/Q值/價格標 typical range）
- 講稿{{待填}} | 預期提問{{待填}}

#### P5 資料 + 與過去方法差別 ⏱2.5'
- 🎯 (a) ONT bulk Tumor/Normal；(b) 白地 = **read 層甲基結構檢定 + normal-anchored cis-test + somatic-subclone 三者無人組合**
- 🔴 禁區：不講「對手只用二代 / 缺檢定」
- 🖼【需製】方法對照表 or【示意】白地定位圖
- 講稿{{待填}} | 預期提問{{待填}}

### 二、方法段 P6–P11（~14 min）— 深度集中，真貢獻也最易被電

#### P6 方法總覽 ⏱2'
- 🎯 **L0→L1→L2→L3 一張流程圖**：前置→longphase-S→HP家族→sSNV共現樹→CN→甲基（先全景再拆解）
- 🎯 **資料特性 → 方法能力（兌現 P4 的 ONT 承諾）**：①長讀單分子→一 read 多 sSNV→**共現**(L1) ②Tumor/Normal→somatic vs germline 可分(L0) ③5mC 同讀→甲基不需額外實驗(L3)
- 🖼【需製】手刻流程圖（零依賴 SVG）— **全場最重要導覽圖**；流程圖旁標「資料特性」三箭頭
- 講稿{{待填}} | 預期提問{{待填}}

#### P7 L0 前置 + HP 家族分類 ⏱2'
- 🎯 (a) longphase-S 出 read HP + 帶 somatic 的 read；(b) **per-HP-家族分樹 = 分 allelic vs clonal = 命脈**；(c) ✅7100 區、3992（56.2%）multi-HP
- ➕ longphase-S 為何非循環（germline-anchored、somatic 獨立、零甲基）→ 備援
- 🖼【示意】HP 家族分群（region 拆 HP1 家族樹 + HP2 家族樹）
- 講稿{{待填}} | 預期提問{{待填}}

#### P8 L1 sSNV 共現建樹（核心，可拆 2 頁）⏱3'
- 🎯 (a) 同 read 多 sSNV **2×2 共現 → 互斥/巢狀/共連**；(b) 布爾超立方體 **有向有根群組斯坦納樹枚舉全最小樹**；(c) **N 樹 = M 形狀**（枚舉非最佳化、非 rank）
- 🎯 **取樣邊界（誠實框架，先講）**：①**read span → regional partition**（只能連一 read 長度內的 sSNV = scope 上限、非 genome-wide，同時擋禁區①）②coverage 是否觀察到中間群 → determined vs ambiguous 是**取樣邊界非生物不可知**③partial-read（覆蓋 sSNV1&3 非2 仍供 1-3 共現）→ Pe'er IDP 救援
- ➕ 計算複雜度定位（多項式島 vs NP-hard）→ **Q&A 強備援**
- 🖼【示意】超立方體 + 2×2 → 樹；【實據】真實 region 枚舉（工作站截圖）
- 講稿{{待填}} | 預期提問{{待填}}

#### P9 L3 甲基輔助 ①corroborate + ②attachment ⏱3'
- 🎯 (a) HP 排除 germline-ASM → 殘餘 cis-clone/subclone；(b) 甲基**一致性佐證**已定樹（Δβ `⚠核實`）；(c) ambiguous attachment 只給**候選提示**（four-gamete 測過→誠實界定，不排序）
- 🔴 校準 A + B 全在此頁；用詞「佐證/一致」非「排序/確認」
- 🖼【實據】Δβ 對齊圖 `⚠核實`；【示意】HP 排除 germline 軸
- 講稿{{待填}} | 預期提問{{待填}}

#### P10 L3 甲基輔助 ③單位點 ⏱2'
- 🎯 無 sSNV 共現時，**單位點甲基雙峰旗標候選**隱藏 subclone（Q5 `⚠核實` 7 候選，需正交確認）
- ➕ 單位點 ASM 為何「最不可解釋」（無錨 + C>T 破 CpG）
- 🖼【實據】單位點雙峰分布 `⚠核實`
- 講稿{{待填}} | 預期提問{{待填}}

#### P11 完整標註狀態總表 ⏱2'
- 🎯 各狀態處理：**all-det 32.6% `⚠核實` / has-ambiguous / capped / recurrence 43→送 CN**（見樹也見林四層）
- 🎯 **determinacy 是 sampling-dependent**：ambiguous ≠ 生物不可知，多為取樣未捕捉中間群 → 更深 read / single-cell 可再解（preempt 口委「為何只 32.6%」）
- ➕ capped（>4 hidden）+ GCAP=8 截斷誠實標
- 🖼【實據】狀態分佈表（✅ 7100 / 12475 / 43 / V1-V7 已驗）
- 講稿{{待填}} | 預期提問{{待填}}

### 三、結果段 P12–P17（~15 min）— 兩個高潮各歸各訊號源

#### P12 HCC1395 aggregate ⏱2'
- 🎯 7100 區、multi-HP 56.2%、all-det 32.6% `⚠核實`、avg 樹 3.61 `⚠核實`
- 🖼【實據】比例長條（聚合已驗）
- 講稿{{待填}} | 預期提問{{待填}}

#### P13 甲基 exemplar 高潮（chr2:18M）🎯甲基高潮 ⏱3.5'
- 🎯 單案完整走 L0–L3，**證方法正確**：甲基分群↔subclone 譜系顯著、↔het allele NS
- ➕ 第二獨立窗（18096269）甲基↔somatic 等位 V 高
- 🖼【既有/實據】chr2:18M 熱圖 `⚠核實`（V/F/reads）
- 講稿{{待填}} | 預期提問{{待填}}

#### P14 甲基 corroborate / Q5 候選 ⏱2.5'
- 🎯 cn-clean 高一致位點 Δβ `⚠核實`；甲基 = 有界佐證
- ➕ CN-gain confound 如何排除（rate-by-CN-state 非 count）
- 🖼【實據】個案圖 `⚠核實`
- 講稿{{待填}} | 預期提問{{待填}}

#### P15 樹枚舉工作站 = 遺傳重建高潮 🎯遺傳高潮 ⏱3'
- 🎯 **1 樹=1 形狀（data 夠）vs N 樹=M 形狀（不足→列全候選+量化不確定）**；「定不出來即答案」是正解非 bug
- 🔴 格式：**2 張靜態截圖**，live 工作站留 Q&A 備援（勿現場點擊）
- 🖼【既有】工作站（今日驗證可用；**L3 甲基 pending → 此頁只講遺傳**）
- 講稿{{待填}} | 預期提問{{待填}}

#### P16 跨樣本一致性 ⏱2.5'
- 🎯 (a) 同株一致（HCC1395 vs DORADO `⚠核實` 56.2/55.2、32.6/34.3）= **方法可重現/basecaller 不變**；(b) 異株有別（COLO829/H2009/HCC1954 比例 `⚠核實`）
- 🔴 **只比比例、無 p 值**（真 n=7 pseudoreplication）；絕對數受深度混淆
- 🖼【實據】7 樣本比例表 `⚠核實`
- 講稿{{待填}} | 預期提問{{待填}}

#### P17 誠實限制 ⏱2'
- 🎯 ⭐3 單樣本 / **regional partition ≠ genome-wide tree** / 甲基有界 / somatic-cis 需 matched-normal / single-pipeline 自我參照
- 🖼【需製】限制卡（純文字）
- 講稿{{待填}} | 預期提問{{待填}}

### 四、收尾段 P18–P20（~5 min）

#### P18 貢獻 ⏱2'
- 🎯 ①read-level bulk ONT **角色分離**重建 ②無監督 read×read 結構檢定 + normal-baseline cis-test ③**甲基角色精確界定**（含測過 tie-break）④跨樣本 + basecaller 不變性
- 🖼【需製】貢獻條列卡
- 講稿{{待填}} | 預期提問{{待填}}

#### P19 討論 + 未來 ⏱2'
- 🎯 matched-normal cis-control（B軸）/ single-cell 正交 / COLO829 補 normal / 7 樣本擴展
- ➕ 為何能把 ⭐3→⭐4（somatic-cis 疑似升級真 subclone）
- 🖼【需製】路線圖
- 講稿{{待填}} | 預期提問{{待填}}

#### P20 參考 + 致謝 ⏱1'
- 🎯 純文字
- 講稿{{待填}} | 預期提問{{待填}}

---

## 3. 圖片工作清單

| 類別 | 頁 | 圖 | 狀態 |
|---|---|---|---|
| **需新製** | P6 | L0–L3 手刻流程圖 | ⬜ |
| | P8 | 超立方體 + 2×2 共現 → 樹 示意 | ⬜ |
| | P9 | HP 排除 germline 軸 + Δβ 佐證 | ⬜（Δβ 需驗）|
| | P10 | 單位點甲基雙峰 | ⬜（需驗）|
| | P5 | 方法對照表 / 白地定位 | ⬜ |
| | P17/P18/P19 | 限制卡 / 貢獻卡 / 路線圖 | ⬜ |
| **既有可重用** | P1–P4 | 標題 / 二次打擊 / 定序表 / 演化圖 | ✅（P2 修字、定序表移 P3）|
| | P13 | chr2:18M 熱圖（research/flagship...）| ✅（數字需驗）|
| | P15 | 樹枚舉工作站截圖 | ✅（今日驗證）|
| | P16 | 7 樣本比例表（census）| ✅（逐樣本需驗）|

---

## 4. 待使用者確認 / 思考的問題（回填於此）

1. 必要/補充切分是否認同（背景快、深度壓方法段）？→ {{待填}}
2. P13 旗艦個案用 chr2:18M？還是換別的 exemplar？→ {{待填}}
3. 需製圖清單優先序？誰先做？→ {{待填}}
4. 先鑽哪頁寫講稿（建議 P6 總覽 + P8 建樹）？→ {{待填}}
5. 其他使用者想到的細節問題：→ {{待填}}

---

## 5. 下一階段工作流

1. 使用者 review 本檔 → 於 §4 + 各頁 `預期提問{{待填}}` 標註問題。
2. 逐頁填 `講稿{{待填}}`（演講者 3 句話）+ `預期提問{{待填}}`（口委可能問 + 你的答）。
3. 個案/跨樣本數字動筆前**重新 grep 核實**（§13.0；⚠核實 標記處）。
4. 需製圖依 §3 清單製作（image-gen / matplotlib）。
5. 講稿 + 圖齊 → 產 PPTX（/pptx-build --from-draft）或 standalone HTML preview。
