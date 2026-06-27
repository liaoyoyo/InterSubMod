---
title: 甲基聚類結構 ↔ 標籤關聯 — 完整狀況分類盤點
date: 2026-06-17
build_branch: docs/method-comparison-ism-external-202606 (worktree ism-review-infra)
data_sources: _assets/20260616_readset_provenance/label_breakdown.json, contingency_direction.json, tn_dependency.json, level2_calibration.json, assoc_label_structure.json
status: data-verified (verify_breakdown.py 25/25 PASS)
---

# 甲基聚類結構 ↔ 標籤關聯：完整狀況分類盤點

> 輸入 = caller 輸出 somatic 位點總數（TP+FP）= **30,490 region**（`output/_wg_d1_unified`，paired run）。
> 全部數字由 `label_breakdown.py` 等腳本算出落 JSON，`verify_breakdown.py` 25/25 一致性 PASS，再注入本文（不手打）。

## §0 一頁總結（先看這裡）

**三種狀況**（甲基聚類結構 ↔ a-priori 標籤 的三種關係）：

| 狀況 | 定義（驗證步驟） | 位點數 | 佔全體 |
|---|---|--:|--:|
| **① 對齊**（一結構→一標籤）| 聚類成功(PassedGating) 且 cluster×label 列聯表顯著(GlobalP≤0.05) | **22497** | 73.78% |
| **② 一標籤多群**（多結構→一標籤）| 單一 HP 內甲基再分≥2 群（WithinHP_CleanMultigroup）| **8077** | 26.49% |
| **③ 一群多標籤**（一結構→多標籤）| 聚類成功但 cluster 不追隨標籤（PassedGating 且 GlobalP>0.05）| **510** | 1.67% |

⚠ **三狀況彼此重疊、非互斥**（如「對齊 ∩ 一標籤多群」= 6180）。乾淨不重疊的 partition 是 §2 的 **6 格表（Σ=30,490）**。

**🔴 read-set 關鍵修正（C++ 源碼逐軸證實，非 agent 回報）：**

| 軸 | read-set | 對應狀況 |
|---|---|---|
| 基礎聚類 + cluster×label 對齊（GlobalP）| **PAIRED tumor+normal**（GlobalTest 無 is_tumor 過濾）| 狀況①③ |
| WithinHP 子結構 + Subclone Δβ | **TUMOR-ONLY**（`!is_tumor continue`）| 狀況② |
| Germline ASM Δβ | **NORMAL-only**（定義 germline）| — |
| Somatic/HP Residual Δβ | **NEEDS BOTH**（tumor−normal）| — |

→ 你最在意的 **subclone 相關方向（狀況②）本來就是 tumor-only**；但 **「對齊」這個基礎分類是 paired**（混了 normal reads），不是純 tumor-only 訊號。

**Level-2 read-level 抽樣校準（300/300 對齊位點，固定 seed 42）：**
- 映射驗證 **89.7%**（read-level 重建對拍 C++ 真值）→ 重建可信。
- 標籤類別歸屬命中 **92.6%**；軸（hp/allele）命中 **73.7%**。
- **tumor-only 對齊保持**：可測 159 個對齊位點中 **51.6% 仍保持**（限 tumor reads，Cramér's V≥0.3）；median V 從 paired **0.665** 掉到 tumor-only **0.319** → 對齊強度約半數來自 normal reads。

---

## §0.5 以位點為單位：有多少有訊號、多少沒訊號

> 「有訊號 vs 無訊號」**高度依賴定義**——因 `SampleASM_Sig`（tumor vs normal 整體甲基差）**飽和命中 29433（96.53%）**，把它算進去就近乎全部有訊號（無意義）。下表由嚴到寬列 5 種定義；**回答「甲基聚類結構」的問題請看定義 A**。

| 定義 | 有訊號 | 無訊號 | 說明 |
|---|--:|--:|---|
| A 結構訊號（聚類 或 within-HP）★推薦 | 24765（81.22%）| 5725（18.78%）| 純甲基聚類結構 |
| B 結構＋標籤對齊（對齊） | 22497（73.78%）| 7993（26.22%）| 結構且追隨標籤 |
| C 結構 或 扣SampleASM的Δβ | 27480（90.13%）| 3010（9.87%）| +germline/residual/subclone Δβ |
| D Significant（C++ 統一 verdict） | 28327（92.91%）| 2163（7.09%）| pipeline 保留判定 |
| E 寬（任一軸含 SampleASM）⚠ | 30144（98.87%）| 346（1.13%）| SampleASM 飽和→近全部 |

- **★ 推薦定義 A（純甲基聚類結構）**：**24765（81.22%）有結構訊號，5725（18.78%）完全無結構訊號**。
- 無結構的 5725 個**皆 ≥20 valid reads（非低覆蓋）= 真平坦**，且等同 §2 六格的「無群×單群」（與 partition 自洽）。
- ⚠ 定義 E（98.9% 有訊號）是 SampleASM 飽和灌出來的假象，**勿當「幾乎所有位點都有甲基訊號」引用**。

---

## §0.6 只用 tumor 資訊：有多少有訊號（嚴→寬）

> 只算 `!is_tumor continue` 程式碼證實的**純 tumor 軸**（WithinHP 結構 + Subclone/AltSubclone Δβ）；對齊/基礎聚類是 paired 故**不算**在此（見 §6）。4 層巢狀（T1⊆T2⊆T3⊆T4）。

| 層（嚴→寬）| 有 tumor 訊號 | 無 | 加入的軸 |
|---|--:|--:|---|
| T1 嚴格 · WithinHP PATTERN（silhouette 驗） | 354（1.16%）| 30136（98.84%）| tumor reads 在單一 HP 內形成 ≥2 乾淨甲基群 |
| T2 · WithinHP CleanMultigroup（PATTERN∪LEVEL） | 8077（26.49%）| 22413（73.51%）| ＋ mean-β 雙峰 level |
| T3 · ＋Subclone Δβ（germline vs carrier） | 11123（36.48%）| 19367（63.52%）| ＋ HP 內 germline vs somatic-carrier 甲基差 |
| T4 寬 · 任一 tumor 軸（＋AltSubclone） | 11764（38.58%）| 18726（61.42%）| ＋ within-family ALT/REF |

- **最嚴（T1 WithinHP PATTERN，silhouette 驗）只有 354（1.16%）**；**最寬（T4 任一 tumor 軸）＝ 11764（38.58%）**。
- 🔴 **對比 §0.5（paired）**：結構訊號 **24765（81.22%）→ tumor-only 最寬僅 38.58%**（約砍半，呼應 §6.4 Level-2 median V 0.665→0.319）→ 可偵測結構約半數來自 normal reads。
- 無 tumor 訊號的 18726 個中 **18606（佔無訊號 99.36%）是 tumor reads ≥20 足夠但真稀疏**（非低覆蓋）；僅 carrier-limited 414（1.36%）因低 VAF→somatic-ALT reads<3 使 subclone 軸測不了。
- 對拍：T4 ＝ 11764 ＝ tn_dependency T_only_signal_total（一致）。

---

## §0.7 ⭐ tumor 自身可偵測 vs 需 normal 歸因 — 互斥分割（差級，修正版）

> **先看邏輯（你問的重點）**：
> - **§0.5 / §0.6 = 包含關係（nested）**：T1⊆T2⊆T3⊆T4，每層含更嚴的層；數字**不可相加**。
> - **本節 §0.7 = 互斥分割（partition / 差級）**：每位點剛好落 ①②③④ **一格**；四格**加總 = 30,490**（可相加）。
>
> **軸定義（精確 — 修正前版誤標）**：
> - **T ＝ tumor 自身可偵測結構**（WithinHP / Subclone / AltSubclone Δβ，純 tumor reads）。⚠ 只表「tumor 看到結構」，**不代表能判定 somatic**——判 somatic 一定需要 normal。
> - **N ＝ 可歸因的 normal 軸**：germline-ASM（差異在 normal ＝ germline）**或** tumor−normal 殘差（tumor-specific）。**扣掉飽和的 SampleASM**。
> - ⚠ **N 內兩子型意義相反**：germline ＝「差異在 normal，故**非** tumor-specific」；residual ＝「**tumor-specific**，經相減確認」。故下方各格拆 germline / residual。

**2×2 互斥分割**

| | 無可歸因 normal 軸 | 有可歸因 normal 軸（N）| 列計 |
|---|--:|--:|--:|
| **T：tumor 可偵測結構 ✓** | ① 4698（15.4%）| ② 7066（23.2%）| **11764（38.6%）** |
| **無 tumor 結構 ✗** | ④ 9544（31.3%）| ③ 9182（30.1%）| 18726（61.4%）|
| **欄計** | 14242（46.7%）| **16248（53.3%）** | 30,490 |

**四格意義（修正）：**
- **① tumor 結構「未歸因」— 4698（15.4%）**：tumor 偵測到結構，但 germline-ASM 與殘差軸**皆不顯著** → 無法判 germline 還是 somatic（未定）。⚠ **不是「normal 不增訊息」**——這些位點 **97% 仍有 SampleASM** 整體 tumor-vs-normal 差異，只是不具歸因力。（＝ tn T2a）
- **② tumor 結構 ＋ normal 可歸因 — 7066（23.2%）**：拆 → **residual-only 3987（tumor-specific）／ germline-only 1431（差異在 normal）／ both 1648**。
- **③ 僅 normal 軸有 — 9182（30.1%）**：tumor 自身無結構，差異只在 normal 軸現形 → **residual-only 4528 ／ germline-only 2724 ／ both 1930**。
- **④ 無結構/無歸因軸 — 9544（31.3%）**：⚠ **不是「皆無」**——**94% 仍有 SampleASM** 整體差異，只是不具結構/歸因力。

**邊際與整體歸因：**
- **tumor 自身可偵測結構 T = 11764（38.6%）**；**可歸因 normal 軸 N = 16248（53.3%）**；兩者**重疊 = ② 7066**（故不可直接相加）。
- 整體歸因（重疊計）：**germline（差異在 normal）7733（25.4%）／ tumor-specific 殘差 12093（39.7%）**。

> ⚠ **SampleASM（tumor vs normal 整體差）在四格都命中 94–98%** → 無區辨力，故排除在 N 外。若算進去「需 normal」會膨脹到 97.6%（假象，**勿引用**）。

---

## §1 方法 — 甲基聚類怎麼算、怎麼驗（8 問）

<details><summary><b>Q1 如何算出甲基聚類？</b></summary>

- 每個 somatic 位點取 ±5000bp 窗口內的 reads，建 **read × CpG 甲基矩陣**（raw 機率 = ML/255）。
- 兩 read 距離 = **BERNOULLI**：`delta(p,q)=p(1−q)+(1−p)q`，每 CpG 權重 `w(p)=2|p−0.5|`（極端甲基權重高），`Dist=Σ(w·delta)/Σw`，C_min=3，nan→SKIP。
- 距離矩陣 → **UPGMA average-linkage** 階層聚類（`linkage_matrix.csv`）。**完全無監督**（標籤不進入）。
</details>

<details><summary><b>Q2 如何決定分幾群？</b></summary>

- `find_optimal_clusters` 掃 k∈[2, min(6, n/2)]，取 **最大 silhouette** 的 k（`optimal_k`）。
- 若分群不平衡（有群 < max(3,n/20)），且 k<5 → 嘗試 k+1（unbalanced bump，RegionProcessor.cpp:2475-2490）。
</details>

<details><summary><b>Q3 如何驗證與標籤一致性？（這就是「對齊」）</b></summary>

- 建 **cluster × a-priori 標籤 列聯表**（rows=無監督 cluster，cols=標籤），三軸各一張：allele(ALT/REF)、HP-family(HP1/HP2)、HP-fine(HP1/HP1-1/HP2/HP2-1)。
- 每張表算 **chi-square + Cramér's V（走 Cochran reliability gate，小期望數→壓 0）+ Fisher-Freeman-Halton**。
- `GlobalP = min(p_alt, p_hp, p_hp_family)`、`CramersV = max(可靠 V)` → 一個位點可**同時**在多軸對齊（**三軸非互斥**）。
- 全體 Fisher 顯著（any 軸）：22632（74.23%）。
- ⚠ HP-fine 對齊檢定**只用 4 群**（HP1/HP1-1/HP2/HP2-1）；HP3、unphase 只當組成計數，不進對齊檢定。
</details>

<details><summary><b>Q4 如何驗證「一群多標籤」（一結構→多標籤）？</b></summary>

聚類成功（PassedGating）**但** cluster×label 不顯著（GlobalP>0.05）= cluster 不追隨標籤、同一群混多種標籤 = **狀況③**（510，1.67%）。
</details>

<details><summary><b>Q5 如何驗證「多群一標籤」（多結構→一標籤）？</b></summary>

在**單一 germline HP 內**（tumor reads）再做甲基分群，兩種偵測：
- **PATTERN**：within-HP BERNOULLI UPGMA，silhouette≥0.5 + 平衡 → `WithinHP_NGroups≥2`。嚴格，全體 354。
- **LEVEL**：per-read 平均 β 1D 雙峰（gap>0.15 & varexpl>0.5）→ `WithinHP_LevelBimodal`。全體 7910。
- **CleanMultigroup = PATTERN ∪ LEVEL = 8077** = 狀況②。
</details>

<details><summary><b>Q6 如何只觀察 tumor 的數據？</b></summary>

C++ 中以 `if (!read_list[i].is_tumor) continue` 只取 tumor reads 的軸（**真 tumor-only**）：WithinHP 子結構（狀況②）、SubcloneDbeta、AltSubcloneDbeta、tumor-intrinsic 再聚類。註解明寫「independent of normal」（:1168）。
</details>

<details><summary><b>Q7 如何同時觀察 tumor 與 normal？</b></summary>

- **基礎聚類 + 對齊（狀況①③）本身就用 tumor+normal 混合 reads**（GlobalTest 無 is_tumor 過濾）。
- **needs-normal 軸**：GermlineAsmDbeta（normal HP1 vs HP2）、SomaticResidual / HP_Residual（tumor−normal 殘差）、SampleASM。
</details>

<details><summary><b>Q8 各分類的定義與「可能是」？</b></summary>

見 §3–§5 每狀況的 interpretation（窮舉可能 + 證據錨定）。原則：**推論可做，但列全可能；除非證據排除才不列；不可把單一可能寫成已確認**。
</details>

---

## §2 狀況盤點 — 乾淨 partition（6 格，Σ=30,490）

> 純度軸（對齊 / 一群多標籤 / 無群）× within-HP 軸（單群 / 一標籤多群）。**六格互斥且加總 = 30,490**。

| | 單群 | 一標籤多群 | 列小計 |
|---|--:|--:|--:|
| **對齊**（狀況①）| 16317（53.52%）| 6180（20.27%）| 22497 |
| **一群多標籤**（狀況③）| 371 | 139 | 510 |
| **無群**（覆蓋不足/未過 gating）| 5725 | 1758 | 7483 |

- **狀況①對齊** = 對齊兩格 = 22497；其中**最乾淨 1:1（對齊×單群）= 16317**。
- **狀況②一標籤多群** = 一標籤多群整欄（含對齊/一群多標籤/無群下的）= 8077。
- **狀況③一群多標籤** = 一群多標籤列 = 510。

---

## §3 狀況① 對齊 — 標籤細分（一結構→一標籤）

#### 狀況① 對齊（n=22497，佔全體 73.78%）

**軸 A — 以 somatic tag（6 類，主導 read 組成；對齊檢定僅用前 4 群）**

| tag | 位點數 | 總比例 /30490 | 區域比例 /本狀況 |
|---|--:|--:|--:|
| HP1 | 6059 | 19.87% | 26.93% |
| HP2 | 6083 | 19.95% | 27.04% |
| HP1-1 | 4134 | 13.56% | 18.38% |
| HP2-1 | 4063 | 13.33% | 18.06% |
| unphase | 1910 | 6.26% | 8.49% |
| HP3 | 248 | 0.81% | 1.1% |

> 帶 somatic carrier（HP1-1/HP2-1>0）：22222（98.78%）；純 germline：49（0.22%）

**軸 B — 以 phase（HP1-family / HP2-family，majority）**

| family | 位點數 | 區域比例 |
|---|--:|--:|
| HP1-family | 11130 | 49.47% |
| HP2-family | 10981 | 48.81% |
| balanced(=) | 202 | 0.9% |
| none | 184 | 0.82% |

**軸 C — 以 REF/ALT（近似 majority；真值見 Level-2）**

| allele | 位點數 | 區域比例 |
|---|--:|--:|
| REF-majority | 14498 | 64.44% |
| ALT-majority | 7762 | 34.5% |
| balanced(=) | 237 | 1.05% |
| none | 0 | 0.0% |

> **DominantLabel（主導軸 partition）**：hp=19153（85.14%） / allele=3344（14.86%）


**字面意義**：甲基距離分群沿 a-priori 標籤(HP/allele)分離 = 結構追隨標籤

**可能是（窮舉，未排除者全列）：**

- **germline ASM (甲基差異在 normal 已存在, 非腫瘤特異)** — 錨：`GermlineAsmDbeta_Sig=6295/22497 (27.98%)` ｜ 狀態：部分位點已標記, 非全部
- **cis-ASM (somatic 變異誘發的等位基因特異甲基)** — 錨：`需 normal-anchored cis-control 才能與 germline ASM 區分` ｜ 狀態：未排除, 待驗
- **tumor-specific somatic 甲基改變(單倍型上)** — 錨：`HP_Residual/SomaticResidual_Sig=9207/22497 (40.93%)` ｜ 狀態：部分位點已標記
- **true subclone (甲基與 haplotype/allele 標籤共分離)** — 錨：`需 multi-sSNV linkage + CCF (chr2:18M 法)` ｜ 狀態：未排除, 待驗
- **technical (覆蓋不均/strand bias 造成的表面相關)** — 錨：`pure_germline(無carrier)僅 49 (0.22%); DominantLabel=allele 3344 (14.86%)` ｜ 狀態：低權重但未排除

> ↓ 降權／提醒：未定(normal軸皆未顯著)=9836 (43.72%) → 約 4 成連 germline/somatic 都分不開, 不可貿然歸任何一類


---

## §4 狀況② 一標籤多群 — 標籤細分（多結構→一標籤；🟢 tumor-only 軸）

#### 狀況② 一標籤多群（n=8077，佔全體 26.49%）

**軸 A — 以 somatic tag（6 類，主導 read 組成；對齊檢定僅用前 4 群）**

| tag | 位點數 | 總比例 /30490 | 區域比例 /本狀況 |
|---|--:|--:|--:|
| HP1 | 2440 | 8.0% | 30.21% |
| HP2 | 2322 | 7.62% | 28.75% |
| HP1-1 | 1384 | 4.54% | 17.14% |
| HP2-1 | 1298 | 4.26% | 16.07% |
| unphase | 583 | 1.91% | 7.22% |
| HP3 | 50 | 0.16% | 0.62% |

> 帶 somatic carrier（HP1-1/HP2-1>0）：8073（99.95%）；純 germline：4（0.05%）

**軸 B — 以 phase（HP1-family / HP2-family，majority）**

| family | 位點數 | 區域比例 |
|---|--:|--:|
| HP1-family | 4095 | 50.7% |
| HP2-family | 3900 | 48.29% |
| balanced(=) | 82 | 1.02% |
| none | 0 | 0.0% |

**軸 C — 以 REF/ALT（近似 majority；真值見 Level-2）**

| allele | 位點數 | 區域比例 |
|---|--:|--:|
| REF-majority | 5612 | 69.48% |
| ALT-majority | 2365 | 29.28% |
| balanced(=) | 100 | 1.24% |
| none | 0 | 0.0% |

> **DominantLabel（主導軸 partition）**：hp=6666（82.53%） / allele=1411（17.47%）


**字面意義**：單一 HP 內甲基再分≥2 群 (subclone 最相關方向, 但須分 PATTERN/LEVEL)

**可能是（窮舉，未排除者全列）：**

- **true subclone (HP 內 somatic-carrier vs germline read 甲基不同)** — 錨：`Subclone/AltSubcloneDbeta_Sig=1786/8077 (22.11%); PATTERN(silhouette驗)=354` ｜ 狀態：PATTERN 子集較可信
- **epiallele / 隨機表觀多型(非克隆)** — 錨：`高甲基熵無克隆結構; 待 entropy/epipoly 檢` ｜ 狀態：未排除
- **LEVEL artifact (整體甲基高低位移, 非真 pattern 子結構)** — 錨：`LEVEL(mean-β雙峰)=7910/8077 (97.93%) 為大宗, 低信心` ｜ 狀態：大宗但低信心
- **cis-ASM 造成 HP 內等位甲基差** — 錨：`需 cis-control` ｜ 狀態：未排除
- **LOH 驅動的覆蓋偏移** — 錨：`Potential_LOH=505/8077 (6.25%)` ｜ 狀態：部分

> ↓ 降權／提醒：PATTERN(嚴格 silhouette 驗)僅 354; 其餘多為 LEVEL → 多群不等於真 subclone


---

## §5 狀況③ 一群多標籤 — 標籤細分（一結構→多標籤）

#### 狀況③ 一群多標籤（n=510，佔全體 1.67%）

**軸 A — 以 somatic tag（6 類，主導 read 組成；對齊檢定僅用前 4 群）**

| tag | 位點數 | 總比例 /30490 | 區域比例 /本狀況 |
|---|--:|--:|--:|
| HP1 | 172 | 0.56% | 33.73% |
| HP2 | 117 | 0.38% | 22.94% |
| HP1-1 | 63 | 0.21% | 12.35% |
| HP2-1 | 63 | 0.21% | 12.35% |
| unphase | 86 | 0.28% | 16.86% |
| HP3 | 9 | 0.03% | 1.76% |

> 帶 somatic carrier（HP1-1/HP2-1>0）：498（97.65%）；純 germline：3（0.59%）

**軸 B — 以 phase（HP1-family / HP2-family，majority）**

| family | 位點數 | 區域比例 |
|---|--:|--:|
| HP1-family | 262 | 51.37% |
| HP2-family | 233 | 45.69% |
| balanced(=) | 8 | 1.57% |
| none | 7 | 1.37% |

**軸 C — 以 REF/ALT（近似 majority；真值見 Level-2）**

| allele | 位點數 | 區域比例 |
|---|--:|--:|
| REF-majority | 361 | 70.78% |
| ALT-majority | 144 | 28.24% |
| balanced(=) | 5 | 0.98% |
| none | 0 | 0.0% |

> **DominantLabel（主導軸 partition）**：hp=404（79.22%） / allele=106（20.78%）


**字面意義**：甲基 cluster 不追隨標籤 (cluster⊥label, 同群混多標籤)

**可能是（窮舉，未排除者全列）：**

- **LOH (單一 HP 近消失, reads 無法依 HP 乾淨切)** — 錨：`Potential_LOH=93/510 (18.24%)` ｜ 狀態：部分
- **甲基結構由 HP/allele 以外的軸驅動(其他表觀/批次/strand)** — 錨：`需逐位點查 strand/批次` ｜ 狀態：未排除
- **覆蓋不均使某標籤太稀無法對齊** — 錨：`NReadsValid<15 =0/510 (0.0%)` ｜ 狀態：部分
- **true 生物 epiallele (與基因型獨立)** — 錨：`需 cis-control + entropy` ｜ 狀態：未排除

> ↓ 降權／提醒：狀況3 稀少(510, 1.67%); 先前歸因多為 LOH/覆蓋, 真乾淨 epiallele 罕見


---

## §6 只看 tumor vs 加上 normal 的差異

### 6.1 read-set provenance（哪些軸用哪些 read；C++ 源碼證實）

| 類別 | 軸 | 證據（行號）|
|---|---|---|
| **PAIRED** | 基礎聚類 / 對齊 GlobalP / DominantLabel | GlobalTest.cpp:149-304 無 is_tumor 過濾 |
| **TUMOR-ONLY** | WithinHP 子結構（狀況②）/ Subclone / AltSubclone Δβ | RegionProcessor.cpp:1172/1075/1089/1107 `!is_tumor continue` |
| **NORMAL-only** | GermlineAsmDbeta | :1064 只取 normal |
| **NEEDS BOTH** | SomaticResidual / HP_Residual / SampleASM | :1041 tumor−normal |

### 6.2 訊號分層（tn_dependency：tumor 自足 vs 需 normal）

| 層 | 意義 | 位點數 | 佔全體 |
|---|---|--:|--:|
| T1 | tumor 偵測 + normal 確認 somatic | 5635 | 18.48% |
| T2a | tumor 偵測但 germline/somatic 未分（需 normal）| 4698 | 15.41% |
| T2b | tumor 偵測但 normal 顯 germline | 1431 | 4.69% |
| T3 | 僅 normal 軸偵測 somatic（tumor 單看不足）| 6458 | 21.18% |
| T4 | 純 germline（normal）| 2724 | 8.93% |
| T5 | 無上述標記訊號 | 9544 | 31.3% |

> **純 tumor 自足（無任何 normal 軸）僅 149（0.49%）**；需 normal 才有的訊號 18156（59.55%）。
> HP-甲基判別力：normal 較強 14268（70.73%）> tumor 較強 5898（29.24%）→ HP 甲基**多是 germline**，tumor 更強才是 tumor-intrinsic 候選。

### 6.3 對齊位點加上 normal 的再歸因（22497 個對齊位點）

| 再歸因 | 位點數 | 佔對齊 |
|---|--:|--:|
| germline ASM only（normal 也有，非 somatic）| 3454 | 15.35% |
| somatic-specific only（扣 normal 仍顯著）| 6366 | 28.3% |
| both flagged | 2841 | 12.63% |
| neither（未定）| 9836 | 43.72% |

### 6.4 Level-2 read-level 校準（tumor-only 對齊是否成立）

限 tumor reads 重算 cluster×hp_fine 關聯（在 paired cluster 內）：可測 159 中 **82 個（51.6%）保持**（V≥0.3）；median V paired **0.665** → tumor-only **0.319**。
**→ 約半數對齊位點在 tumor reads 自身仍有中等訊號，但整體強度約半數來自 normal reads。**

---

## §7 誠實邊界（看數字必讀）

1. **「是訊號」≠「是 subclone」**：每個狀況只代表「通過了哪些驗證 → 字面意義」，**不等於**已確認 subclone 身分。生物身分需 §3–§5「可能是」逐一驗（cis-control / multi-sSNV linkage+CCF / CN 軸 / 跨樣本）。
2. **Level-1 是「歸屬」非「已驗證 cell 對應」**：per-site CSV 無 cluster×label 逐 cell 計數；標籤細分用 DominantLabel + 組成 argmax 歸屬。Level-2 抽樣校準：類別命中 92.6%、軸命中 73.7%、映射驗證 89.7%。
3. **三軸非互斥**：一位點可同時在 hp-fine / hp-family / allele 對齊；DominantLabel 是其中主導軸。
4. **HP-fine 對齊檢定僅 4 群**；HP3 / unphase 只組成計數。
5. **對齊基礎是 paired**（混 normal reads）；真正 tumor-only 的是狀況②與 subclone Δβ 軸。
6. **單一 paired run**（`_wg_d1_unified`）；跨樣本復現另計。

---
*方法 appendix（列聯表/方向/PERMANOVA double-dip 細節）→ `InterSubMod/docs/methodology/20260617_structure_label_association_explained_01.standalone.html`*
*數據來源 JSON：`_assets/20260616_readset_provenance/{label_breakdown,contingency_direction,tn_dependency,level2_calibration}.json`；驗證 `verify_breakdown.py` 25/25 PASS。*
