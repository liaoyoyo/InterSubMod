# InterSubMod

**用 ONT 長讀長的「單分子體細胞突變共現」重建腫瘤亞群譜系。**

[English →](README.md) · [解釋中心 →](https://github.com/liaoyoyo/InterSubMod/wiki) · [怎麼跑 →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

一顆腫瘤不是單一細胞群，而是好幾群帶著不同突變組合的細胞混在一起。
要理解抗藥性與病程發展，就需要知道**哪個突變先發生**、**哪些突變住在同一群細胞裡**。

短讀長測序只看得到每個位點各自的變異等位頻率，而要從這些邊際頻率反推聯合結構，
在數學上是**已被證明無解**的反卷積問題。ONT 長讀長改變了問題本身：
一條分子可以同時跨過好幾個體細胞突變，於是
「這兩個突變在不在同一個細胞譜系裡」從**頻率推論**變成了**直接觀測**。

這個 repo 把這個想法從頭到尾實作出來——並且同樣重要的是，
它對「證據到哪裡就用完了」講得很清楚。

---

## 系統架構

![系統架構](docs/images/architecture-overview.png)

由下而上五層。**只有一部分現在真的跑得動**——見下方[誠實狀態表](#誠實狀態表)。

| 層 | 部件 | 做什麼 |
|---|---|---|
| 原始資料 | tumour / normal BAM、參考基因組 | 帶甲基化標記的 ONT 長讀長 |
| 上游工具 | Dorado · ClairS · LongPhase-S · SAVANA | 鹼基判讀＋甲基化、體細胞突變、單倍型標記、拷貝數 |
| 中介契約 | **HP/PS sidecar TSV** | 每條 read 一列；兩支引擎唯一 byte-level 一致的介面 |
| 兩支引擎 | **`inter_sub_mod`**（C++17）· `longlineage`（C++17） | per-region 統計 · per-read artefact |
| 呈現層 | Python 分析 ＋ standalone HTML | 圖表、漏斗、互動判讀工作站 |

---

## 核心方法論斷言

這個專案最重要的設計決定，是**什麼東西有資格驅動重建**——而甲基化被刻意排除在外。

![為什麼甲基化不能驅動重建](docs/images/methylation-circularity.png)

當你在某個位點看到兩群甲基化模式不同的 read，至少有四種可能成因：
germline 等位特異性甲基化、雜合性缺失造成的解遮蔽、拷貝數劑量效應、以及真正的譜系差異。
**單一 bulk 樣本無法區分這四者。**
因此用甲基化去「確認」某群 read 是亞群，前提是你已經知道亞群的歸屬——
而那正是待證明的結論本身。

所以骨幹是**同一條物理分子上的體細胞突變共現**——它不依賴任何待推論的標籤，
因此**非循環**。甲基化被保留為嚴格的 **bounded-auxiliary（有界輔助）**訊號：
它在樹已由遺傳證據定好**之後**才計算，只做註記，**動不了任何一條邊**。

> 實測數據支持這個克制：811 個可評估的甲基化單元中，
> 只有 **3 個（0.37%）**達到穩健關聯。若當初拿它當骨幹，會發現它幾乎沒有訊號可用。

---

## 實際跑出來的結果

![全 7 樣本 funnel](docs/images/funnel-7samples.png)

**7 個資料集（6 個生物樣本）、chr1–22** 的 canonical 結果，2026-08-01 凍結。
每一層都可稽核，且算術自洽
（39,648 + 23,858 = 63,506；+ 8,449 = 71,955；+ 3,224 + 45 = 75,224；+ 10,717 = 85,941；+ 13,014 = 98,955）。

| 項目 | 數值 |
|---|---|
| sSNV 資料列 | 469,849 |
| 單點無共現夥伴（無法建樹） | 170,131（66.52%） |
| 帶突變的分析單元 | 85,941 |
| 因搜尋節點上限而主動放棄 | 10,717（12.47%） |
| 可用 read-AF 排序的單元 | 71,955 |
| 收斂到**單一 rooted-unlabeled 拓撲** | 63,506（可排序者的 **88.26%**） |
| **確認的細胞亞群** | **0** — 見下 |

> **88.26% 的正確讀法。**它的意思是：*在已經可排序的 71,955 個單元中*，
> 有 88.26% 收斂到單一樹形。這是一個 **model-conditional 的圖形統計**，
> **不是**「腫瘤演化史已經解出 88%」——全部突變中有三分之二早在上游就以孤立單點流失了。

---

## 能力邊界

這個專案輸出**機器可讀的宣稱邊界**。canonical 結果的欄位裡明白寫著
`technical_all_pass = true` 但 `validation_evidence_eligible = false`：
所有雜湊都對、所有測試都過，然而系統自己宣告這批結果**還不能作為驗證證據**。

<table>
<tr><th width="50%">✔ 可以主張</th><th width="50%">✘ 明文禁止主張</th></tr>
<tr valign="top"><td>

- 嚴格 read-linked 的**局部**結構
- 家族完整時的完整最小候選家族
- 允許遞迴的 Hamming-1 父子樹候選
- **未經 CN/LOH 校正**的 deterministic read-AF 排序
- exact rooted-unlabeled 拓撲普查
- pattern-conditioned 的甲基化**關聯**

</td><td>

- **確認的細胞層級 clone / subclone**
- CN/LOH 校正過的 CCF
- 校準過的 likelihood 或 posterior
- 把 read cluster 當成 cell clone
- 用甲基化差異重排拓撲
- 把 basecaller 重跑當成生物學重複

</td></tr>
</table>

**單一 bulk 樣本只能 characterize（描述），不能 confirm（確認）。**
要突破這個天花板需要 single-cell 或 multi-region 資料——
這是**資訊論上的界限，不是實作缺口**。推論出的內部節點一律標為 `inferred`。

---

## 誠實狀態表

以下每一項都是**實際跑過**確認的，不是照文件宣稱抄的。

| 部件 | 狀態 | 說明 |
|---|:---:|---|
| `inter_sub_mod` | ✅ 可跑 | 最小指令 **2.9 秒**完成、exit 0 |
| C++ 測試套件 | ✅ 全過 | **265 tests / 38 suites**，2.06 秒 |
| 7 樣本 HP/PS sidecar | ✅ 齊全 | 7/7 PASS |
| 分層樹枚舉 solver | ✅ 可跑 | **這才是產出上面那些數字的路徑** |
| `longlineage preflight` | ✅ 可跑 | 驗證 8 角色 manifest |
| `longlineage dataset-gate` | ⚠️ 受限 | 唯一能出科學結果的入口，但**硬編碼綁死單一資料集** |
| `longlineage run` / `probe` | 🔴 被鎖 | 依設計回 `KernelBlocked` |
| `longlineage` 拓撲產出 | 🔴 0 個單元 | 見[下方說明](#關於-longlineage) |
| 輸出帶標籤的 BAM | 🔴 不支援 | 兩支引擎都不能寫 BAM；LongLineage 更在契約層禁止 |
| 串起兩支引擎的單一腳本 | 🔴 沒有 | 目前是兩條各自獨立的線 |

### 關於 LongLineage

`longlineage` 是同一套科學想法的 **clean-room 工業化重寫**
（**git 根 commit 完全不同**——不是 fork）。它的契約設計相當好：
每個 artefact 都被 schema 鎖死並附 SHA-256 收據，
`topology_unit` 更把「解到什麼程度」拆成四個獨立狀態欄位，每種放棄都有具名理由。

但**在真實資料上它目前輸出 0 個拓撲單元**，而這**不是 bug**：

![LongLineage 漏斗](docs/images/longlineage-funnel.png)

它的流程把拓撲建構**閘控在甲基化分群的穩定性上**，
於是 79,687 個位點中有 66,836 個（83.9%）根本不會產生任何一列共現資料。
這與上面的方法論立場**直接衝突**，任何跨引擎比較都必須明講這一點。

這裡的「BLOCKED」指的是**對照驗證的證據尚未存在**——**不是**程式碼沒寫。
M1／M2／topology 的核心都已實作，也被實際執行過。

---

## 快速開始

每一步都附預期輸出的完整流程：[**怎麼跑**](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

![六個步驟](docs/images/howto-six-steps.png)

```bash
# 1. 編譯（repo 內的執行檔是 STALE 的，務必重新編譯）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 2. 驗證編譯結果   -> 應印出 "265 tests from 38 test suites ... PASSED"
./build/bin/run_tests

# 3. Python 依賴
pip install -r requirements.txt

# 4. 單一位點試跑（約 3 秒）
./build/bin/inter_sub_mod \
  --tumor-bam data/bam/HCC1395/tumor.bam \
  --reference data/ref/hg38.fa \
  --vcf       one_snv.vcf \
  --output-dir out_min
```

<details>
<summary><b>三個一定會咬到你的地方</b>（皆已對照原始碼確認）</summary>

- **`--threads` 預設是 16 不是 1。** help 寫 `Default: 1`，但 `Config.hpp` 的實際預設是 `16`。
  資源估算會整整差 16 倍。
- **`--distance-metric` 會被靜默改成 `NHD`。** `Config.hpp` 宣告 `BERNOULLI`，
  但參數解析器會無條件清空再填入 NHD。**永遠明確指定它。**
- **`tree.nwk` 的葉子是 *read* 不是 clone。** 它是「read 依甲基化相似度」的階層分群樹，
  **不是**亞群演化譜系樹。這是最常被誤讀的輸出。

另外兩個會**靜默出錯**的：`methylation.csv` 第一欄是矩陣列號而非 read 名
（與 read 的綁定完全靠列序，沒有任何 key 校驗）；
`significance_summary.csv` 的**欄數會隨 binary 版本改變**且檔案裡沒有版本欄位——
一律用**欄名**取值，不要用欄號。

</details>

---

## 文件導引

**解釋中心**（`docs/explain/`）是主入口——17 個 standalone 頁面，零外部依賴，可離線開啟。

| | 頁面 | 給誰看 |
|---|---|---|
| 🗺️ | [系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) | **先看這頁**——架構、誠實狀態、funnel、能力邊界 |
| ⚙️ | [InterSubMod 部件](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) | 3 個輸入、8 個階段、17 種輸出檔（附真實 header）、9 條陷阱 |
| 🧬 | [LongLineage 部件](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) | 子命令狀態、M1→M2→topology 鏈、artefact 契約 |
| 🔬 | [上游與資料](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) | Dorado / ClairS / LongPhase-S / SAVANA、sidecar 格式、7 樣本 |
| 📊 | [分析與呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) | 該用哪些腳本、拒絕渲染的 HTML 生成器 |
| ▶️ | [怎麼跑](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run) | 六個步驟，每步附驗收條件 |

> **關於上面的連結。**這些連結指向專案 **Wiki**，GitHub 會原生渲染。
> 同樣的內容也以完全自足的 HTML 頁面存在於 `docs/explain/`——排版更豐富、有可展開的細節，
> 但 GitHub 對 `.html` 只會顯示原始碼，所以那些請 clone 後開啟（或在 `/docs` 啟用 GitHub Pages）。

第 01–10 頁涵蓋科學方法本身（名詞地基、ISM 方法、甲基化讀取與篩選、
真實案例逐步拆解、三統計分工、能力 vs 敘述）。

---

## 兩個值得借用的設計

這個 repo 裡有兩個模式可以推廣到基因體學以外。

**1 · 串流取代落地。** 7 個樣本的 haplotagged BAM 合計 **1.67 TiB**。
做法改成把標記後的串流導進具名管道，即時抽成 9 欄的 sidecar——
**5.83 GiB，縮小 287 倍**——因為分析真正需要的只有
「哪條 read、在哪裡、屬於哪個單倍型」。
序列與品質字串佔了 99% 以上的體積，卻是 0% 的用途。

![上游工具鏈](docs/images/upstream-toolchain.png)

**2 · 讓捏造在結構上不可能。** 報告生成器吃一份宣告式 spec，
當必填指標缺失時**直接拒絕渲染（exit 3）**——它不會填破折號、也不留空白。
缺少的數字沒辦法被悄悄裝扮成存在的數字。

![拒絕渲染設計](docs/images/workstation-refuse-design.png)

---

## 環境需求

- **C++17** 編譯器、CMake ≥ 3.14、OpenMP
- **htslib**（BAM/VCF I/O）、Eigen
- **Python ≥ 3.9**（`requirements.txt`）；注意少數嚴格切割腳本需要 **3.10**
- 參考 FASTA **必須有 `.fai`**（缺索引是硬性錯誤）
- 輸入 BAM 必須帶 `MM`／`ML` 甲基化標記（沒有的 read 會被丟棄）

## 目錄結構

```
include/core/ src/core/     C++17 分析引擎
scripts/                    Python 分析入口
docs/explain/               解釋中心（17 個 standalone 頁面）
docs/images/                本 README 與 wiki 使用的圖
tools/                      渲染、QA 與抽取工具
state/                      研究 cycle 狀態機
```

## 本文件的狀態

本 README 中的每個數字、指令與檔案格式，都在 **2026-08-06** 經由實際執行指令與閱讀原始碼驗證。
圖片由 `tools/extract_svg_for_github.py` 從 `docs/explain/` 產生——
要改圖請改上游解釋頁後重跑，不要直接編輯圖片。

**誠實標註的已知缺口**：拷貝數目前是 `NOT_INTEGRATED`；
LongLineage 的 7 樣本執行時間與記憶體上界**從未實測**，
且其自身文件明文禁止由單一樣本外推。
