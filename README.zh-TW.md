# InterSubMod

**用 ONT 長讀長的單分子體細胞突變共現，重建局部突變狀態候選結構。**

[English →](README.md) · [**文件網站 →**](https://liaoyoyo.github.io/InterSubMod/) · [Wiki →](https://github.com/liaoyoyo/InterSubMod/wiki) · [怎麼跑 →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

腫瘤可能包含多個帶著不同突變組合的細胞群。理解突變順序與細胞共屬性很重要，
但**未加細胞條碼的 bulk 長讀長不會直接觀察到這兩個生物學量**。

若輸入只有各位點的邊際變異等位頻率，沒有 linkage 或額外模型假設，可能有多個聯合結構
產生相同邊際值，因此無法只從這些邊際值識別唯一聯合結構。ONT 長讀長增加了另一種觀測：
一條分子可以同時跨過好幾個體細胞突變，因此「它們是否共現在同一條物理分子」
變成**直接觀測**；但因為來源細胞未知，細胞共屬性與譜系仍是 model-dependent 推論。

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
| 兩支引擎 | **`inter_sub_mod`**（C++17）· `longlineage`（C++17） | per-region 甲基化／統計 · 版本限定的 per-read artefact |
| 呈現層 | Python 分析 ＋ standalone HTML | 圖表、漏斗、互動判讀工作站 |

---

## 核心方法論斷言

這個專案最重要的設計決定，是**什麼東西有資格驅動重建**——而甲基化被刻意排除在外。

![為什麼甲基化不能驅動重建](docs/images/methylation-circularity.png)

當你在某個位點看到兩群甲基化模式不同的 read，至少有四種可能成因：
germline 等位特異性甲基化、雜合性缺失造成的解遮蔽、拷貝數劑量效應、以及真正的譜系差異。
**以目前 single-bulk measurement set、且沒有 orthogonal data 或額外假設時，無法識別這四者。**
若要用甲基化獨立「確認」細胞亞群，必須另有外部細胞歸屬；依突變定義群組後觀察到
一致甲基化，只是 concordant association，不是獨立確認。

所以骨幹是**同一條物理分子上的體細胞突變共現**——它不依賴任何待推論的標籤，
因此**非循環**。甲基化被保留為嚴格的 **bounded-auxiliary（有界輔助）**訊號：
它在遺傳候選結構固定**之後**才計算，只做註記，**動不了任何一條邊**。

> 實測數據支持這個克制：811 個可評估的甲基化單元中，
> 只有 **3 個（0.37%）**支持穩健的 pattern-conditioned association。這個低 yield 支持
> 不用甲基化挑拓撲，但不是對所有可能甲基訊號的總檢定。

---

## 實際跑出來的結果

![全 7 樣本 funnel](docs/images/funnel-7samples.png)

**7 個資料集（6 個生物樣本）、chr1–22** 的 canonical 結果，2026-08-01 凍結。
每一層都可稽核，且算術自洽
（39,648 + 23,858 = 63,506；+ 8,449 = 71,955；+ 3,224 + 45 = 75,224；+ 10,717 = 85,941；+ 13,014 = 98,955）。

| 項目 | 數值 |
|---|---|
| sSNV 資料列 | 469,849 |
| `k=1` strict read-linkage components | 170,131 / 255,752 strict components（66.52%） |
| 帶突變的分析單元 | 85,941 |
| 因搜尋節點上限而主動放棄 | 10,717（12.47%） |
| 可用 read-AF 排序的單元 | 71,955 |
| 收斂到單一 **rooted-unlabeled 數學拓撲 signature** | 63,506（可排序者的 **88.26%**） |
| **確認的細胞亞群** | **0** — 見下 |

> **88.26% 的正確讀法。**它的意思是：*在已經可排序的 71,955 個單元中*，
> 有 88.26% 在 frozen recurrence-allowed model 下收斂到單一數學形狀。這是一個
> **model-conditional 圖形統計**，**不是**「腫瘤演化史已解出 88%」。另一個不同 grain
> 的數字是 170,131 / 255,752 strict components（66.52%）為 `k=1`；相對 469,849 筆
> sSNV dataset records，170,131 是 36.21%。這些分母不能互換。

---

## 能力邊界

這個專案輸出**機器可讀的宣稱邊界**。canonical `cohort_receipt.json` 與
`summary/all7_summary.json`（不是 `authority_manifest.json` 頂層）明白寫著
`technical_all_pass = true` 但 `validation_evidence_eligible = false`：
所有雜湊都對、所有測試都過，然而系統自己宣告這批結果**還不能作為驗證證據**。

<table>
<tr><th width="50%">✔ 可以主張</th><th width="50%">✘ 明文禁止主張</th></tr>
<tr valign="top"><td>

- 嚴格 read-linked 的**局部**結構
- 家族完整時的完整最小候選家族
- 局部、允許 recurrence 的 Hamming-1 candidate arborescences
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

下表是**版本限定**的驗證狀態；「已驗證」只代表所列 artifact 在指定版本實跑或核對，
不是對每一條公開指令與每一個檔案的 blanket guarantee。

| 部件 | 狀態 | 說明 |
|---|:---:|---|
| `inter_sub_mod` | ✅ 可跑 | tracked core `73afaeac` 的 fresh build/run audit；命令收據見下方稽核 |
| C++ 測試套件 | ✅ 全過 | 2026-08-12 fresh run：**270 tests / 39 suites**，CTest **270/270** |
| 7 樣本 HP/PS sidecar | ✅ 齊全 | 7/7 PASS |
| research exact-PS topology solver | ✅ 可跑 | 產生上方 exact-PS funnel 的獨立 research solver，不是 `inter_sub_mod` |
| `longlineage preflight` | ✅ 可跑 | 驗證 8 角色 manifest |
| `longlineage dataset-gate` | ⚠️ 受限 | 唯一能出科學結果的入口，但**硬編碼綁死單一資料集** |
| `longlineage run` / `probe` | 🔴 被鎖 | 依設計回 `KernelBlocked` |
| `longlineage` 拓撲產出 | 🔴 0 個單元 | 見[下方說明](#關於-longlineage) |
| 輸出帶標籤的 BAM | ✅ 可跑 | `inter_sub_mod` 不寫 BAM；LongLineage public main `583e03e` 已建出 `longlineage-tag-bam`（`CMakeLists.txt:221`），2026-08-17 併入 |
| 串起兩支引擎的單一腳本 | ⚠️ 只有單引擎的 | LongLineage 有 `scripts/run_sample.sh`（partition → tagged BAM）。尚無單一腳本同時串起 LongLineage **與** `inter_sub_mod`；執行順序記載於 [LongLineage README](https://github.com/liaoyoyo/LongLineage#與-intersubmod-的關係先讀這節再決定要不要裝) |

### 關於 LongLineage

`longlineage` 是同一套科學想法的 **clean-room 工業化重寫**
（**git 根 commit 完全不同**——不是 fork）。它的契約設計相當好：
每個 artefact 都被 schema 鎖死並附 SHA-256 收據，
`topology_unit` 更把「解到什麼程度」拆成四個獨立狀態欄位，每種放棄都有具名理由。

在 frozen **HCC1395 dataset-gate receipt** 中它輸出 0 個 topology units；
此結果不能外推到所有 LongLineage real-data run：

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

# 2. 驗證編譯結果   -> 2026-08-12 稽核為 270 tests / 39 suites，全數通過
./build/bin/run_tests

# 3. Python 依賴
pip install -r requirements.txt

# 4. 先用完全合成的 fixture 跑一次，確認建置能端到端運作。
#    需要 pysam + samtools + bgzip + tabix，產生的輸入約 14 KB。
#    2026-08-17 實測：exit 0、2 個 region、244 個 CpG、每 region 24 條 read。
python3 scripts/make_synthetic_fixture.py
bash tests/fixtures/synthetic/RUN.sh

# 5. 接著才提供你自己有授權且已建立索引的輸入。repo 不附**真實**的
#    BAM/FASTA/VCF；以下仍是 placeholder，不是可直接複製的資料路徑。
./build/bin/inter_sub_mod \
  --tumor-bam /path/to/tumor.mm_ml.bam \
  --reference /path/to/reference.fa \
  --vcf       /path/to/candidates.vcf \
  --output-dir out_min
```

> 合成 fixture **沒有任何生物學意義**。它唯一的用途是給你一個已知會成功的基準，
> 這樣當你自己的資料跑不出東西時，才能區分「我的資料沒有訊號」與「我的環境壞了」。
> 見 [`tests/fixtures/synthetic/README.md`](tests/fixtures/synthetic/README.md) ——
> 裡面也記錄了一個靜默失敗的坑：`inter_sub_mod` 只認 MM 標籤的 `C+m?` flavor，
> 其他 flavor 會得到 `Total CpG sites found: 0` 而 **exit code 仍是 0**。

### 和 LongLineage 一起跑

兩個引擎各自獨立建置、CMake 互不引用，但執行期以**單一方向**銜接：

```
上游（不屬於這兩個 repo）：dorado（MM/ML）→ 比對 → LongPhase haplotag → somatic VCF
      ↓
LongLineage   scripts/run_sample.sh   →  寫入 BAM aux tag lc:Z lu:Z lv:Z ls:A
      ↓
InterSubMod   inter_sub_mod --tumor-bam <tagged.bam> ...
```

**InterSubMod 沒有 LongLineage 也完全跑得起來。** 輸入 BAM 若沒經過 `longlineage-tag-bam`，
lineage 軸會被直接跳過，而不是把所有 read 併成一個大群組
（見 `include/core/Stats.hpp` 的 lineage 軸註解）。反方向 —— LongLineage 在
`src/compat/regional_crosswalk.cpp` 讀取 InterSubMod 的 manifest —— 是**驗證用的
crosswalk，不是相依**，所以不要以為必須先跑 InterSubMod。

搭檔 repository 與完整執行順序：
[**LongLineage**](https://github.com/liaoyoyo/LongLineage#與-intersubmod-的關係先讀這節再決定要不要裝)。
確認本 README 記載的版本是否仍等於 LongLineage 公開 `main`：

```bash
python3 scripts/handoff/check_companion_version.py
```

<details>
<summary><b>三個一定會咬到你的地方</b>（皆已對照原始碼確認）</summary>

- **`--threads` 預設是 16 不是 1。** help 寫 `Default: 1`，但 `Config.hpp` 的實際預設是 `16`。
  資源估算會整整差 16 倍。
- **`--distance-metric` 會被靜默改成 `NHD`。** `Config.hpp` 宣告 `BERNOULLI`，
  但參數解析器會無條件清空再填入 NHD。**永遠明確指定它。**
- **`tree.nwk` 的葉子是 *read* 不是 clone。** 它是「read 依甲基化相似度」的階層分群樹，
  **不是**亞群演化譜系樹。這是最常被誤讀的輸出。

另外兩個版本敏感點：`methylation.csv` 第一欄是矩陣列號而非 read 名
（與 read 的綁定完全靠列序，沒有 key 校驗）。在 audited core `73afaeac`，
`significance_summary.csv` 有 **199 欄**，包含 `VerificationSchemaVersion=2` 與
`RegionStratificationSchemaVersion=1`；它們是 component schema 欄，不是全檔 layout version。
請用**欄名**取值並釘定 producing commit。

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

> **同一份內容有兩種讀法。**上表連到 **Wiki**，GitHub 原生渲染、最適合快速瀏覽。
> 同樣的內容也以完全自足的 HTML 形式發佈在 **[liaoyoyo.github.io/InterSubMod](https://liaoyoyo.github.io/InterSubMod/)**。
> `docs/explain/` 是 editorial upstream；Wiki 是人工同步的衍生版本，且發布是獨立步驟。
> 在稽核的 Pages deploy `fbdf7c7`，17 個 standalone 頁面含 **37 個 inline `<svg>` elements**；
> 計數命令見 correction receipt。這個版本限定的 element count 不等於 37 張語意獨立圖，
> 也不表示 Wiki 與 Pages byte-identical。

第 01–10 頁涵蓋科學方法本身（名詞地基、ISM 方法、甲基化讀取與篩選、
真實案例逐步拆解、三統計分工、能力 vs 敘述）。

---

## 兩個值得借用的設計

這個 repo 裡有兩個模式可以推廣到基因體學以外。

**1 · 串流取代落地。** 標記後的串流可導進具名管道，即時抽成目前的 9 欄 sidecar contract。
七個經稽核的 sidecars 合計 **6,256,168,164 bytes（5.83 GiB）**。先前顯示的
**1.67 TiB** tagged-BAM total 沒有已提交的七檔 path／exact bytes／hash／compression receipt，
因此標為 **UNVERIFIED**，也不宣稱任何縮減倍率。目前 sidecar 不保留 `SEQ`／`QUAL`；
它們的 byte share 與 downstream utility 尚未做 field-level census。

![上游工具鏈](docs/images/upstream-toolchain.png)

**2 · 對缺少的必填指標 fail closed。** 經稽核的工作站生成器吃一份宣告式 spec，
當已宣告的必填指標缺失時**直接拒絕渲染（exit 3）**。這可防止那些已宣告欄位被靜默略過，
但不會驗證內容真偽、分母或來源，也不能偵測未宣告的 optional fields。

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

## 本文件的驗證範圍

本 README 受
[2026-08-12 公開 claim 稽核](docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md)
約束。Frozen exact-PS 數字已對 authority manifest 與 denominator registry；tracked C++ core
也 fresh build/test。可是公開 quick-start 資料**沒有隨 repo 提供**，GitHub About／Wiki／Pages
有各自發布狀態，LongLineage capability 也只適用上面釘定的 commit，因此本文件不宣稱
「每個數字與指令都已驗證」。圖片由 `tools/extract_svg_for_github.py` 從 `docs/explain/`
抽取；應修改來源頁後重生，而不是直接編圖片。

| artefact／claim 家族 | 驗證身分與日期 | 可重跑檢查與結果 | 適用範圍與已知失敗 |
|---|---|---|---|
| Frozen exact-PS funnel | Frozen authority artefacts，2026-08-12 重查 | manifest／hash census 加獨立分母重算；精確計數皆重現 | 只限 frozen 7-dataset analysis；不是 `inter_sub_mod` CLI，也不能識別 cellular clones |
| Tracked C++ core | `73afaeac-dirty`、GCC 11.4.0、htslib 1.18；2026-08-12 執行 | `cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release`、build、直接 GoogleTest 與 CTest：build exit 0；270 tests / 39 suites；CTest 270/270 | 經稽核的 C++／CMake build inputs 與 remote feature `ddd8909` byte-equivalent；實跑身分仍是 local-dirty，不是 clean-checkout release certification |
| 公開 quick start | Current source，2026-08-13 檢查 | Build/test 路徑可重跑；分析命令明示由使用者提供輸入 placeholder | Repo 未附公開 tumor BAM／reference／VCF fixture，因此不宣稱 end-to-end biological result 已重現 |
| LongLineage 狀態 | Public main `583e03e`（2026-08-17）、frozen HCC1395 evidence | 三條 `agent/public-preview-*` 分支與 tag-bam／solver commit 併入公開 `main`，於隔離 worktree 重建並重跑測試：build exit 0、**ctest 49/49** | `scripts/ci/check_public_preview_gate.sh HEAD` 仍回報 **FAIL，5 個未結案 blocker**（授權審查未結案、4 筆 `NO_GO` 來源、21 筆授權未核准、11 個相依項 `NOASSERTION`、歷史中 4 個 blob 含開發機絕對路徑）。可讀可建置，**不是**已完成授權清算的釋出版。尚無 7 樣本 runtime／memory 或多資料集 topology 驗證 |
| GitHub surfaces | 2026-08-13 本地來源已校正；GitHub About 為 `RESOLVED_LIVE` | Live API 重查，加上本輪 P0/P1/P2 claim guards 與 source checks | Default `main`、Wiki 與 Pages 仍是先前部署 bytes，尚待發布，不得描述成 live 已校正 |

完整命令與實際輸出保存在
[command receipts](research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md)；
逐 claim 修正狀態記錄於
[P0 correction cycle](research/20260813_public_docs_p0_correction/00_INDEX.md)；
[2026-08-13 公開介面更新循環](research/20260813_intersubmod_public_surfaces_refresh/00_INDEX.md)
則以[遠端狀態收據](research/20260813_intersubmod_public_surfaces_refresh/remote_state_receipt.md)
明確分開本地修正與 live 發布。

**誠實標註的已知缺口**：拷貝數目前是 `NOT_INTEGRATED`；
LongLineage 的 7 樣本執行時間與記憶體上界**從未實測**，
且其自身文件明文禁止由單一樣本外推。
