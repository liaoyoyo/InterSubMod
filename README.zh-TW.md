# InterSubMod

**研究交接第一入口：**[2026-08-13 完整研究資料與軟體交接 snapshot](docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md)——負責導航與治理；科學數值 authority 仍是 2026-08-01 frozen bundle，發布／release gates 仍 fail-closed。

**用 ONT 長讀長的單分子體細胞突變共現，重建局部突變狀態候選結構。**

[English →](README.md) · [**文件網站 →**](https://liaoyoyo.github.io/InterSubMod/) · [Wiki →](https://github.com/liaoyoyo/InterSubMod/wiki) · [怎麼跑 →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

腫瘤可能包含多個帶著不同突變組合的細胞群。理解突變順序與細胞共屬性很重要，
但**未加細胞條碼的 bulk 長讀長不會直接觀察到這兩個生物學量**。

短讀長測序只看得到每個位點各自的變異等位頻率；若**只靠 per-locus 邊際 VAF，
且沒有 linkage 或額外模型假設**，聯合結構在此逆問題下不可識別。ONT 長讀長改變了問題本身：
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
| read 標籤表示 | **BAM HP/PS aux tags** · **HP/PS sidecar TSV** | `inter_sub_mod` 從 BAM aux tags 讀 HP/PS；exact-PS／LongLineage 使用 sidecar。兩者是 commit-scoped 的平行契約，不是兩支引擎直接串接的執行期介面。 |
| 兩支引擎 | **`inter_sub_mod`**（C++17）· `longlineage`（C++17） | per-region 甲基化／統計 · 版本限定的 per-read artefact |
| 研究／呈現層 | commit-pinned Python research solver · validator/builder · standalone HTML | Python solver 只有在自帶 input/receipt 時才是 science producer；公開 builder 與 HTML 只呈現 validated data，不暗中重算 science |

---

## 核心方法論斷言

這個專案最重要的設計決定，是**什麼東西有資格驅動重建**——而甲基化被刻意排除在外。

![為什麼甲基化不能驅動重建](docs/images/methylation-circularity.png)

當你在某個位點看到兩群甲基化模式不同的 read，至少有四種可能成因：
germline 等位特異性甲基化、雜合性缺失造成的解遮蔽、拷貝數劑量效應、以及真正的譜系差異。
**在目前 single-bulk 測量集合下，若無 orthogonal data 或額外假設，這四者不可識別。**
若要用甲基化獨立「確認」細胞亞群，必須另有外部細胞歸屬；依突變定義群組後觀察到
一致甲基化，只是 concordant association，不是獨立確認。

所以骨幹是**同一條物理分子上的候選體細胞 allele 共現**——它只在「不使用甲基化衍生標籤」
這個限定下**非循環**，仍依賴上游 variant calling、alignment、basecalling 與 haplotag 假設。甲基化被保留為嚴格的 **bounded-auxiliary（有界輔助）**訊號：
它在遺傳候選結構固定**之後**才計算，只做註記，**動不了任何一條邊**。

> 實測數據支持這個克制：811 個可評估的甲基化單元中，
> 只有 **3 個（0.37%）**支持穩健的 pattern-conditioned association。這個低 yield 支持
> 不用甲基化挑拓撲，但不是對所有可能甲基訊號的總檢定。

---

## 實際跑出來的結果

![7 個 technical datasets／6 個 biological IDs funnel](docs/images/funnel-7samples.png)

**7 個資料集（6 個生物樣本）、chr1–22** 的 canonical 結果，2026-08-01 凍結。
每一層都可稽核，且算術自洽
（39,648 + 23,858 = 63,506；+ 8,449 = 71,955；+ 3,224 + 45 = 75,224；+ 10,717 = 85,941；+ 13,014 = 98,955）。

| 項目 | 數值 |
|---|---|
| sSNV 資料列 | 469,849 |
| `k=1` strict read-linkage components | 170,131 / 255,752 strict components（66.5219%） |
| 帶突變的分析單元 | 85,941 |
| 因搜尋節點上限而主動放棄 | 10,717（12.47%） |
| 可用 read-AF 排序的單元 | 71,955 |
| frozen model 對 rankable candidate units 指派單一 **rooted-unlabeled 候選形狀 signature** | 63,506 / 71,955 rankable candidate units（**88.2579%**） |
| **確認的細胞亞群** | **0** — 見下 |
| **確認的線性祖先關係** | **0** — model-selected 候選形狀不是生物 ancestry truth |

> **88.2579% 的正確讀法。**它的意思是：*在已經可排序的 71,955 個單元中*，
> frozen recurrence-allowed model 對其中 63,506 個（88.2579%）指派單一數學候選形狀。這是一個
> **model-conditional 圖形統計**，**不是**「腫瘤演化史已解出 88%」。另一個不同 grain
> 的數字是 170,131 / 255,752 strict components（66.5219%）為 `k=1`；component census
> 無法推得 isolated dataset-record percentage。

---

## 能力邊界

這個專案輸出**機器可讀的宣稱邊界**。凍結的 `canonical/cohort_receipt.json`
承載 `technical_all_pass = true`，而 `summary/all7_summary.json` 承載
`validation_evidence_eligible = false`（它們不是 `authority_manifest.json` 的 top-level 欄位）。
Frozen cohort receipts 記錄的具名 checks 對其引用的 frozen artifacts 皆通過；這不是 current-source
certification，也不是 production／release gate PASS，因此這批結果**還不能作為驗證證據**。

<table>
<tr><th width="50%">✔ 可以主張</th><th width="50%">✘ 明文禁止主張</th></tr>
<tr valign="top"><td>

- 嚴格 read-linked 的**局部**結構
- 家族完整時的完整最小候選家族
- 局部、允許 recurrence 的 Hamming-1 candidate arborescences
- **未經 CN/LOH 校正**的 deterministic read-AF 排序
- frozen-model rooted-unlabeled 候選形狀普查
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

**在目前單一 bulk 的 observation／model，且尚未整合 CN／LOH 的條件下，不能確認
cellular subclone 或 linear ancestry。** 額外的獨立證據（例如 single-cell、multi-region、
orthogonal copy-number 或 purity measurement）可能提高識別性；本專案不宣稱只有某一種
assay 才能突破限制。推論出的內部節點一律標為 `inferred`。

---

## 誠實狀態表

下表是**版本限定**的驗證狀態；「已驗證」只代表所列 artifact 在指定版本實跑或核對，
不是對每一條公開指令與每一個檔案的 blanket guarantee。

| 部件 | 狀態 | 說明 |
|---|:---:|---|
| `inter_sub_mod` | ⚠️ historical dirty audit 通過 | `73afaeac-dirty` 已 build/run，C++/CMake inputs 與 release baseline `ddd8909a` byte-equivalent；這不是 clean-checkout release certification |
| C++ 測試套件 | ⚠️ historical dirty audit 通過 | locked audit 為 0 failure，但 release 仍須 exact-commit、repo-external clean build；test/suite count 必須動態產生 |
| 7 個 technical datasets／6 個 biological IDs 的 HP/PS sidecar | ✅ 齊全 | 7/7 technical datasets PASS |
| research exact-PS topology solver | ✅ 可跑 | 產生上方 exact-PS funnel 的獨立 research solver，不是 `inter_sub_mod` |
| `longlineage preflight` | ✅ 可跑 | 驗證 8 角色 manifest |
| `longlineage dataset-gate` | ⚠️ 受限 | 在 candidate `b9aaa12` 的 frozen HCC1395 receipt 中，唯一觀察到能產生 research artifacts 的介面；**硬編碼綁死單一資料集**，不是 production science validation |
| `longlineage run` / `probe` | 🔴 被鎖 | 依設計回 `KernelBlocked` |
| `longlineage` 拓撲產出 | 🔴 一份 frozen receipt 為 0 個 candidate units | 只限 candidate `b9aaa12` 的 HCC1395 dataset-gate；見[下方說明](#關於-longlineage) |
| 輸出帶標籤的 BAM | ⚠️ 版本限定 | `inter_sub_mod` 不寫 BAM；LongLineage private baseline/main snapshot `5daf50f` 不支援，private public-preview candidate `b9aaa12` 提供 `longlineage-tag-bam`（NOT_READY／non-production） |
| 串起兩支引擎的單一腳本 | 🔴 沒有 | 目前是兩條各自獨立的線 |

### 關於 LongLineage

LongLineage 目前是 **PRIVATE research-preview candidate**。Source-origin、license 與 dependency
audit 尚待完成，因此本交接不先宣稱 clean-room 起源或可公開 redistribution。Immutable
`b9aaa12` 是 public-preview candidate、不是已發布 build；狀態為 `NOT_READY`／non-production，
`P3/P4/P5/P7/P8` 仍 **BLOCKED**。Production `run` 刻意 fail-closed；目前沒有 production tag
或 GitHub Release。`5daf50f` 是 private baseline/main capability snapshot。其契約以 schema 與
SHA-256 receipt 鎖定 artifacts，並記錄具名 abstention；這些契約本身不證明生物正確性。

在 candidate `b9aaa12` 的 frozen **HCC1395 dataset-gate receipt** 中，它輸出 0 個 candidate topology units；
此結果不能外推到所有 LongLineage real-data run：

![LongLineage 漏斗](docs/images/longlineage-funnel.png)

它的流程把拓撲建構**閘控在甲基化分群的穩定性上**，
於是 79,687 個位點中有 66,836 個（83.9%）根本不會產生任何一列共現資料。
這與上面的方法論立場**直接衝突**，任何跨引擎比較都必須明講這一點。

這裡的「BLOCKED」指 `release_attestation.json` 中的 P3/P4/P5，以及 phase ledger／public-safety receipt 中的 P7/P8 具名 gates 尚未通過，
範圍包含 parity／validation、source-origin、license、dependency 與 release-safety。
部分 M1／M2／topology kernels 曾由一份具名 `dataset-gate` receipt 執行；這不等於
全部核心路徑、production entry 或公開 redistribution 已完成。

---

## 快速開始

每一步都附預期輸出的完整流程：[**怎麼跑**](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

![六個步驟](docs/images/howto-six-steps.png)

```bash
# 1. 從 exact checkout 在 clean build directory 編譯（build outputs 不入版控）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 2. 驗證編譯結果   -> 由本次輸出動態取得 test/suite count；要求 exit 0、0 failure
./build/bin/run_tests

# 3. hash-locked Python 3.10 驗收環境
REPO_ROOT="$(pwd -P)"
PYTHON_ENV="$(mktemp -d "${TMPDIR:-/tmp}/ism-python.XXXXXXXX")/venv"
python3.10 -m venv "$PYTHON_ENV"
"$PYTHON_ENV/bin/python" -m pip install --require-hashes \
  --requirement "$REPO_ROOT/requirements-ci.lock"
export PATH="$PYTHON_ENV/bin:$PATH"

# 4. 執行 tracked tiny synthetic fixture（只屬 DEMO/SMOKE，不是科學驗證）。
./scripts/handoff/run_tiny_public_e2e.sh --repo-root "$PWD" --jobs 4

# 5. 真實分析請自行提供有授權且已建立索引的輸入。
./build/bin/inter_sub_mod \
  --tumor-bam /path/to/tumor.mm_ml.bam \
  --reference /path/to/reference.fa \
  --vcf       /path/to/candidates.vcf \
  --output-dir out_min
```

tiny fixture 追蹤 synthetic FASTA／VCF／SAM source，並在新的暫存目錄建立 BAM與index。
PASS receipt 只證明 clone→build→run→schema 的接線可用；scope 是 **DEMO**，不可支持
生物結論、benchmark或release science claim。

<details>
<summary><b>三個一定會咬到你的地方</b>（皆已對照原始碼確認）</summary>

- **`--threads` 預設是 16 不是 1。** help 寫 `Default: 1`，但 `Config.hpp` 的實際預設是 `16`。
  資源估算會整整差 16 倍。
- **未指定 `--distance-metric` 時的 effective CLI default 是 `NHD`。** parser 以 CLI list
  取代 `Config.hpp` initializer；明示的 `NHD`、`BERNOULLI` 或重複多值都會保留。若指定多個
  metric，只有第一個驅動 clustering，後續 metric 只多輸出 distance matrix，所以順序會改變分群結果。
- **`tree.nwk` 的葉子是 *read* 不是 clone。** 它是「read 依甲基化相似度」的階層分群樹，
  **不是**亞群演化譜系樹。這是最常被誤讀的輸出。

另外兩個版本敏感點：`methylation.csv` 第一欄是矩陣列號而非 read 名
（與 read 的綁定完全靠列序，沒有 key 校驗）。Frozen release baseline `ddd8909a` 的
`significance_summary.csv` source header 有 **199 欄**，包含 `VerificationSchemaVersion=2` 與
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
| 🔬 | [上游與資料](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) | Dorado / ClairS / LongPhase-S / SAVANA、sidecar 格式、7 個 technical datasets／6 個 biological IDs |
| 📊 | [分析與呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) | 該用哪些腳本、拒絕渲染的 HTML 生成器 |
| ▶️ | [怎麼跑](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run) | 六個步驟，每步附驗收條件 |

> **同一份內容有兩種讀法。**上表連到 **Wiki**，GitHub 原生渲染、最適合快速瀏覽。
> 同樣的內容也以完全自足的 HTML 形式發佈在 **[liaoyoyo.github.io/InterSubMod](https://liaoyoyo.github.io/InterSubMod/)**
> ——排版更豐富、有可展開的細節區塊；locked 2026-08-12 deploy 共有 **37 個 inline
> SVG elements**（以 DOM `<svg>` 元素計數）。`docs/explain` 是 editorial upstream；
> Wiki 是手動同步的 derivative，可能 drift。發布是獨立步驟，未核對
> deployed commit 前不要假設兩者 byte-identical。

第 01–10 頁涵蓋科學方法本身（名詞地基、ISM 方法、甲基化讀取與篩選、
真實案例逐步拆解、三統計分工、能力 vs 敘述）。

---

## 兩個值得借用的設計

這個 repo 裡有兩個模式可以推廣到基因體學以外。

**1 · 串流是可支援的 workflow mode。**凍結的 9 欄 sidecar 合計 **5.83 GiB**，
保留該 frozen downstream contract 所要求的全部欄位。FIFO mode 可避免在該次 run
新增落地 tagged BAM，不代表歷史上從未存在 tagged BAM。2026-08-13 盤點已找到
7 個歷史 PRE-FIX `paired_full` BAM，合計 1,840,983,466,353 bytes（1.67436 TiB），
另有 7 個 `paired_pileup`；14 個總計 3,709,322,840,333 bytes。Exact paths/bytes/mtimes
與 sampled-content set identity 已知，但未讀全檔 SHA-256，也未完成 generation-level
correspondence，全為 `PARTIAL`/`NON_FINAL`。7 個 `paired_full` stat bytes 除以現行 7 個
sidecar stat bytes 為 294.2669×；這只是跨世代 storage-footprint quotient，不是因果壓縮效果、
byte-equivalent replacement 或 frozen authority。舊 287× 數字是錯的。

![上游工具鏈](docs/images/upstream-toolchain.png)

**2 · 對缺少的必填指標 fail closed。** 經稽核的工作站生成器吃一份宣告式 spec，
當必填指標缺失時**直接拒絕渲染（exit 3）**。這項保護只涵蓋該 generator 與已宣告欄位，
不代表所有公開數字或文件都會自動被驗證。

![拒絕渲染設計](docs/images/workstation-refuse-design.png)

---

## 環境需求

- **C++17** 編譯器、CMake ≥ 3.14、OpenMP
- **htslib**（BAM/VCF I/O）、Eigen
- **samtools**（建立與索引 tracked tiny synthetic fixture）
- **Python 3.10** 驗收環境：以 `--require-hashes` 安裝 `requirements-ci.lock`；
  `requirements.txt` 是 lock 的輸入／一般繪圖依賴清單，不屬 handoff acceptance
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
也 fresh build/test。Public repo **已追蹤 tiny synthetic DEMO fixture**，但真實 tumour BAM／reference／VCF
沒有隨 repo 提供；GitHub About／Wiki／Pages
有各自發布狀態，LongLineage capability 也只適用上面釘定的 commit，因此本文件不宣稱
「每個數字與指令都已驗證」。圖片由 `tools/extract_svg_for_github.py` 從 `docs/explain/`
抽取；應修改來源頁後重生，而不是直接編圖片。

| artefact／claim 家族 | 驗證身分與日期 | 可重跑檢查與結果 | 適用範圍與已知失敗 |
|---|---|---|---|
| Frozen exact-PS funnel | Frozen authority artefacts，2026-08-12 重查 | manifest／hash census 加獨立分母重算；精確計數皆重現 | 只限 frozen 7-dataset analysis；不是 `inter_sub_mod` CLI，也不能識別 cellular clones |
| Frozen release-baseline source | `ddd8909a` | Tracked source header 為 199 欄；current test/suite count 必須由 exact-commit clean run 動態產生 | 歷史 `73afaeac-dirty` audit 的 C++／CMake inputs 與 baseline byte-equivalent，0-failure receipt 只作歷史證據，不是 clean-checkout release certification |
| 公開 quick start | Current source，2026-08-13 檢查 | Build/test 加 tracked tiny synthetic fixture 可重跑 clone→build→run→schema 接線 | Fixture 僅為 DEMO；repo 未附真實 tumor BAM／reference／VCF，因此不宣稱生物結果已重現 |
| LongLineage 狀態 | Private repo、public-preview candidate `b9aaa12`、private baseline snapshot `5daf50f`、frozen HCC1395 evidence | Source／commit 比對加單一資料集 frozen artefact 稽核 | `NOT_READY`；P3/P4/P5/P7/P8 BLOCKED；tagged-BAM 只屬 candidate；尚無 7 technical datasets／6 biological IDs cohort 的 runtime／memory 或多資料集 topology 驗證 |
| GitHub surfaces | Main-repo source 於 2026-08-13 完成本地修正 | 本 correction cycle 的 P0 claim registry 與 source checks | About、Wiki、Pages 與 default branch 發布是獨立外部動作；只改 source 不等於已部署 |

完整命令與實際輸出保存在
[command receipts](research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md)；
逐 claim 修正狀態記錄於
[P0 correction cycle](research/20260813_public_docs_p0_correction/00_INDEX.md)。

**誠實標註的已知缺口**：CN/LOH 未整合進 frozen candidate reconstruction，現有推論也未做
CN/LOH correction；InterSubMod 選填的 LOH BED 僅供 annotation／stratification，不能當作
allele-specific CN、purity、multiplicity、CCF 或 cellular lineage inference。拷貝數整合狀態為
`NOT_INTEGRATED`；
LongLineage 的 7 個 technical datasets／6 個 biological IDs cohort 執行時間與記憶體上界**從未實測**，
且其自身文件明文禁止由單一樣本外推。
