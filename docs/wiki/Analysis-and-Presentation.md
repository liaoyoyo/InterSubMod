# Analysis & Presentation 分析與呈現層
[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> 部件分冊 · 第 15 頁 · Python 分析層 + HTML 呈現層

C++ 只吐檔案，**把資料變成看得懂的圖與表的是這一層**。這也是實驗室同學實際上花最多時間的地方。本頁挑出真正該用的入口，並誠實說明這層目前的碎片化問題。

| 數字 | 意義 |
|---|---|
| **2,147** | deploy `fbdf7c7` 的全 repo Python 檔（以 `find . -type f -name '*.py'` 計數） |
| **59 / 117** | ⚠️ 讀 CSV 但無任何欄位檢查 |
| **182** | 🔴 硬編了別的分支路徑的腳本 |
| **1** | ✅ 最推薦的第一支腳本（見下） |

---

## 01 先跑這一支 — 重算 historical 35,332-site pipeline 指標

如果你只想重驗 historical 35,332-site script 與它的具名 inputs，跑這個；它不是全系統 health check。

```bash
REPO_ROOT="<REPO_ROOT>"
cd "$REPO_ROOT"
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/verify_pipeline_numbers.py
```

**它會做什麼**：重算 2026-06-27 教材所用的 **historical 35,332-site pipeline** 指標，並與該版文件宣稱值比對。它**不會**驗證 2026-07-24 exact-PS funnel、LongLineage、儲存量、程式碼檔案數或目前測試數。實跑輸出摘要：

```text
sSNV 總數 35,332（TP 30,490 / FP 4,842）      ✓
共現連上   21,554  ·  訊號不足 5,458  ·  孤立 8,320   （加總 = 35,332 ✓）
```

✅ **實跑 exit 0** — 這只證明 historical 35,332-site 指標與其輸入相符，不能外推成「整份 Wiki／所有版本已驗證」。

| Claim family | 應使用的分開驗證來源 | `verify_pipeline_numbers.py` 是否覆蓋 |
|---|---|---|
| historical 35,332-site pipeline | `verify_pipeline_numbers.py` | ✅ |
| 2026-07-24 exact-PS funnel | `docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`、`authority_manifest.json` 與對應 solver receipts | ❌ |
| LongLineage 功能／版本 | 對指定 commit 分別跑 `preflight`／`dataset-gate`／binary inventory | ❌ |
| 儲存量、程式碼數、測試數 | 在目標 commit／機器上分別做檔案統計、source inventory 與 fresh test run | ❌ |

> **建議的閱讀順序**
> 先跑上面這支確認 historical 35,332-site script 與具名 inputs 在本環境相符 → 再看 §02 挑你需要的分析腳本 → 最後用 §03 的工作站 HTML 看結果。其他 claim family 必須分別驗證。

---

## 02 主力腳本 — 各自產出什麼

從 deploy `fbdf7c7` 的 2,147 個 Python 檔中挑出真正是「入口」的那幾支；此數字不代表 production entry-point 數。

| 腳本 | 產出什麼、回答什麼問題 |
|---|---|
| `verify_pipeline_numbers.py`<br>✅ **先跑這個** | 重算 historical 35,332-site pipeline 的 TP／FP／共現／訊號不足／孤立等指標；**不覆蓋** exact-PS、LongLineage、儲存量、code count 或 test count。 |
| `sm_linkage_genomewide.py` | 建立**全基因組突變共現地圖** —— 這是整套方法的骨幹。對每一對距離 50 kb 內的突變，統計有幾條 read 同時看到它們、四種等位組合各幾條、判成巢狀／互斥／共連鎖／獨立哪一種關係。<br>*實測產物 26.5 MB：53,094 對配對 + 35,332 筆位點帳本。* |
| `build_strict_ps_hp_regions.py` | 用**嚴格規則**切出可建樹的區域：同一定相組、同一 germline 單倍型內、至少 3 條分子同時確定兩端 —— 才連一條邊。<br>🔴 **距離只記錄、不參與連邊**，這是刻意的方法學選擇（避免用幾何鄰近冒充分子連鎖）。 |
| `build_layered_per_sample.py` | 產生 7 technical datasets／6 biological IDs 的**分層工作站 HTML**（呈現層主成品）。 |
| `build_exact_ps_layered_workstation.py` | 4,741 行的完整工作站建構器。支援 `--verify-only` —— **只驗證上游契約還成立、不重建**，適合用來確認資料沒漂移。 |
| `build_workstation.py`<br>**通用** | 317 行的**通用**工作站生成器：吃一份宣告式 spec，吐互動判讀 HTML。見 §03。 |

<details>
<summary>關鍵中間產物 — 想直接分析資料時該讀哪些檔</summary>

| 檔案 | 內容 |
|---|---|
| `per_sSNV_census.tsv` | 每個突變一列的帳本：來源、正常組織的支持數、是否確認為體細胞突變、變異頻率、50 kb 內有幾個潛在夥伴、實際連上幾個、最高共讀深度、拷貝數狀態。 |
| `regions.tsv` | 在 frozen recurrence-allowed revision 中，每個共現區域一列：座標、跨度、突變數、local candidate graph／arborescence 的形狀摘要、節點與相容／矛盾計數。這不是 biological tree、ancestry 或 clone truth，且尚未做 allele-specific CN／LOH correction。 |
| `funnel_census_HCC1395.json` | **所有百分比的分母來源**。六層漏斗，每層丟了多少、為什麼丟，並**自我檢查各層加總等於總數**。要引用任何比例前先看這個檔。 |

</details>

---

## 03 HTML 工作站 — 人真正看的東西

產出的是**零外部依賴、可離線開啟**的單一 HTML 檔，可以直接寄給別人。

### 圖 1 · 工作站生成器的「拒絕渲染」設計 — 防止必填欄位被靜默省略

![workstation-refuse-design](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/ddd8909a838318d8a77969313e9561c8ff9d01c2/docs/images/workstation-refuse-design.png)

> 這個設計來自一次真實事故：曾有報告把「預期數字」當成真實結果寫進 HTML，而分析其實沒跑完。此 generator 能防止已宣告必填欄位被靜默省略；**它不保證數值來源或科學解讀正確**。
>
> 🔴 但目前這個通用生成器**只有 2 個真正的複用者** —— 實務上大家用的是 4,741 行的專用版本。這是這一層最明顯的技術債。

圖中流程（逐字對應）：

| 階段 | 內容 |
|---|---|
| **宣告式 spec** | 列出每個待判讀項目與它的必填指標清單 |
| **生成器逐項檢查** | 每個項目是否都具備 spec 宣告的必填指標？（317 行，無外部相依） |
| ✅ **齊全 → 渲染互動 HTML** | 人工判讀存入瀏覽器本機儲存、可匯出 JSON/CSV；含來源標記、修正歷程、點圖放大 |
| 🔴 **缺任一必填 → 退出碼 3，拒絕渲染** | 不是填破折號、不是留空白 —— 是整個拒絕產出（本輪以自建最小 spec 實測確認） |

⚠️ **為什麼要這樣設計**
如果缺資料時填一個破折號，報告看起來仍然完整，讀的人不會發現那一格其實沒有依據。
🔴 「拒絕產出」阻止已宣告必填欄位缺失時仍生成完整外觀；它不涵蓋未宣告欄位，也不保證科學正確性。

### 現成的工作站成品

| 檔案 | 大小 | 備註 |
|---|---|---|
| `HCC1937.html` | 14.6 MB | ✅ **示範給人看時先開這個** —— 最小、開得動 |
| `HCC1395.html` | 35.5 MB | 主力樣本 |
| `COLO829.html` | 41.5 MB | 有外部真值可對照 |
| `H1437.html` | 72.3 MB | |
| `H2009.html` | 188.2 MB | 🔴 **瀏覽器開這個會很慢且吃大量記憶體** |

7 個工作站合計約 393 MB，**不進 git**（該目錄有自己的忽略規則）。每個頁面的 HTML 標頭內嵌了上游各收據的 SHA-256 當作防漂移印記。

### 兩種主要的圖

| 圖種 | 怎麼看 |
|---|---|
| **雙面板熱圖** | 左邊是 read × CpG 甲基化矩陣（藍＝未甲基／紅＝甲基／灰＝無資料），右邊是同一批 read 的兩兩距離矩陣（暗＝相近／亮＝相遠）。左側固定側欄依序是 **群 ｜ 單倍型 ｜ 突變型 ｜ 腫瘤/正常 ｜ 正反股**。<br>**要看的是**：同一群 read 是不是同時在「甲基化型態」與「帶不帶突變」上都一致。 |
| **全基因組染色體圖** | 22 條染色體橫條，用紅色深淺表示突變密度，疊上著絲點位置。<br>**用途**：目視確認密集區是不是落在已知的假象好發區。 |

---

## 04 這一層的四個坑（誠實揭露）

這是全系統目前最碎片化的一層，以下問題都是實測確認的。

### 🔴 坑一：預設的 python3 版本會讓某些腳本直接崩

預設 `python3` 是 **3.9.12**，但部分腳本需要 **3.10**。照著說明打 `python3 scripts/build_strict_ps_hp_regions.py` 會在 import 階段就崩，**而且錯誤訊息指向 dataclass 而不是版本問題** —— 極容易誤以為程式壞了。

**解法**：改用 `/usr/bin/python3.10`。
*（`run_layered_v4_strict.py` 有自動偵測會幫你處理，但直接呼叫子腳本時沒有這層保護。）*

### 🔴 坑二：182 支腳本硬編了「另一個分支」的路徑

這些腳本硬編了一個 git worktree 的絕對路徑，而**那個 worktree 目前 checkout 在另一個分支上**（不是主工作分支）。

也就是說**這些圖與中間結果實際上是從另一個分支的檔案樹讀出來的**。那個 worktree 一旦被移除或切換分支，這 182 支腳本會全部靜默讀到不同內容，或直接找不到檔案。

### 🔴 坑三：名字像入口、其實是函式庫

`scripts/ism_heatmap_std.py` 看起來像可執行的入口，但它**沒有 main、沒有參數解析** —— 直接跑它**什麼都不會發生**（它只定義常數與函式）。

真正畫圖的是 import 它的那些 renderer。它是熱圖規格的**單一真實來源**（側欄順序、配色），有專門的稽核腳本檢查其他頁面有沒有偏離這個配色。

### ⚠️ 坑四：資料讀取普遍沒有欄位檢查

`scripts/` 下 **291** 個 Python 檔（deploy `fbdf7c7`；含 generated/test helpers）中，117 個會讀 CSV／TSV，其中 **59 個**原始碼裡**沒有任何 schema 或欄位檢查**。不同母體不可直接用「超過一半」解讀。

配合 [ISM 分冊](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine)提到的「`significance_summary.csv` 欄數跨版本不一致」，這代表**拿舊資料餵新腳本時很可能靜默讀到錯的欄位**。

**自保做法**：自己寫分析時一律**用欄名取值**，不要用欄號。

---

## 本頁的驗證方式

- **腳本可用性**：主力腳本的 `--help` 本輪實跑取 exit code；`verify_pipeline_numbers.py` 完整跑完並讀回 historical 35,332-site 輸出。此驗證不涵蓋其他 claim family。
- **拒絕渲染行為**：以自建的最小 spec（故意缺一個必填指標）實測，確認退出碼為 3。
- **檔案大小與腳本計數**：實際 `ls` 與遞迴統計，非估計值。
- **python 版本問題**：以預設 `python3` 實跑重現崩潰，再以 3.10 實跑確認正常。

**路徑**：`InterSubMod/docs/explain/15_python-html-layer.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06

---

[← 上一頁：Upstream & Data](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [回 System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [下一頁：How to Run →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)
