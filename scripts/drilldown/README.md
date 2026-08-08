# drilldown — 多層下鑽觀察儀表板

從全基因組 sSNV 一路下鑽到單一位點，在最細的一層同時看到「該區域的演化分支」
與「甲基的分群狀況」，並用守恆等式自檢整體結果。

完整的資料格式與界線見
`InterSubMod/docs/references/20260808_lineage_pipeline_interface_contract_01.md`。

---

## 快速開始

### 1. 設定站點路徑（一次就好）

```bash
cp drilldown.paths.example.json drilldown.paths.json   # 專案根，已 gitignore
```
```json
{
  "topology_root": "/path/to/<ROUND>/samples",
  "ism_root":      "/path/to/ism_run_root",
  "lineage_root":  "/path/to/lineage_out/<SAMPLE>"
}
```

也可改用環境變數 `DD_TOPOLOGY_ROOT` / `DD_ISM_ROOT` / `DD_LINEAGE_ROOT`，
或每次用 CLI 旗標指定。優先序：**CLI > 環境變數 > 設定檔**。

### 2. 先探測，確認資料齊不齊（不產檔）

```bash
python3 scripts/build_drilldown_dashboard.py --sample HCC1395 --probe-only
```

印出能力矩陣：每一層的狀態、探測停在哪一段、與硬核心的連結率、輸入的 sha256。
**硬核心只有 `topology.jsonl`**；缺它 `exit 3`，其餘缺了只會讓對應面板降級。

### 3. 產頁

```bash
python3 scripts/build_drilldown_dashboard.py --sample HCC1395 --out <DIR>
```

加 `--bake-panels all` 會另外產 16,302 個位點的甲基雙面板 PNG（約 133 MB、
24 核平行約 30 分鐘）。`--bake-panels 200` 可先產一小批試看。

### 4. 開啟

```bash
xdg-open <DIR>/index.html
```

⚠ **必須整個資料夾一起** —— `index.html` 需要同層的 `data/`、`panels/`、`igv/`。
分片走 `<script src>` 而非 `fetch`（Chrome 對 `file://` 的 fetch 一律阻擋）。

---

## 常用旗標

| 旗標 | 說明 |
|---|---|
| `--probe-only` | 只探測不產檔。**每次換資料先跑這個** |
| `--bake-panels {0\|N\|all}` | 甲基雙面板 PNG：不產／前 N 個／全部 |
| `--figs-mode {copy\|link\|ref}` | 預設 `copy`（輸出夾自足可搬走）；`link` 用 symlink，搬移會斷 |
| `--annot-dir <DIR>` | 註釋 drop-in 夾，預設 `<out>/annotations/` |

---

## 註釋 drop-in

把檔案丟進 `<out>/annotations/`，重跑就自動變成篩選維度，**不用改任何程式碼**：

| 格式 | 需要的欄位 |
|---|---|
| `.bed` | `chrom start end [name]` |
| SAVANA / CN segment `.tsv` | `chromosome start end copyNumber [minorAlleleCopyNumber]` |
| 位點表 `.tsv` / `.csv` | 含 chrom + pos（容忍 `chr`/`#chrom`、`pos`/`position` 等寫法） |

檔名去副檔名即維度名（`cgc_genes.bed` → 維度「cgc_genes」）。
解析失敗**不會靜默略過** —— 能力矩陣逐檔說明讀到什麼、命中幾個、為什麼失敗。

---

## 模組結構

```
capability.py      能力探測框架（S0 存在 / S1 結構 / S2 體量 / S3 連結率）
selfcheck.py       12 條守恆等式
cooccurrence.py    維度交集矩陣（Pearson 殘差、n<30 抑制）
sources/           各資料層的探測與載入，每支都可獨立失敗
  topology.py        硬核心
  mlhp.py            read pattern
  ism.py             甲基統計
  strict_edges.py    read 連鎖證據（門檻 what-if 的來源）
  lca_ab.py          LCA 增益 A/B
  annotation.py      drop-in
panels/            甲基雙面板 / IGV 圖的產生（零依賴 PNG 編碼器）
emit/              payload 組裝與 HTML shell
assets/            CSS/JS，build 時 inline
```

**設計原則**：能力缺席就降級並在頁面上寫明原因，不整份拒繪；
所有比例標分母；「未檢定」是獨立的第三態，不等於「不顯著」。

---

## 依賴

- Python 3.9+，**只用 stdlib**（無 numpy / pandas / matplotlib）
- 色盤 import `InterSubMod/scripts/ism_heatmap_std.py`（視覺規約 SoT）
- 驗證用 playwright（選用，只在跑 QA 截圖時需要）

---

## 已知限制

見介面契約文件的「已知限制」章節。最重要的三條：

1. **ISM 覆蓋 81.6% 且非隨機** —— 缺的是被 ISM 自己的 coverage/CpG gate 濾掉的，
   任何「拓撲 × 甲基」的結論都要聲明母體
2. **「甲基自身分群」軸是循環論證**（自檢 C10 會標紅）—— p 值不可當證據
3. **lineage vertex ≠ subclone**，**read 比例 ≠ CCF**
