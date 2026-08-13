<!--
建立時間: 2026-08-05
狀態: HISTORICAL_SUPERSEDED（僅代表 2026-08-05 當時快照曾被驗證）
目標: 工程交接包的單一入口 — 30 秒結論、閱讀順序、第一天該做什麼
處理範圍: 資料層 / 程式層 / 輸出層 / 收斂驗收；補 20260801 科學交接包未涵蓋的工程面
關聯檔案:
  - InterSubMod/docs/handoff/20260805_系統交接與驗收_01/handoff_manifest.json
  - InterSubMod/docs/handoff/20260805_系統交接與驗收_01/acceptance_receipt.json
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/README.md（科學層 parent）
驗證方式: 2026-08-05 實跑 clean build + C++ 258 tests + 19 artifact hash + 資料量測
-->

> [!CAUTION]
> **HISTORICAL / SUPERSEDED — 只用於追溯 2026-08-05 當時觀察。** 本包中的磁碟餘量、測試數、dirty-file 數、CLI/schema 欄數與 readiness 判定都不是 2026-08-13 現況，也不得覆蓋 frozen science authority。請從[2026-08-13 完整研究交接](../20260813_完整研究資料與軟體交接_01/00_INDEX.md)進入；原始 JSON receipts 保留為歷史證據，未被改寫。

# InterSubMod 系統交接與驗收包（2026-08-05）

> **30 秒結論**
>
> 科學證據鏈**完好**（19/19 artifact 雜湊四天後仍逐位元一致），程式**可從零編譯**（exit 0），
> C++ 測試 **258/258 全過**。但**現在還不能正式交接**：44% 的 Python 測試因環境缺 pytest 跑不了、
> 磁碟只剩 617 GB 無法 rerun、580 項未提交讓版本無法指認、新人文件路徑是錯的。
> 這四項**全是工程問題，不是科學問題**，估計 1–2 個工作天可全部修完。

---

## 1. 這個包補了什麼

已經有一份 **2026-08-01 的科學交接包**
（`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/`），
它涵蓋科學結論、claim ceiling、分母定義、authority artifacts。

**本包補的是它沒有的工程面**：

| 面向 | 20260801（科學） | 本包（工程） |
|---|---|---|
| 結論可不可信 | ✅ | — |
| claim 邊界 | ✅ | — |
| 數字的分母 | ✅ | — |
| **怎麼 build、怎麼跑** | ✗ | ✅ |
| **資料在哪、多大、哪些可刪** | ✗ | ✅ |
| **輸出目錄語意與信任度** | ✗ | ✅ |
| **驗收怎麼判定通過** | ✗ | ✅ |
| **收斂待辦的執行順序** | 有 P0 大方向 | ✅ 細到可勾選 |

> **規則**：科學數字一律回 20260801 那包查，本包**不建立第二份科學數值 SoT**。

---

## 2. 閱讀順序

### 如果你是接手者（要能跑起來、能續開發）

1. 本 README（你在這）
2. **[00_交接總覽與驗收.standalone.html](00_交接總覽與驗收.standalone.html)** — 主視覺文件，含四層架構與驗收儀表板的 SVG 圖，離線單檔可讀
3. [02_程式層.md](02_程式層.md) — build、執行、程式地圖、環境缺陷 ← **第一天最重要**
4. [01_資料層.md](01_資料層.md) — 資料在哪、多大、哪些可信
5. [03_輸出層.md](03_輸出層.md) — 輸出怎麼看、receipt 體系
6. [04_收斂驗收清單.md](04_收斂驗收清單.md) — 7 個 gate + 9 步收斂順序
7. 然後才讀 [20260801 科學交接包](../20260801_exactPS_readAF_CNV_AI交接_01/README.md)

### 如果你是驗收者（教授／口試委員）

1. 本 README §3「目前狀態一覽」
2. [00_交接總覽與驗收.standalone.html](00_交接總覽與驗收.standalone.html) 的驗收儀表板段落
3. [04_收斂驗收清單.md](04_收斂驗收清單.md) §6「給驗收者的三個提問」
4. 科學結論看 [20260801 交接手冊](../20260801_exactPS_readAF_CNV_AI交接_01/20260801_研究狀態與AI交接手冊_01.md)

### 如果你是新接手的 AI

先讀 `handoff_manifest.json` 與 `acceptance_receipt.json`（機器可讀），
再讀 20260801 的 `authority_manifest.json`。不要只讀某一段散文就自行推論狀態。

---

## 3. 目前狀態一覽（2026-08-05 實測）

| Gate | 判準 | 狀態 | 實測值 |
|---|---|:---:|---|
| **G1** 可編譯 | 從零 clean build 成功 | ✅ | configure/build 皆 exit 0，6/6 binary |
| **G2** 證據鏈完整 | authority artifact 雜湊一致 | ✅ | **19/19 MATCH** |
| **G3** 測試可執行 | 全部測試可跑且結果相符 | ❌ | C++ 258/258 ✅；Python **63/143 跑不了** |
| **G4** 資料可定位 | 路徑/大小/信任度清楚 | ⚠ | 7/7 BAM 在（1.67 TiB）；但 manifest completeness 欄位過期 |
| **G5** 有作業餘裕 | 磁碟 ≥2 TB 可 rerun | ❌ | **僅剩 617 GB（99% 已用）** |
| **G6** 新人可上手 | 照文件能一天上手 | ❌ | QUICKSTART 路徑寫成 `/big8_disk`（實際 `/big7_disk`） |
| **G7** 版本可指認 | 有不可變版本識別 | ❌ | **580 項未提交**；frozen binary `not_proven_byte_reproducible` |

---

## 4. 第一天該做什麼

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

# ① 確認能 build（約 10-20 分鐘，含 jemalloc）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)

# ② 確認 C++ 測試全過（應為 258 tests PASSED）
./build/bin/run_tests

# ③ 確認證據鏈仍完整（19/19 應全部 MATCH）
#    做法見 04_收斂驗收清單.md §1 G2

# ④ 確認磁碟餘量（跑任何長 pipeline 之前必做）
df -h /big7_disk
```

然後讀 `02_程式層.md` §6 的「已知缺陷」五條，理解為什麼還不能正式交接。

---

## 5. 三個最容易在交接時出事的點

1. **frozen binary 在別的磁碟。**
   `/bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af`
   （SHA-256 `ba13ccf2...33b2e`）產生了論文核心數字，卻**不在 repo 也不在 /big7_disk**。
   搬遷或清理時極易遺失。

2. **資料不能複製。**
   光 `paired_full` 的 7 個 tagged BAM 就 **1.67 TiB**，而磁碟只剩 617 GB。
   交接的是路徑 + checksum + 重生方法，**不是檔案**。

3. **technical PASS ≠ 科學完成。**
   最新 release 同時是 `report_build_status=PASS` 和
   `scientific_completeness_status=INCOMPLETE_WITH_ABSTAIN`，兩者不矛盾。
   `confirmed cellular subclone = 0`、`linear ancestry = 0` 的 claim ceiling 不變。

---

## 6. 本包的驗證方式

所有數字都來自 2026-08-05 的實際執行，不是清單或估計：

| 驗收項 | 方法 |
|---|---|
| 編譯 | 新建 `build_handoff_verify/` 做 from-scratch build（不動既有 `build/`） |
| 證據完整性 | 逐檔重算 SHA-256 與 `authority_manifest.json` 比對 |
| C++ 測試 | 執行本輪 build 產出的 `run_tests`（非舊 binary） |
| Python 測試 | `python3 -m unittest` 跑 20260801 release 測試 |
| 資料規模 | `ls -la` 逐一量測 7 個 tagged BAM |
| 磁碟 | `df -h` |

完整機器可讀結果：[`acceptance_receipt.json`](acceptance_receipt.json)

---

## 7. 更新規則

1. 本包的驗收數字**只能由重跑更新**，不可手改。
2. 每次重跑後更新 `acceptance_receipt.json` 的 `probe_date` 與對應 check。
3. Gate 由 FAIL 轉 PASS 時，必須同時在 `handoff_manifest.json` 的
   `convergence_order` 標記該步驟 `status: DONE_<date>` 並附實測結果。
4. 科學數字有更新時，**改 20260801 那包**，本包只引用不複製。
