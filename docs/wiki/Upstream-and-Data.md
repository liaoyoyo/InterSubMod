# Upstream Toolchain 上游工具鏈與資料
[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> **兩支 C++ 的輸入是怎麼來的？**
>
> ISM 與 LongLineage 都不是從原始訊號開始的 —— 它們吃的是**已經被四個外部工具處理過**的資料。這一層**不是本專案寫的程式**，但它決定了所有下游分析的品質上限，而且藏著幾個會讓人得出相反結論的陷阱。

---

## 四樣輸入、四個來源

| 輸入 | 產生者 | 它回答什麼 |
|---|---|---|
| **甲基化訊號**<br>`MM`／`ML` tag | Dorado basecaller | 這條 read 上每個 CpG 有沒有甲基化、機率多高。**資料進來時就已經在 BAM 裡了。** |
| **體細胞突變清單**<br>VCF | ClairS（配對模式） | 哪些位置是腫瘤才有、正常組織沒有的突變。 |
| **單倍型標籤**<br>`HP`／`PS` | LongPhase-S | 這條 read 來自父源還是母源那一套染色體。**本交接的 frozen FIFO workflow 抽成 sidecar，該次 run 不新增落地 BAM；歷史 PRE-FIX workflows 曾落地 tagged BAM。** |
| **拷貝數**<br>CN 區段 | SAVANA（僅 cna） | 哪些區段被複製或缺失。⚠️ **目前狀態是 `NOT_INTEGRATED`，尚未接入主線。** |

| 指標 | 數值 |
|---|---|
| ✅ 樣本 sidecar 全部生產完成 | **7 / 7** |
| 7 technical datasets／6 biological IDs 的 sidecar 總量（壓縮後） | **5.83 GiB** |
| 7 `paired_full` + 7 `paired_pileup` PRE-FIX tagged BAM | `paired_full` = **1,840,983,466,353 bytes（1.67436 TiB）**；14 個總計 **3,709,322,840,333 bytes**。Paths/bytes/mtimes 與 sampled set identity known；full-file hashes/provenance correspondence 未閉合，全為 `PARTIAL`/`NON_FINAL` |
| 7 `paired_full` stat bytes ÷ 7 sidecar stat bytes | **294.2669×**；僅為跨世代 storage-footprint quotient，不是因果 compression/reduction、byte-equivalent replacement 或 frozen authority；舊 287× 錯誤 |

---

## 01 · 完整前處理鏈（8 個步驟）

從原始訊號到兩支 C++ 吃得下的東西。步驟 ① 通常在資料交付時就已完成。

![upstream-toolchain](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/a34b0cb96a8ef247c5a6f423d46b2920c7e541aa/docs/images/upstream-toolchain.png)

> **圖 1 · 上游前處理鏈，以及 FIFO 如何避免該次 run 新增落地 BAM**
>
> sidecar 把「下游 frozen contract 要求的欄位」與 BAM 容器分開；可說保留該九欄 contract 所需欄位，不能說保留所有可能有用資訊。7 個 paired_full 與 7 個 sidecar 的 exact stat bytes quotient 是 294.2669×，只代表跨世代 footprint；不是壓縮率或內容等價，舊 287× 無效。

### 八個步驟逐項

| 步驟 | 做什麼 | 細節 |
|---|---|---|
| ① Dorado basecalling | 電訊號 → 鹼基 + 甲基化機率 | 產出 `MM` / `ML` tag |
| ② germline 定相 | Clair3 叫變異 + LongPhase 定相 | 在 normal BAM 上做 |
| ③ ClairS 體細胞突變 | 腫瘤 + 正常配對比較 | `134,122` 筆 |
| ④ VCF 正規化 | 只改一個欄位型別，不動任何 FILTER | 守恆稽核：進出筆數必須相同 |
| ⑤ **LongPhase-S — 這裡產生 `HP` / `PS` 標籤** | 以 germline 單倍型為骨架，判斷每條 read 屬於 HP1 還是 HP2；若該 read 同時帶體細胞變異，再細分成 `1-1` / `2-1`（germline 骨架上的體細胞支系）或 `3` | 🔴 執行參數：12 執行緒、MAPQ ≥ 20、標記補充比對、輸出體細胞 VCF |
| ⑥ **FIFO mode：tagged BAM 寫進具名管道，避免該次 run 新增落地 BAM** | 見下方展開 | ✅ |
| ⑦ FILTER 雙向重校正稽核 | 比對正規化前後的 VCF，逐位點統計狀態轉移；並確認除 FILTER 外其他欄位未被改動 | ⚠️ 見下方 §03 — 這裡有個容易算錯的地方 |
| ⑧ SAVANA 拷貝數 | 在雜合位點數等位比例，擬合純度與倍性；刻意只跑 `cna` 子命令 | 🔴 跑完整流程會吃到 400–488 GB 記憶體近乎 OOM |

### ⑥ 具名管道（FIFO）這一步

```text
LongPhase-S 寫出                具名管道          另一程序邊讀邊抽              sidecar TSV
tagged BAM 串流       ───►       .fifo      ───►  每條 read 抽 9 個欄位   ───►   5.83 GiB
（若落地 = 每 technical dataset 上百 GB）                 CIGAR 只存 8 位元組摘要        （7 technical datasets 合計）
```

- 這個 FIFO mode 可讓該次 run 不新增落地 BAM；不代表歷史 canonical 從未有 tagged BAM。
- 接著驗證：列數必須等於比對數、不得有未知的 `HP` 值、座標必須仍然有序，任一不過就中止。

### 為什麼要大費周章用管道，而不是存 tagged BAM？

- frozen sidecar 為 **5.83 GiB**（壓縮後）／13.98 GiB（未壓縮）。
- 2026-08-13 盤點識別 7 個 PRE-FIX `paired_full` BAM，exact paths/bytes/mtimes 合計 **1,840,983,466,353 bytes（1.67436 TiB）**；另有 7 個 `paired_pileup`，14 個總計 **3,709,322,840,333 bytes**。
- 7 個 historical `paired_full` 有 sampled-content identity set SHA-256，但不是 7 個 full-file SHA-256；producer/generation correspondence 仍未逐檔閉合，故保持 `PARTIAL`/`NON_FINAL`。
- `paired_full` exact stat bytes ÷ 7 sidecar exact stat bytes（6,256,168,164）= **294.2669×**；這只是跨世代 footprint quotient，不是壓縮效果，舊 287× 錯誤。
- SEQ/QUAL未納入九欄 sidecar contract；在可重現 byte decomposition 前，不宣稱其占比 >99%或用途為0%。
- ⚠️ 代價：無法直接用 IGV 之類的工具看標籤（見下方 §04）。

---

## 02 · sidecar 長什麼樣 — exact-PS／LongLineage 的 commit-scoped contract

InterSubMod 從 BAM aux tags 讀 HP/PS；exact-PS／LongLineage 使用 sidecar。兩者是平行
provenance contracts，**不是**兩支引擎直接串接的執行期介面。LongLineage preview 的
tagged-BAM writer 只可在明示 branch/commit 下宣稱，不能由這份 sidecar 推成 supported E2E。

```text
#CHROM	START0	END0	QNAME	FLAG	MAPQ	CIGAR_B2	HP	PS
```

| 欄位 | 內容 |
|---|---|
| `QNAME` | read 的原始名稱，實測是明文 UUID（如 `ec1d98d5-5fcf-45fc-9691-7641e8731385`）。 |
| `CIGAR_B2` | 比對形狀的 **8 位元組摘要**（不是完整 CIGAR）。省空間，但仍足以當比對用的鍵。 |
| `HP` | 單倍型標籤。值域：`1`／`2`（germline）、`1-1`／`2-1`（帶體細胞變異的支系）、`3`、`.`（未定相）。 |
| `PS` | 定相組編號 —— 同一組內的單倍型標籤才可互相比較。 |

### 7 個 technical datasets／6 個 biological IDs 的實際規模

| 資料集 | sidecar 列數 | 備註 |
|---|---|---|
| HCC1395 | `40,859,727` | 主力樣本，1.43 GB（壓縮後） |
| HCC1395_DORADO | `40,033,094` | **同一個生物樣本的另一次 basecalling**，不是生物學重複 |
| COLO829 | `8,255,461` | 有外部真值可對照 |
| H1437 | `14,434,968` | 純度 0.95 |
| H2009 | `21,690,297` | 純度 0.95 |
| HCC1937 · HCC1954 | `—` | HCC1937 的 BAM 是最大的（416 GiB） |

✅ **7/7 PASS** — 全部生產完成並通過驗證（2026-07-11 起跑，隔日 01:31 全部驗證通過）。**7 個生產資料集 = 6 個生物樣本**（HCC1395 有兩套 basecalling）。

---

## 03 · 三個會讓你算錯數字的陷阱

### 🔴 ● 陷阱一：突變數量會對不上，因為重校正是**雙向**的

LongPhase-S 不只把「原本不合格」的位點救回來，**也會把原本合格的降級**。如果誤以為重校正只是「加分」，就會把重校正後的合格集當成原本合格集的**超集** —— 那是錯的。

> HCC1395 實測：救回 **4,592** 個、降級 **5,528** 個 → 淨變化 **−936**（113,997 → 113,061）。
> 用錯集合，骨幹突變數就會對不上。七個樣本的救回／降級數各不相同。

### 🔴 ● 陷阱二：`HP` 標籤有兩種資料型別，用錯 parser 會靜默丟資料

同樣叫 `HP` 的欄位：LongPhase-S 寫的是**字串**（`HP:Z:1-1`），另一個工具寫的是**整數**（`HP:i:11`），而且整數編碼是 `11↔"1-1"`、`21↔"2-1"`、`33↔"3"` 這種非直覺的映射。

> **用錯 parser 不會報錯** —— 只會讓所有「帶體細胞變異的支系」read 消失在分組之外，看起來就像**「這個樣本沒有亞群訊號」**。這是最容易得出相反結論的一個坑。

### 🔴 ● 陷阱三：`1-1`／`2-1` 不等於「已確認的亞群」

這兩個標籤的定義本身**就是用體細胞變異切出來的** —— 它們是「在 germline 骨架上再被體細胞證據細分出來的 read 群」。

> 因此**拿它去「驗證」體細胞變異的分群，就是循環論證**。ISM 甚至內建一個開關，專門把這些標籤降級成「未定相」來避免這個循環。
>
> **正確用法**：HP 是**鑑別器**（幫忙排除「其實只是來自不同那套染色體」的假象），**不是確認器**。

<details>
<summary>另外兩個較次要但值得知道的限制</summary>

⚠️ **比對用的鍵不含序列與品質**：sidecar 與 BAM 對回去時用的是「read 名 + 染色體 + 起訖 + FLAG + CIGAR 摘要」。設計文件明寫這只能在特定條件全部成立時使用 —— 若上游 BAM 換版本、重新比對、或 read 名重複，同一把鍵可能對到多筆。HCC1395 實測確實有 6,594,614 筆重複的完全相同比對列，**但目前沒有衝突**。

⚠️ **用的 LongPhase-S 是自訂修改版**，未經上游審查。原版在某些情況會直接中止整個程序，修改版改成印警告後繼續並保留該位點。這改變了低品質位點的處置語意。HCC1395 出現 2 次該警告，其餘 6 個資料集為 0。收據自己標註了「證據是有界的」。

</details>

---

## 04 · 如果想要「帶標籤的 BAM」怎麼辦？

這是被問最多的需求之一，答案有點反直覺。

> **現況**
>
> `inter_sub_mod` 沒有 BAM writer。LongLineage 則必須標明版本：private baseline/main snapshot **`5daf50f`** 沒有 tagged-BAM writer，但 private public-preview candidate **`b9aaa12`** 已含 `longlineage-tag-bam`。Candidate 仍是 `NOT_READY`／non-production，`P3/P4/P5/P7/P8` BLOCKED，能力尚未發布。現存生產檔案究竟由 LongPhase-S 或哪個 revision 產生，必須依 dataset provenance／receipt 判定。

✅ **好消息是技術上完全可行**，因為比對用的鍵是確定性的正向計算：

- sidecar 裡的 `QNAME` 是**明文**的 UUID
- LongLineage 的 `read_id` 是它的 **SHA-256 雜湊**
- 所以讀 BAM 時對每條 read 名重算一次雜湊就能對上，**不需要任何反查表**

若以 private baseline `5daf50f` 行為作參考，可行做法仍是獨立匯出工具（讀凍結結果 + 原始 BAM，重算雜湊做比對）；若評估 candidate `b9aaa12` 的 `longlineage-tag-bam`，必須另做該 revision 的 schema、provenance 與輸出驗證，不能把 candidate 能力寫成已公開或 production 可用。

🔴 ● **磁碟規劃仍是 blocker**：14 個 PRE-FIX BAM 合計 3,709,322,840,333 bytes，但全檔 hash、byte decomposition、generation correspondence 與 production eligibility 都未驗。實務上建議只針對感興趣區間產生子集。

---

## 本頁的驗證方式

- **sidecar 格式**：實際 zcat 取得表頭與資料行，與 LongLineage 原始碼中的常數逐位元組比對一致。
- **cohort 狀態**：實際讀取生產狀態表，7/7 technical datasets 皆為 PASS；對應 6 個 biological IDs。
- **FILTER 轉移數**：取自稽核 JSON，並驗算 4,592 − 5,528 = −936 與前後總數相符。
- **磁碟數字**：實際檔案位元組數加總。

## 誠實標註

- ⚠️ 步驟 ① basecalling **不是本專案執行的**，資料進來時已完成，本頁描述取自資料集的來源記錄。
- ⚠️ SAVANA 的拷貝數結果目前狀態為 `NOT_INTEGRATED` —— **尚未接入主線分析**，因此所有現有結論都是「未經拷貝數校正」的。

---

**來源**：`InterSubMod/docs/explain/14_upstream-data.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06

[← 上一頁：LongLineage 部件](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [回系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [下一頁：分析與呈現層 →](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation)
