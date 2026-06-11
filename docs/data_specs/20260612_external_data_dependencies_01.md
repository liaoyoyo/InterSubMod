<!--
建立時間: 2026-06-12
報告類型: 外部數據依賴契約（External Data Dependency Contract）
任務類型: 數據可尋性 + SPOF 防護（改進 ④）
狀態: current（所有 ls -l 經 2026-06-12 fresh 驗證）
data_sources: 主回合 2026-06-12 ls -l /big8_disk/data/{6 samples}/.../BL.bam fresh 驗證輸出
驗證方式: 逐檔 ls -l（owner/size/mtime/symlink target）；甲基 tag 狀態為 dir-name 推論，標 未確認
-->

# 外部數據依賴契約 — matched normal BAM（big8 跨帳號）

> **本檔職責**：把「不在我方寫權範圍、但研究關鍵」的外部數據檔登錄成契約，記 owner/path/size/last-seen，並標 SPOF 風險與啟用前必驗項。
> **為何需要**：V10 copy-clean 陰性對照 + cis-test 依賴 matched normal 甲基；這些檔全在他人帳號（`zhenyu112`）的唯讀 NFS，對方清理/搬移即斷鏈且**我方無權挽救**。

---

## L0 一眼結論

- **6 樣本 matched normal BAM 全部存在**，但全是 **`zhenyu112` 帳號**的檔、big8 唯讀 NFS、mtime 都是 ~1 年前（2025-02~03）。
- **🔴 跨樣本 G-A blocker**：**只有 HCC1395 確認有甲基（5mCG+5hmCG）normal**。其他 5 樣本的 normal 是 plain `_Normal_ONT`（4 個還是 symlink 指到第三處），**甲基 tag 未確認** → 跨 6 樣本 V10（matched-normal 甲基 copy-clean）對 **5/6 樣本卡住**，須先驗 MM tag 或重 basecall。
- **🔴 雙重 SPOF**：4/6 是 symlink → `/big8_disk/Google_somatic_data/bams/`（另一目錄樹），原檔或 symlink 任一動就斷。

## L1 依賴清單（fresh 驗證 2026-06-12）

| 樣本 | normal BAM 路徑 | 形式 | owner | size/mtime | 甲基(MM/ML)? |
|---|---|---|---|---|---|
| **HCC1395** | `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam` | 實檔 | zhenyu112 | 149.4 GB / 2025-02-25 | ✅ **5mCG+5hmCG**（V10 已用） |
| HCC1937 | `/big8_disk/data/HCC1937/ONT/HCC1937BL.bam` | **symlink** → `/big8_disk/Google_somatic_data/bams/HCC1937/HCC1937_Normal_ONT.GRCh38.sorted.bam` | zhenyu112 | 2025-03-04 | ⚠ **未確認**（名為 _Normal_ONT，非 5mCG dir） |
| HCC1954 | `/big8_disk/data/HCC1954/ONT/HCC1954BL.bam` | **symlink** → `.../bams/HCC1954/HCC1954_Normal_ONT.GRCh38.sorted.bam` | zhenyu112 | 2025-03-04 | ⚠ 未確認 |
| H1437 | `/big8_disk/data/H1437/ONT/H1437BL.bam` | **symlink** → `.../bams/H1437/H1437_Normal_ONT.GRCh38.sorted.bam` | zhenyu112 | 2025-03-04 | ⚠ 未確認 |
| H2009 | `/big8_disk/data/H2009/ONT/H2009BL.bam` | **symlink** → `.../bams/H2009/H2009_Normal_ONT.GRCh38.sorted.bam` | zhenyu112 | 2025-03-04 | ⚠ 未確認 |
| COLO829 | `/big8_disk/data/COLO829/ONT_R10/COLO829BL.bam` | 實檔 | zhenyu112 | 141.6 GB / 2025-03-02 | ⚠ 未確認（R10 flowcell） |

> 其他 5 樣本**皆無** `*5mCG*` / `*5khz*` 甲基變體 normal 目錄（2026-06-12 ls 確認為空）。

## L2 風險與意涵

1. **SPOF（單點失效）**：6 檔全 `zhenyu112` 擁有、我方唯讀。對方清理/改權限/搬 symlink target → V10/cis 對照整條斷且無法自救。
2. **跨樣本 tier 解鎖前置**：foundation doc §5 G-A「跨 6 樣本重現 V10」假設資產已備，但**甲基 normal 只有 1/6 確認** → G-A 真正可跑前需先解這個資料 gap。這是 tier ⭐3→⭐4 的隱藏 blocker。
3. **symlink 脆弱**：4/6 經 symlink 跳到 `Google_somatic_data/`，多一層失效點。

## L3 啟用前必驗 / 防護動作（建議，未執行）

- **驗甲基 tag**（5 個未確認 normal，每個跑一次）：
  ```bash
  samtools view /big8_disk/data/HCC1937/ONT/HCC1937BL.bam 2>/dev/null | head -2000 | grep -m1 -oE 'M[Mm]:Z:[^\t]+' || echo "NO_MM_TAG(無甲基)"
  ```
  有 `MM:Z:` → 可用於 V10；無 → 該樣本 G-A 需改用 tumor-only 路徑或排除。
- **偵測上游變更**：對 6 檔記 `stat -c '%s %Y'`（size+mtime）快照，定期比對；real-file 兩個可加 `samtools quickcheck` + 首尾 block md5。
- **對 PI / zhenyu112 知會**：這 6 檔是 paper 關鍵證據，建議勿動；或評估在 big7 留唯讀快照（HCC1395BL 149GB、COLO829BL 141GB → 容量需評估，big7 cleanup 後有空間）。

---

> **維護**：研究啟動跑 G-A 前**先讀本檔**，確認甲基 tag 已驗、SPOF 未變。新增其他外部依賴（reference genome / 第三方 callset）追加進本表。
