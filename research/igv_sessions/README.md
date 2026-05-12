<!--
建立時間: 2026-05-01 00:00
最新更新: 2026-05-12（v6_germline_absent_audit.xml 因 sed-replace bug archive；改用 v5_v6_compare_with_paired.xml 並列）
目標: 說明 research/igv_sessions/ 目錄內所有 IGV session（XML）用途與管理規範
處理範圍: IGV session、phase-context VCF、TP/FP VCF、audit marker BED 與 site-level manifest
關聯檔案:
  - research/igv_sessions/v5_purity_compare_with_paired.xml
  - research/igv_sessions/v5_v6_compare_with_paired.xml
  - research/igv_sessions/_archive/v6_germline_absent_audit_BROKEN_replaced_v5.xml
  - research/igv_sessions/annotations/site_layer_manifest.tsv
-->

# IGV Sessions 管理規範

## 目錄定位

**固定路徑**：`InterSubMod/research/igv_sessions/`

所有跨 baseline / V3F / V5 / V6 / paired 的 IGV audit session XML 都收在此目錄。**不要分散到其他位置** — 統一管理 audit infrastructure（phase VCF / TP/FP marker / LOH/GE BED / colored audit sites）。

## 命名規範

`{version}_{purpose}_{scope}.xml`

範例：
- `v5_all.xml` — V5 全 sample 基礎 session
- `v5_all_4versions.xml` — V5 baseline/V2b/V3F/V5 四版並列
- `v5_color_HP.xml` — V5 HP tag 著色版
- `v5_purity_compare_with_paired.xml` — V5 + purity sim + paired ground truth
- `v5_v6_compare_with_paired.xml` — **新增** V5 + V6 並列 + purity sim + paired（V6 germline-absent revert vs V5 audit）
- `v6_all_5versions.xml` — **新增（2026-05-12）** 以 `v5_all_4versions.xml` 為基底**加入** V6 panel；baseline/V2b/V3F/V5/V6 五版並列（乾淨版，無 purity sim BAM 雜訊）
- `_archive/v6_germline_absent_audit_BROKEN_replaced_v5.xml` — **已棄用**：前一個 agent 用 `sed s|V5|V6|g` 整檔取代 V5 路徑為 V6，導致 V5 track 消失、V6 BAM 重複 load；改由 `v5_v6_compare_with_paired.xml` 取代

## Session 清單（依時序）

| Session XML | 用途 | BAM tracks | VCF context | BED markers |
|---|---|---|---|---|
| `v5_all.xml` | 基礎 V5 + baseline | 2 | — | — |
| `v5_all_4versions.xml` | 四版（baseline/V2b/V3F/V5）並列 | 4 | — | — |
| `v5_all_4versions_color_MOD.xml` | 上加 MOD 著色 | 4 | — | — |
| `v5_all_color_MOD.xml` | V5 + MOD 著色 | 2 | — | — |
| `v5_color_HP.xml` | V5 HP tag 著色 | 2 | — | — |
| `v5_purity_compare.xml` | V5 + purity simulation | 4 | — | — |
| `v5_purity_compare_with_paired.xml` | 上加 paired ground truth + 全 audit context | 6 | 7 | 6 |
| **`v5_v6_compare_with_paired.xml`** | **V5 + V6 並列 + purity sim + paired，audit germline-absent revert 對 V5 的差異** | 7 | 7 | 6 |
| **`v6_all_5versions.xml`** | **5 版並列乾淨版（baseline/V2b/V3F/V5/V6 + tumor untagged + normal），給 PPT slide 04a/04b/16 使用** | 7 | 3 | 2 |
| `_archive/v6_germline_absent_audit_BROKEN_replaced_v5.xml` | **棄用**（sed-replace bug，V5 path 被整檔覆寫成 V6） | 6 | 7 | 6 |

## 全 session 共同 audit layer（已建好的 annotation 資源）

詳見 `annotations/` 子目錄：

- `audit_sites_all_colored.bed` — 全特殊位點合併（含 TP/FP/V5max/SelfPhasing 顏色）
- `audit_sites_TP.bed` / `audit_sites_FP.bed` — paired 層級 TP/FP 分類
- `audit_sites_V5max.bed` — V5 改變最大位點
- `audit_sites_SelfPhasing.bed` — self-phasing 觀察位點（含 SP1/SP2/SP3）
- `site_layer_manifest.tsv` — 每位點 TP/FP VCF 命中 + phase anchor 摘要
- `lp_s_normal_phase_context_15sites.vcf` / `lp_to_*_phase_context_15sites.vcf` — paired/TO 各版 phase backbone

## V6 (V5+V6 並列) Session 啟用流程

> **⚠ 教訓（2026-05-12）**：第一版 `v6_germline_absent_audit.xml` 用 `sed s|V5|V6|g` 整檔替換，導致 V5 path 完全消失、V6 BAM 重複 load；該檔已 archive 到 `_archive/`。新版改為**加入** V6 panel 而非取代 V5（V5 與 V6 並列）。**未來建立 V7+ session 時，必須複製整個 V5 panel block 並只修改新 panel 內的 path/attributeKey/name；不可用全檔 sed-replace**。

1. **建立**：複製最近並列 audit session（推薦 `v5_v6_compare_with_paired.xml`）→ `v{N-1}_v{N}_compare_with_paired.xml`；在 `<Resources>` 與 BAM Panel section **加入** V{N} entry（不替換現有）
2. **載入**：開 IGV → File → Open Session → 選新 XML
3. **batch snapshot**（headless）：IGV 用 `--batch` 模式 + script 指定位點 + snapshot
4. **產出**：PNG 入 `by_HP_v{N-1}v{N}/` 目錄（區分純 V{N-1} 與並列）
5. **記錄**：在本 README 加 entry + 在產出目錄加 `_session_note.md` 連結回此 session
6. **驗證**：`xmllint --noout` + `grep -c <V{N-1}_path>` 與 `grep -c <V{N}_path>` 各應 ≥3（Resource + Coverage + reads）

## 後續 version 規範（V7 / V8 ...）

未來新版本：
- 複製最近 audit-complete session 為基礎
- 命名 `v{N}_{purpose}_audit.xml`
- 在本 README 加 entry
- 不刪除舊 V{N-1} session（保留 reproducibility）

---

# IGV V5 Purity Compare With Paired Session

## 目的

此 session 用於人工確認 HCC1395 LongPhase-TO baseline / V5 在不同 purity 條件下的 haplotag 結果，並與 paired LongPhase-S ground truth 對照。新增的 track 目標是同時檢查三件事：

1. phase graph / phase VCF 在特殊位點附近是否有 anchor 與 `GT|PS` 資訊。
2. BAM read 的 `HP` tag 是否與 paired truth、baseline、V5 一致。
3. 每個特殊位點在 LongPhase-S / LongPhase-TO 的 TP、FP 層級判斷是否一致。

## Track 分層

| Panel | Track | 用途 |
|---|---|---|
| `VariantPanel` | ClairS-TO@0.93 / ClairS-TO@0.6 | 原始 tumor-only candidate VCF |
| `PhaseContextPanel` | LP-S normal phase context | paired normal phase anchor 參考 |
| `PhaseContextPanel` | LP-TO BL/V5 phase context | tumor-only baseline / V5 phase anchor 對照 |
| `TPFPVariantPanel` | LP-S paired TP/FP VCF | paired truth 層級的 TP/FP 參考 |
| `TPFPVariantPanel` | LP-TO PASS TP/FP VCF | tumor-only 層級的 TP/FP 結果 |
| `AuditSiteMarkerPanel` | all colored / TP / FP / V5max / self-phasing markers | 特殊位點定位與分類顏色 |
| `BedPanel` | LOH / GE BED | phase block 與 genome event 背景 |
| BAM panels | paired、baseline、V5、purity 0.6 BAM | 用 `HP` tag 檢查 read-level tag 結果 |

## 顏色規則

| 分類 | 顏色 | RGB |
|---|---|---|
| TP | 藍色 | `37,99,235` |
| FP | 紅色 | `220,38,38` |
| V5max | 紫色 | `124,58,237` |
| Self-phasing | 橘色 | `234,88,12` |

## 判讀順序

1. 先看 `Audit markers: all colored sites` 定位特殊位點與分類。
2. 看 `TPFPVariantPanel`：確認該位點是否出現在 LP-S paired TP/FP 或 LP-TO TP/FP VCF。
3. 看 `PhaseContextPanel`：確認特殊位點附近是否有 phased anchor、`GT` 是否為 phased genotype、`PS` 是否落在合理 phase set。
4. 看 BAM panels：以 `colorByTag=HP` 與 `groupByOption=PHASE` 檢查 paired、baseline、V5 的 read-level haplotag 是否一致。
5. 交叉看 `BedPanel` 的 LOH/GE track，避免把 LOH 或 genome event 背景誤判為 tagging 改善。

## 重要限制

- `LP-TO V5@0.93 phase context` 使用 `pononly_v2b/tumor_phased.vcf` 作為 phase backbone；目前 V5 somatic fallback run 主要產出 tagged BAM，未看到獨立 V5 phase VCF。
- `site_layer_manifest.tsv` 的 `*_phase_n1kb` 是特殊位點正負 1 kb 內 phase-context VCF record 數，代表可用於人工審查的 phase anchor 密度，不等於該位點一定被當作 anchor。
- V5max 與 self-phasing marker 主要是 tag 行為觀察點，部分 marker 沒有 ref/alt 欄位，因此不會在 TP/FP VCF exact-match 欄位標成 `Y`。

## 產物

| 檔案 | 說明 |
|---|---|
| `v5_purity_compare_with_paired.xml` | IGV session 主檔 |
| `annotations/audit_sites_all_colored.bed` | 全部特殊位點合併 colored BED |
| `annotations/audit_sites_TP.bed` | TP 特殊位點 |
| `annotations/audit_sites_FP.bed` | FP 特殊位點 |
| `annotations/audit_sites_V5max.bed` | V5 改變最大位點 |
| `annotations/audit_sites_SelfPhasing.bed` | self-phasing 觀察位點 |
| `annotations/site_layer_manifest.tsv` | 每個特殊位點的 TP/FP VCF 命中與 phase-context anchor 摘要 |
