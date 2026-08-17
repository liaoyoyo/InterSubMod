<!--
建立時間: 2026-08-13 09:54 +08:00
目標: 在把 frozen 7-dataset BAM/HP-PS manifest 接入一頁式 dashboard 前，審核身分、分母、漂移與 claim ceiling
處理範圍: 7 datasets / 6 biological samples；只做 bounded metadata、sampled chunks、BAI、header、quickcheck 與既有 validation receipt readback；不掃完整 BAM、不重跑 LongPhase-S/ISM
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_metric_contract.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/implementation-notes.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json
-->

# Pre-decision audit：多 BAM manifest intake 與 producer receipt panel

> **Verdict：GO_WITH_BOUNDED_IDENTITY（80/100）— 可以建立 manifest-driven snapshot 與工程 QC panel；production/live、full-BAM identity、primary-read HP/PS QC 與科學驗證仍是 NO-GO。**

## §0 Cynefin domain gate

- **Domain：Complicated with one portability edge case。** 既有 `storage_identity_v1`、BAI full SHA、sidecar validation 與 dashboard artifact 都有可重複方法；但 mount remap 使 `st_dev` 漂移，必須把 strict metadata 與 payload evidence 分欄。
- **本輪安全邊界：** 只讀 7 BAM 的 stat/header/EOF、每 BAM 三個 1 MiB chunks、BAI 與小型 receipts；不做 flagstat、mosdepth、read N50 或 MM/ML 全量掃描。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| Frozen source manifest 宣告 7 datasets / 6 biological samples，sample set 與 topology CSV 一致 | ✓ | L1 | `input_manifest.snapshot.json` + `results/cohort_topology_metrics.csv` |
| 7/7 BAM 有 sampled-chunk identity，7/7 BAI 有 full SHA256 | ✓ | L1 | source manifest `samples[].alignment_payload.storage_identity_v1` |
| HCC1395 quick pilot：三 chunks、BAI SHA、validation receipt SHA 均吻合 | ✓ | L1 | 2026-08-13 shell readback |
| 7/7 current payload sampled chunks + BAI identity 與 frozen values 相符 | ✓ | L1 | 2026-08-13 bounded all-sample readback |
| 7/7 strict storage identity 完全相符 | ✗ | L1 | frozen `st_dev=59`，current `st_dev=60`；其他 payload fields 相符 |
| 7/7 sidecar native validation 為 all-region PASS、missing/extra/conflict/unknown HP=0 | ✓ | L1 | `samples/*/sidecar_validation.json` |
| BAM full-content SHA256、reference FASTA/hash、primary-mapped denominator | □ | — | 目前 manifest/receipt 未提供 |
| depth、read N50、MM/ML completeness、KDE、truth F1 | □ | — | 本輪禁止以 0 或 proxy 補值 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | manifest → bounded identity receipt → dashboard 是可稽核資料產品契約 |
| 觀察支撐 | 20 | 7/7 direct manifest/receipt；HCC pilot 與全 7 bounded payload readback |
| 機制清晰度 | 20 | strict identity 失敗只由 `st_dev` remount drift 解釋，chunks/size/mtime/inode/index均吻合 |
| 反例風險 | 10 | sampled chunks 不是 full hash；producer all-alignment denominator不等於 primary-read B5 |
| 所需資源 | 10 | schema、normalizer、artifact wiring、browser QA 約 1–6 小時 |
| **TOTAL** | **80 / 100** | **GO_WITH_BOUNDED_IDENTITY** |

**Falsifier observable：** 任一 dataset 的 sample set、chunk SHA、BAI SHA、validation receipt SHA、quickcheck 或 sidecar conservation 失敗，就不得標 `PAYLOAD_SAMPLE_MATCH`／`PRODUCER_RECEIPT_PASS`，dashboard build 必須 fail closed 或顯示 FAIL；不得由其餘 6 份插補。

**Reality-test 三反例：**

1. `st_dev` 因 remount 改變時，不能誤寫成 payload corrupted；同時也不能把 strict identity 寫成 PASS。
2. `total_tagged/total_alignment` 的 denominator 包含 primary、secondary、supplementary alignment records；不得標成 canonical B5 primary-read HP/PS rate。
3. BAM 路徑存在與 quickcheck PASS 不代表 reference、depth、MM/ML、KDE 或 truth benchmark ready。

## §3 Assumption map

| Assumption | Importance | Known | Action |
|---|---|---|---|
| frozen sample set 與 topology dataset identity 完全一致 | HIGH | KNOWN | builder 逐列 assert |
| `st_dev` 是 payload identity 必要條件 | HIGH | KNOWN FALSE for portability | 分開 strict metadata / payload sample status；不改寫 upstream manifest |
| sampled chunks 可等同 BAM full SHA | HIGH | KNOWN FALSE | `full_content_sha256=null`、claim 明示 bounded |
| producer HP/PS receipt可支持 all-alignment tag coverage | HIGH | KNOWN | 顯示 exact numerator/denominator與population |
| producer HP/PS receipt可支持 primary-read B5 | HIGH | KNOWN FALSE | B5 保持 `NOT_COLLECTED` |
| 全量 QC 可在本輪無成本完成 | LOW | KNOWN FALSE | 長計算另立 job/receipt cycle |

## §4 Quick pilot（已執行）

1. HCC1395 `samtools quickcheck -v` → exit 0、stdout/stderr 空。
2. 讀 BAM first/middle/last 各 1 MiB → SHA256 全部吻合 frozen manifest；BAI full SHA 與 validation receipt SHA 亦吻合。
3. 對 7 datasets 重算 `storage_identity_v1` → payload chunks/index/size/mtime/inode 7/7 match；strict identity 0/7，唯一差異欄位為 `st_dev` 與衍生 `identity_sha256`。

**Checkpoint：PASS for bounded payload intake；FAIL for strict metadata identity。** 因此 UI 必須同時顯示 `7/7 payload sample match` 與 `0/7 strict metadata match (mount-device drift)`。

## §5 Gap diagnosis

| Missing / gap | Impact | Effort | Priority |
|---|---:|---:|---:|
| Manifest schema + normalizer + validator | HIGH | 1–2 h | P0 |
| strict metadata vs sampled payload split | HIGH | <1 h | P0 |
| dashboard opportunity/read-tag panels + selector wiring | HIGH | 1–2 h | P0 |
| reference FASTA/hash + RG/reference compatibility contract | HIGH | 1–3 h | P1 |
| primary-read HP/PS denominator、mapping/depth/MM-ML/KDE/truth receipts | HIGH | >6 h | P2（另立長計算 cycle） |

## §6 Evidence conflict scan

- Repository `InterSubMod/MEMORY.md` 不存在，明示 unavailable；未以外部記憶補寫。
- `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` 是甲基假說 NEGATIVE，與本輪工程 intake 無直接衝突，但限制任何 methylation benefit claim。
- `InterSubMod/docs/CURRENT_FOCUS.md` 明示 latest HP/PS 應來自 same-run sidecar、BAM embedded HP 不作權威；本輪沿用該 producer contract。
- 舊 dashboard contract 將 B1/B5 標 `NOT_COLLECTED`；新直接 receipt 可解除「input manifest 完全未知」與「producer all-alignment tag status 未知」，但不能解除 full BAM SHA、reference 或 primary-read B5 gate。

## §7 Decision path與紅隊

- **最強反方 failure mode 1：** 把 7 月 frozen manifest 當 current strict PASS，遮蔽 `st_dev` drift。防線：雙欄狀態 + 7/7 current readback receipt。
- **最強反方 failure mode 2：** 將 all-alignment tag rate誤稱 primary read QC。防線：dataset 名稱、圖標題、tooltip、method/footer固定寫 denominator。
- **最強反方 failure mode 3：** 看到 BAM present 後把 mapping/depth/MM-ML/truth 缺值變成 0。防線：schema `null + NOT_COLLECTED`；dashboard percent validator拒絕假 0 proxy。
- **Concluded overlap：** 無重新開啟 methylation filter NEGATIVE；本輪只呈現工程 provenance/availability。
- **Verdict：GO_WITH_BOUNDED_IDENTITY。** Red-team 通過，因每個 failure mode都有 fail-closed observable。
- **Decision lock：Y。** 不覆寫 upstream manifest/BAM/receipts；本輪僅新增 normalized snapshot、schema、script、dashboard artifact/HTML與驗證文件。
