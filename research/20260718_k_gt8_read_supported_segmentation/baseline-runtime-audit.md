<!--
建立時間: 2026-07-18 Asia/Taipei
目標: 稽核舊 HCC1395 densest-8 canonical downstream pipeline 的歷史 runtime 來源，建立新 k>8 read-supported segmentation full run 的可追溯比較基線
處理範圍: 只讀檢查 terminal-success layered v3 run、HCC1395 sample timestamps、LongPhase-S upstream /usr/bin/time 紀錄、輸入 identity 與公平比較限制
關聯檔案:
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/pre-decision-audit.md
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/implementation-notes.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/
異動限制: 本文件為只讀稽核紀錄；未修改 full-run scripts、canonical outputs 或歷史 receipts
-->

# HCC1395 舊 densest-8 runtime 基線稽核

> **TL;DR：舊 HCC1395 densest-8 downstream 歷史基線為約 84.78 分鐘，但它是原始 run 目錄的 filesystem birth-timestamp proxy，不是單一 `/usr/bin/time`。高精度時間差為 5,086.484135 秒，即 84.774736 分鐘；79.20 / 4.11 / 0.39 / 1.08 分鐘分別對應五個 mlhp parts、layered reconstruction、region view、site ledger + manifest。LongPhase-S upstream 另有正式 `/usr/bin/time -v`：2:27:01、peak RSS 9,921,380 KiB、exit 0。**

## 1. 稽核結論與權威 run

歷史基線使用下列 terminal-success run：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/
```

選擇此 run，而不是同系列 v2、v3 或 v4，原因是：

- `state_events/000006_SUCCEEDED.json`：
  - `state = "SUCCEEDED"`
  - `reason = "verification and provenance readback passed"`
- `verification_summary.json`：
  - `all_pass = true`
  - `dataset_count = 7`
  - `biological_sample_count = 6`
- HCC1395 child worker 在 `children.json` 中為 `returncode = 0`。

因此，約 84.78 分鐘應標示為：

> **HCC1395 在 terminal-success v5 run 內，從 sample worker 第一個 mlhp log 建立到 sample `output_manifest.json` 建立的歷史 downstream wall proxy。**

它不是：

- 全部 7 datasets 的整個 run wall time；
- preflight + worker + verifier 的總時間；
- LongPhase-S upstream 的時間；
- 只針對 densest-8 選點演算法本身的 CPU 時間。

## 2. 原始路徑與 timestamp 邊界

HCC1395 sample 根目錄：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/HCC1395/
```

本稽核使用「輸出檔案 birth timestamp」作連續階段邊界：

| 符號 | 檔案 | Birth timestamp（Asia/Taipei） |
|---|---|---|
| `t0` | `mlhp_part_1.log` | `2026-07-13 01:39:24.833540934 +0800` |
| `t1` | `mlhp_part_5.json` | `2026-07-13 02:58:36.541667482 +0800` |
| `t2` | `layered_reconstruction_HCC1395.json` | `2026-07-13 03:02:42.965674045 +0800` |
| `t3` | `layered_region_view_HCC1395.json` | `2026-07-13 03:03:06.293674666 +0800` |
| `t4` | `output_manifest.json` | `2026-07-13 03:04:11.317676398 +0800` |

### 2.1 計算公式

```text
D_mlhp             = (t1 - t0) / 60
D_layered          = (t2 - t1) / 60
D_region           = (t3 - t2) / 60
D_ledger_manifest  = (t4 - t3) / 60
D_total            = (t4 - t0) / 60
                   = D_mlhp + D_layered + D_region + D_ledger_manifest
```

### 2.2 Exact values

| 階段 | Exact seconds | Exact minutes | 報告顯示值 |
|---|---:|---:|---:|
| 五個 mlhp parts，含 densest-8/read linkage 與 part handoff | `4751.708126548` | `79.1951354425` | `79.20 min` |
| Layered reconstruction | `246.424006563` | `4.1070667761` | `4.11 min` |
| Region view | `23.328000621` | `0.3888000104` | `0.39 min` |
| Site ledger + output manifest | `65.024001732` | `1.0837333622` | `1.08 min` |
| **HCC1395 downstream total** | **`5086.484135464`** | **`84.7747355911`** | **約 `84.78 min`** |

四個 exact stage seconds 的守恆：

```text
4751.708126548
+ 246.424006563
+  23.328000621
+  65.024001732
= 5086.484135464 seconds
```

### 2.3 為何歷史值寫 84.78，而 high-resolution total 是 84.774736

若使用 GNU `stat -c '%W'` 的整數秒 birth epoch：

```text
start = 1783877964
end   = 1783883051

(1783883051 - 1783877964) / 60
= 5087 / 60
= 84.783333... minutes
```

因此整數秒 proxy 顯示 `84.78 min`。使用納秒 birth timestamp 則是
`84.7747355911 min`。兩者只差約 `0.516 s`，對本輪分鐘級工程比較沒有實質影響，
但報告不可把 `84.78` 寫成正式 monotonic benchmark 的 exact value。

另外，四個 stage 各自四捨五入後：

```text
79.20 + 4.11 + 0.39 + 1.08 = 84.78 min
```

這也是歷史分項顯示值的來源。正式計算比值時應優先使用未四捨五入的
`5086.484135464 s`。

## 3. 舊 pipeline 執行參數與原始 command

### 3.1 Effective analysis parameters

來源：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/launch_receipt.json
```

與本比較直接相關的值：

| 參數 | 舊 v5 值 |
|---|---:|
| Scope | `chr1`–`chr22` |
| `MAX_SNV` | `8` |
| `MAPQ_MIN` | `20` |
| `BASEQ_MIN` | `0` |
| `MINREAD` | `3` |
| `TIER_R` | `50000` |
| `parallel_parts_per_sample` | `1` |
| `parallel_samples` | `4` |
| `ANALYSIS_TREE_CAP` | `0` |
| `DISPLAY_TREE_CAP` | `32` |

### 3.2 HCC1395 child worker exact argv

來源：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/children.json
```

實際 argv：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/source_bundle/files/core/runner.py \
  _worker \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --sample HCC1395
```

Registry evidence：

```text
label: worker:HCC1395
launched_at_utc: 2026-07-12T17:39:24.652711Z
command_sha256: c3091abc1f0cb44f4c969ce4feba66ac0145e52440c7647c5c57d5031925bbce
returncode: 0
```

同一 batch 還並行執行 `COLO829`、`H1437`、`H2009`；因此 HCC1395 的歷史
wall proxy 可能包含共享 CPU / NFS I/O contention。

### 3.3 Frozen source identities

來源：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/source_bundle/source_bundle_manifest.json
```

| 元件 | Frozen SHA-256 |
|---|---|
| `run_layered_v3.py` runner | `3f036b653c4580e1589f6d4ca60f2862f29575a2a0d3a69c6a6de52bd4583d5f` |
| `sm_multilocus_combinations.py` | `12ae8e9b79fecc7e66266cf39e13b10a81723dd954d1ff98c3ad22e0434e10bc` |
| `layered_tree_reconstruction.py` | `70d28fa4dfe3e69a0ea610e346457a894e7e024ffdb437895343fc66516639d0` |
| `build_region_view.py` | `cd39bb7799f6b62190f9dca2a3afad321d0b384e254ab3fb4a102e594d02d872` |
| `build_ssnv_site_ledger.py` | `cd6eaf955a13e250e18e4f07d10d557e7aa73c5c8d2f4c8aec84ff4f55171ad2` |

## 4. LongPhase-S upstream 正式 runtime

### 4.1 原始證據

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/longphase_s.production.log
```

原始 log 片段：

```text
total process time        : 8821s
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:27:01
Maximum resident set size (kbytes): 9921380
Exit status: 0
```

正式 `/usr/bin/time -v` 數值：

| 指標 | 值 |
|---|---:|
| Wall time | `2:27:01` = `8,821 s` = `147.0167 min` |
| User CPU | `26,720.31 s` |
| System CPU | `1,764.48 s` |
| CPU utilization | `322%` |
| Peak RSS | `9,921,380 KiB` |
| Exit status | `0` |

### 4.2 LongPhase-S exact argv

權威 command receipt：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/producer_capture_receipt_v2.json
```

```bash
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/source_bundle/longphase-s \
  somatic_haplotag \
  -s /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  -b /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam \
  --tumor-snv-file /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.clairs.raw_all.normalized.vcf.gz \
  --tumor-bam-file /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 12 \
  --tagSupplementary \
  -q 20 \
  --output-somatic-vcf \
  -o /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395_production
```

`2:27:01` 是 LongPhase-S producer 的正式 wall。producer 寫入 FIFO，因此會受同步
stream consumer backpressure 影響；但此時間不包含後續 sidecar gzip/index、VCF
recalibration 與 receipt closeout 的完整尾端。

## 5. 新舊輸入 identity 是否一致

舊 v5 frozen input：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/frozen_input_lock.json
```

新 k>8 full runner 使用的 M2 canonical manifest：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/release_contract/input_contract/canonical_manifest.json
```

HCC1395 關鍵 identity 對齊：

| 輸入 | 舊 v5 | 新 M2 | 結果 |
|---|---|---|---|
| Tumor BAM storage identity | `153c5bfaf2db93d06d4d9a9131e32d9bbdb58663985993f1339ae0387e4fe715` | 同左 | 一致 |
| LongPhase-S recalibrated PASS VCF | `bc393dda8eb549d46e9993eccebd4b7c4cd7bf809cde20e071b8fdd8ff6437e6` | 同左 | 一致 |
| Read-tag sidecar | `8e79e89cb668e65b22a8c02aee9dcb9fea81fc8903962c094839d317ebd93dfb` | 同左 | 一致 |
| Sidecar index | `f972c1538c91be36d504b516358faa86a44a36e2e452f0c0e0066e35d88b0d76` | 同左 | 一致 |

因此，若新舊都從已存在的 PASS VCF、tumor BAM 與 read-tag sidecar 開始，
輸入 provenance 足以支持 downstream operational runtime 比較。

## 6. 公平比較判定

### 6.1 可以做的比較

1. **Historical downstream workflow proxy**

   ```text
   runtime_ratio = new_full_wall_seconds / 5086.484135464
   ```

   可描述為「新 full HCC1395 k>8 segmentation 相對舊 densest-8 downstream
   historical wall proxy 的成本倍數」。

2. **最接近的 extraction-stage proxy**

   ```text
   extraction_ratio = new_extraction_wall_seconds / 4751.708126548
   ```

   但仍需標示 base-quality threshold、輸出內容與 concurrency 不同。

3. **新 DP partition 自身 runtime**

   新版有逐染色體 `/usr/bin/time -v`，應獨立報：

   - DP partition wall；
   - CPU；
   - peak RSS；
   - filesystem I/O。

   舊 densest-8 selection 沒有獨立 instrumentation，因此不可計算嚴格的
   `DP speedup over densest-8 selection`。

### 6.2 不可直接宣稱嚴格 speedup 的原因

| Caveat | 舊 v5 | 新 k>8 full | 影響 |
|---|---|---|---|
| 計時工具 | filesystem timestamp proxy | `/usr/bin/time -v` + runner receipts | 精度與指標不同 |
| Sample concurrency | `parallel_samples=4` | HCC1395 sequential per chromosome | CPU / I/O contention 不同 |
| Base quality | `BASEQ_MIN=0` | `BASEQ_MIN=20` | read-call workload 與保留量不同 |
| Threading | 舊 Python worker搭配共享 batch | `samtools_threads=1` | 不可直接作 per-core speed claim |
| 分析範圍 | 全部 8,222 multilocus regions | 全 PASS-sSNV sparse extraction；partition 聚焦 408 個 k>8 components | 工作內容不同 |
| Downstream 終點 | mlhp + tree + region + ledger | extraction + partition | 舊基線包含較多後處理 |
| 舊資源 receipt | 無 per-sample CPU/RSS/I/O | 有逐 stage CPU/RSS/I/O | 只能比較 wall proxy |
| Host load | preflight 時低負載；runtime 全程未監控 | 本輪另有 resource gate/實際 contention | 時代與負載不同 |

### 6.3 Preflight 資源不能代表 runtime 全程

舊 v5：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/process_observation.json
```

Preflight observed：

```text
logical_cpus: 48
load_1m: 1.1953125
load_per_cpu: 0.02490234375
iowait_percent: 2.0010840711
nfs_read_mbps_decimal: 0.0
available_memory_bytes: 516185526272
```

這只是 worker 啟動前的五分鐘 observation；不能推論 HCC1395 的 84.78 分鐘期間
一直維持相同 load 或 NFS traffic。

### 6.4 LongPhase-S upstream 的納入規則

舊 84.78 分與新 full runner 都使用已存在的 LongPhase-S PASS VCF 與 read-tag
sidecar。因此目前 downstream 比較應：

- 舊、新兩邊都**排除** LongPhase-S `2:27:01`；
- 若改成 raw-to-result operational comparison，兩邊都加相同 upstream；
- 不可只替其中一邊加入 `2:27:01`。

舊 producer wall 加舊 downstream proxy，僅作近似參考：

```text
8821.000000 + 5086.484135
= 13907.484135 seconds
= 231.791402 minutes
≈ 3 h 51 min 47 s
```

此加總仍不含 producer 後續所有 gzip/index/recalibration/receipt closeout 尾端，
不可稱完整 raw-to-final exact wall。

## 7. 本次只讀稽核命令

### 7.1 Timestamp inventory

輸入：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/HCC1395/
```

命令：

```bash
base=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/HCC1395

for f in \
  mlhp_part_1.log mlhp_part_1.json \
  mlhp_part_2.log mlhp_part_2.json \
  mlhp_part_3.log mlhp_part_3.json \
  mlhp_part_4.log mlhp_part_4.json \
  mlhp_part_5.log mlhp_part_5.json \
  layered.log layered_reconstruction_HCC1395.json \
  region_view.log layered_region_view_HCC1395.json \
  site_ledger.log ssnv_site_ledger_HCC1395.tsv.gz \
  ssnv_site_ledger_HCC1395.tsv.gz.tbi \
  ssnv_site_ledger_HCC1395.summary.json output_manifest.json
do
  stat --printf='%n\tbirth=%w\tmtime=%y\tctime=%z\n' "$base/$f"
done
```

實際輸出關鍵片段：

```text
mlhp_part_1.log
birth=2026-07-13 01:39:24.833540934 +0800

mlhp_part_5.json
birth=2026-07-13 02:58:36.541667482 +0800

layered_reconstruction_HCC1395.json
birth=2026-07-13 03:02:42.965674045 +0800

layered_region_view_HCC1395.json
birth=2026-07-13 03:03:06.293674666 +0800

output_manifest.json
birth=2026-07-13 03:04:11.317676398 +0800
```

此命令沒有分析輸出檔；結果只顯示於 terminal，並整理至本文件。

### 7.2 Terminal-state verification

命令：

```bash
python - <<'PY'
import json
from pathlib import Path

root = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
)

terminal = json.loads((root / "state_events/000006_SUCCEEDED.json").read_text())
verification = json.loads((root / "verification_summary.json").read_text())

print(terminal["state"], terminal["reason"])
print(
    verification["all_pass"],
    verification["dataset_count"],
    verification["biological_sample_count"],
)
PY
```

實際輸出：

```text
SUCCEEDED verification and provenance readback passed
True 7 6
```

### 7.3 LongPhase-S resource lines

命令：

```bash
rg -n \
  'total process time|Command being timed|Elapsed \(wall clock\)|Maximum resident|Exit status' \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/longphase_s.production.log
```

實際輸出位置：

```text
line 79:  total process time : 8821s
line 100: Command being timed: longphase-s somatic_haplotag ...
line 104: Elapsed (wall clock) time: 2:27:01
line 109: Maximum resident set size (kbytes): 9921380
line 122: Exit status: 0
```

## 8. 報告使用規則

後續 HTML 或數據報告建議固定使用以下措辭：

> 舊 HCC1395 densest-8 downstream historical wall proxy 約 84.78 分鐘
>（terminal-success v5；filesystem birth timestamps；不是正式 monotonic benchmark）。
> 新 full run 使用相同 BAM、LongPhase-S PASS VCF 與 read-tag sidecar identity，
> 因此可比較 operational wall cost；但因 base-quality threshold、sample concurrency、
> instrumentation 與 downstream endpoint 不同，不宣稱嚴格 algorithmic speedup。

若需要論文級效能結論，必須另跑：

1. 同一時間窗；
2. 相同 host load gate；
3. 相同 `MAPQ_MIN` / `BASEQ_MIN`；
4. 相同 sample concurrency；
5. 相同輸入與輸出 endpoint；
6. 兩個方法都由 `/usr/bin/time -v` 包住；
7. 至少三次重複，報 median、range 與 peak RSS。

## 9. 新 full run 的 shared-host operational 環境紀錄

這不是舊 baseline 的一部分，只用來界定本輪新 wall time 的可解讀範圍。

- 2026-07-18 06:36 +08:00 啟動前：48 logical CPUs，load1=25.2，通過預註冊
  `load1≤60` gate。
- 執行政策：chr1–22 sequential；外層 `nice -n 10 ionice -c2 -n7`。
- 2026-07-18 07:31 +08:00 單次 spot check：load average
  `67.42, 65.33, 55.22`；同一個一秒 `iostat -xz` 視窗的 `sdc`
  `%util=100.00`。

實際命令：

```bash
uptime
awk '/MemAvailable/ {print}' /proc/meminfo
df -h /big7_disk | tail -1
iostat -xz 1 2
```

實際輸出關鍵片段：

```text
07:31:38 up 13 days, 21:14, 21 users, load average: 67.42, 65.33, 55.22
MemAvailable: 391346948 kB
/dev/sdc1 42T 39T 716G 99% /big7_disk
sdc ... aqu-sz=15.37 %util=100.00
```

因此，新 full run 的正式 `/usr/bin/time -v` 是在 shared host、低優先序與可觀察
I/O contention 下的實際 operational wall。它適合回答「這次使用者實際等多久」，
不適合單獨證明演算法 speedup。
