<!--
建立時間: 2026-07-11 04:12 +08:00
目標: 以唯讀方式稽核 7-dataset clean LongPhase-S FIFO HP/PS sidecar 與 layered v2 full run 的 process exclusivity、canonical input integrity、CPU/RAM/disk/inode/I/O 容量及 consumer readiness
處理範圍: host=bip7；觀察窗 2026-07-11 04:01-04:12 +08:00；production manifest、既有 production run root、layered manifest、既有 layered resource baseline；不啟動長計算、不終止或修改 process
關聯檔案:
  - InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_longphase_production_sidecars.sh
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh
task_type: B Comprehensive validation / Gate-1 launch-readiness
goals: G4,G5
host: bip7
observation_cutoff: 2026-07-11T04:11:58+08:00
evidence_level: host/process/file observations O-L1; source/manifest facts F-L1; resource projections I-L2
mutation_policy: read-only; no kill, signal, process mutation, or long computation
-->

# Launch Resource Red Team：clean tagging + layered v2

> **TL;DR — 目前對任何「額外啟動」判定 `NO-GO`（45/100）：同一個 7-dataset production FIFO-sidecar run 已在 bip7 執行，兩個 LongPhase-S process 已把 big8 NFS 讀取推到約 116 MB/s，host 短窗 iowait 57–60%；且截至截點只有 2 個 `START`、0/7 sidecar 完成，active layered manifest 仍沒有 sidecar artifact 路徑。** RAM、inode 與 FIFO-sidecar預估容量本身足夠；禁止啟動重複 job，也不可在 production run 尚未完成時疊加 layered full run。（影響：高；信心：高）

用 **Pre-Mortem + Assumption/Falsifier Map**：先分「既有 job 是否可被干擾」與「是否可再 launch」，再以 observable threshold 決定 `NO-GO → PROBE → GO`。本判定只服務 resource/launch safety，不取代方法學與 biological claim gate。

## 1. 研究問題、假設與邊界

| 欄位 | 本次定義 |
|---|---|
| 研究問題 | bip7 是否可安全承載 7 datasets 的無 truth flags LongPhase-S FIFO HP/PS sidecar，之後再跑 sidecar-aware layered v2 full run？ |
| 假設 | `SM_PARALLEL_SAMPLES=2` 不會超出 RAM/CPU；FIFO 避免落地完整 tagged BAM；canonical inputs 在 run 期間唯讀；sidecar consumer 可 exact-match alignment identity。 |
| 成功條件 | 無重複/conflicting process；7/7 production validation PASS；input/code identity readback；sidecar projected output ≤60 GiB；bounded consumer probe 的 missing/conflict=0、exact=exposures；之後才可 full run。 |
| 失敗條件 | 同 pipeline 已存在卻再 launch；canonical input 出現 writer/identity drift；FIFO 變 regular BAM；output projection >60 GiB 或 big7 free <500 GiB；consumer fallback 回 BAM HP/PS；任一 sidecar missing/conflict。 |
| 主要指標 | process/PID/command、load/iowait/NFS kB/s、RSS、disk free/inode、input stat digest、sidecar row projection、historical wall/RSS/output baseline。 |
| 主要限制 | point-in-time observation，沒有 filesystem lease；clean-tagging 尚無 completed sample 的實測 `/usr/bin/time -v`；resource estimate 不等於結果正確。 |
| 本次改動 | 只新增本稽核文件；沒有修改 pipeline、manifest、run root 或 process。 |
| 預計輸出 | 本檔 `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/06_launch_resource_redteam.md`。 |

研究啟動 5 問：Thread D=是（HP/PS read-level evidence）；Thread B=否；KDE=不適用；caller AF=不進 resource 判定；長計算/檔案 IO=是，因此本輪只稽核、不啟動。

## 2. Pre-decision verdict

### 2.1 分數與 hard gates

| 維度 | 分數 | 直接觀察 |
|---|---:|---|
| Process exclusivity | 0/20 | 相同 production wrapper 與 2 個 sample jobs 已執行；再 launch 必定重複。 |
| Canonical input stability | 16/20 | 29 個 manifest input 的 follow-symlink stat digest 兩次相同；`lsof` 只見 `r` handles，未見 writer；但沒有 immutable lease/full BAM checksum。 |
| Disk/inode strategy | 16/20 | FIFO sidecar預估可容納；但 big7 已 98%，若誤落完整 BAM 則 1.674 TiB > 932 GiB free。 |
| CPU/RAM/I/O | 8/20 | RAM available 471 GiB；兩 process RSS 合計約 13.6 GiB；但 NFS 約 116 MB/s、host iowait 57–60%、load 70.58/48 logical CPUs。 |
| Consumer readiness | 5/20 | code 已有 exact sidecar key/gate；active layered manifest 仍是 schema 2.0 且 7/7 `read_tag_sidecar=null`，目前必定 preflight fail。 |
| **總分** | **45/100** | **Hard NO-GO for any additional launch**。 |

### 2.2 判定

- **`NO-GO`（現在）**：不可再啟動 production sidecar；不可同時啟動 layered full run。
- **既有 run 處置**：本稽核不 kill、不 signal、不改 process；只把它視為「已在外部進行的 run」，其成功與否由 7/7 artifacts/validation 決定。
- **`PROBE`（下一個可接受動作）**：既有 production run 完成且 7/7 validation/hash readback 後，只跑一個 bounded sidecar consumer probe；不得直接把舊 7/7 layered PASS 當 clean-tagged evidence。
- **Conditional `GO`**：bounded probe 通過 §7 thresholds，且 launch 前 5-minute resource baseline 低於門檻，才可啟動 sidecar-aware 7-dataset layered full run。

**Verdict falsifier**：若在新 observation cutoff 同時看到 (1) 無任何 LongPhase/layered同類 process、(2) production sidecars 7/7 PASS且 code/input hashes readback OK、(3) schema 2.1 layered manifest 7/7 路徑完整、(4) COLO829 bounded consumer probe 全 exact、(5) iowait<20%且 big7 free≥500 GiB，則本次 `NO-GO` 可被推翻為 `GO`。

## 3. Host 與 process observation

### 3.1 Host/time/capacity

**輸入**：bip7 host state。  
**命令**：

```bash
date --iso-8601=seconds
hostname -f
uptime
lscpu | rg '^(Architecture|CPU\(s\)|On-line CPU|Model name|Thread\(s\)|Core\(s\)|Socket\(s\)|NUMA node\(s\))'
free -h
df -hT /big7_disk/liaoyoyo2001/big7_disk_output /big8_disk/data
df -ih /big7_disk/liaoyoyo2001/big7_disk_output /big8_disk/data
```

**實際輸出片段（2026-07-11T04:11:58+08:00）**：

```text
host=bip7
load average: 70.58, 67.89, 50.20
CPU(s)=48; AMD Opteron 6344; 24 physical cores / 48 threads
Mem=503 GiB total; 471 GiB available; swap 891 MiB used
big7=/dev/sdc1 42T total, 932G available, 98% used
big8=NFS 27T total, 845G available, 97% used
big7 inode=481M free (29% used); big8 inode=410M free (3% used)
```

RAM/inode不是 blocker；I/O contention與錯誤輸出模式才是 blocker。

### 3.2 已存在的相同 pipeline

**命令（host PID namespace，唯讀）**：

```bash
ps -eo pid,ppid,user,lstart,etime,stat,psr,%cpu,%mem,rss,args |
  rg -i 'run_longphase_production_sidecars|longphase-s somatic_haplotag|capture_longphase_tagged_bam_sidecar|run_layered_7samples_newbb|sm_multilocus_combinations|layered_tree_reconstruction|verify_layered_v2' |
  rg -v 'rg -i'
```

**實際輸出摘要（04:02:30 +08:00）**：

| PID | 角色 | Sample / command contract | 狀態 |
|---:|---|---|---|
| 656363 | wrapper | `run_longphase_production_sidecars.sh` | active |
| 656465 | LongPhase-S | HCC1395_DORADO, `-t 12 -q 20 --tagSupplementary`, 無 truth flags | active, ~100% CPU, RSS ~6.6 GiB |
| 656470 | LongPhase-S | HCC1395, 同上 | active, ~100% CPU, RSS ~7.0 GiB |
| 656462/656467 | FIFO capture | `capture_longphase_tagged_bam_sidecar.py` | blocked/等待 tagged BAM stream |

Run root：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1/
```

截至 04:11:58，`run_status.tsv` 只有：

```text
HCC1395_DORADO  production_tagging START expected_records=112387;truth_flags=absent
HCC1395         production_tagging START expected_records=113997;truth_flags=absent
```

未觀察到另一個 active layered runner；但 production run 本身就是預定工作，故任何第二次 production launch 都是明確衝突。

### 3.3 Current per-process telemetry

**命令**：

```bash
pidstat -dru -p 656465,656470 1 3
mpstat 1 3
vmstat 1 3
nfsiostat 1 3
```

**實際輸出摘要**：

```text
LongPhase-S HCC1395_DORADO: avg CPU 103.67%, RSS 6,987,428 KiB, read 57,685 kB/s
LongPhase-S HCC1395       : avg CPU 117.33%, RSS 6,638,439 KiB, read 56,976 kB/s
big8 NFS short interval   : ~116,486-116,871 kB/s read, avg RTT ~29 ms, avg exe ~53 ms
host aggregate short span : iowait ~57-60%; 64-65 blocked tasks
```

兩個 sample 已大致吃滿觀察到的 big8 NFS read throughput；把 parallel 從 2 增至 3/4 不會有可信加速依據，反而提高 tail latency 與 failure risk。

## 4. Canonical input writer audit

Canonical production input由：

```text
InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json
```

列出 7 tumor BAM、7 normal BAM、7 caller PASS VCF、7 phased germline VCF與 1 reference（symlink 以 `stat -L` 追到 target）。

**Open-file mode 命令**：

```bash
jq -r '.samples[] | [.tumor_bam,.normal_bam,.caller_pass_vcf,.germline_phased_vcf,.reference] | .[]' \
  InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json |
  sort -u | xargs -r lsof --
```

**觀察**：只有 PIDs 656465/656470 對兩個 normal BAM 的 `3r..14r` read FDs；未見 `w`/`u` writer。其餘 input 在該瞬間沒有 open handle。

**Identity snapshot 命令**：

```bash
jq -r '.samples[] | [.tumor_bam,.normal_bam,.caller_pass_vcf,.germline_phased_vcf,.reference] | .[]' \
  InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json |
  sort -u | xargs -r stat -Lc '%d:%i\t%s\t%Y\t%n' | sha256sum
```

兩次 snapshot 相同，digest：

```text
68151418579ab1273faef4da76dad53e0b7cb92b1370d3ebcb07b64ebbd633b8
```

**結論**：`[O-L1]` 觀察窗內未見 canonical input writer或 size/mtime/inode drift。這是 point-in-time evidence，不是鎖；BAM內容沒有 full checksum，run只 hash BAI與小型VCF/manifest，因此完成時仍須做 storage identity/readback。

## 5. Resource projection

### 5.1 Clean tagging：FIFO HP/PS sidecar，非 persisted BAM

Manifest input BAM（tumor+normal）合計：

```text
2,732,404,335,632 bytes = 2.485 TiB
```

7 個 tumor BAM 由 `samtools idxstats` 得到 mapped alignments：

| Dataset | mapped alignments |
|---|---:|
| HCC1395 | 42,516,868 |
| HCC1395_DORADO | 41,498,156 |
| COLO829 | 8,662,392 |
| H1437 | 15,092,960 |
| H2009 | 22,466,426 |
| HCC1937 | 24,292,931 |
| HCC1954 | 15,966,684 |
| **Total** | **170,496,417** |

Sidecar每 alignment 一列 `CHROM/START/END/QNAME/FLAG/MAPQ/CIGAR digest/HP/PS`。以 100–160 bytes/row 作保守 planning range，未壓縮約 **15.9–25.4 GiB**；script 同時保留 raw TSV 與 bgzip TSV，加 index/VCF/log 後預估 **20–40 GiB**，設定 **60 GiB hard planning cap**。

歷史 canonical tagged BAM合計為 **1.674 TiB**。如果 FIFO 失效而改落 regular BAM，932 GiB free 不足；所以容量判定只對 FIFO-sidecar模式成立。現有 run root 的 `*_production.bam` 實際檔型為 named pipe (`p`)，不是 regular BAM。

歷史 LongPhase logs（含 truth-BED舊 run，只能作 lower-baseline）各 sample時間：

| Dataset | Historical seconds | Hours |
|---|---:|---:|
| HCC1395 | 5,893 | 1.64 |
| HCC1395_DORADO | 8,039 | 2.23 |
| COLO829 | 6,112 | 1.70 |
| H1437 | 6,113 | 1.70 |
| H2009 | 16,820 | 4.67 |
| HCC1937 | 13,451 | 3.74 |
| HCC1954 | 11,040 | 3.07 |

總和 18.75 sample-hours；historical dynamic parallel=2 理想 makespan約 10.7h。Clean run移除 truth-BED、需產/驗完整 sidecar，且目前 NFS已飽和，因此實務 wall projection為 **14–24h**；第一個 completed sample 的 `/usr/bin/time -v` 必須取代此估計。

### 5.2 Sidecar-aware layered full run

Resource baseline採既有、已完成但使用歷史 restricted tags 的 run：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/
```

- wall：2026-07-10 23:25:01 → 2026-07-11 02:57:42，約 **3h32m41s**。
- parallel：2 samples；每 sample 5 chromosome splits（最多約10個 BAM readers）。
- output：**484 MiB**。
- per-split peak RSS：多數 0.14–0.36 GiB，最大約 **0.97 GiB**。
- 7/7 verifier PASS只證此 run 的 engineering artifacts；不能替代 clean sidecar結果。

Sidecar-aware path增加 tabix sidecar lookup與 exact identity accounting，尚無 full-run實測；預估 **4–8h wall、<2 GiB output**。Clean tagging與layered不得重疊，因此 end-to-end planning range約 **18–32h wall**。

## 6. Consumer gap

### 6.1 已存在的 consumer capability

`sm_linkage_genomewide.py` 已讀 `SM_READ_TAG_SIDECAR`，以 `QNAME+reference coordinates+FLAG+CIGAR digest` 作 alignment identity，並輸出：

```text
sidecar_configured
sidecar_exact_matches
sidecar_missing
sidecar_conflicts
```

`run_layered_7samples_newbb.sh` 也已將 sidecar path傳入 MLHP，並要求每 part：

```text
sidecar_missing == 0
sidecar_conflicts == 0
sidecar_exact_matches == alignment_group_exposures
read_tag_source != BAM_HP_PS
```

### 6.2 目前仍不可 launch 的 handoff gap

截至截點，active：

```text
InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json
```

仍是 schema 2.0，7/7 的 `read_tag_sidecar`、`caller_raw_vcf`、`longphase_recalibrated_all_vcf` 等 production artifact fields 尚未指向完成 run。Runner現在會 fail-closed於缺檔，而不是 silently fallback；這是正確安全行為，但表示 **consumer code ready ≠ launch-ready manifest**。

此外，active production run啟動後工作樹仍可能變動；截點執行：

```bash
sha256sum -c /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1/code.sha256
```

四項皆 `OK`，但必須在 run complete後再讀回一次；否則後段 validator可能與 run-start hash不一致。

## 7. 三個 strongest failure modes、falsifier 與 safe threshold

| Rank | Failure mode | 為何最強 | Falsifier | Safe probe / launch threshold |
|---:|---|---|---|---|
| 1 | **Duplicate/overlap造成I/O collapse或兩個run爭同一input/output** | 已直接觀察同一 production run active；兩job已約116 MB/s NFS，host iowait 57–60%，load>logical CPUs。 | 新截點無同類process，run root唯一，5-minute baseline iowait<20%。 | 現有 run期間不再 launch；未來 `SM_PARALLEL_SAMPLES≤2`；禁止與layered重疊；launch前 big8 NFS baseline <80 MB/s。 |
| 2 | **FIFO模式退化／中斷，誤落完整 BAM或停在 START** | Persisted tagged BAM預估1.674 TiB > big7 free 932 GiB；named pipe/capture任一端死亡可留下不完整run。 | 每個 output prefix在執行期皆 `test -p`；結束後只有0-byte consumed FIFO標記、完整sidecar；run 7/7 PASS且無 START-only。 | Sidecar總投影≤60 GiB；big7 free≥500 GiB；任一 regular `*_production.bam` 或 output growth>60 GiB即停在NO-GO，由操作者另行處置。 |
| 3 | **Production sidecar未正確接入layered，或 identity不全卻產tree** | 目前manifest尚無7/7 sidecar paths；這會阻擋run。若未來放鬆gate，可能回退歷史BAM tags而重現truth-BED confound。 | schema 2.1 manifest 7/7 production validation PASS；bounded consumer輸出 exact=exposures、missing=0、conflicts=0、source不是BAM。 | 先用 **COLO829（最小8.66M mapped）** 作1-sample probe：wall≤60min、每split RSS≤2GiB、output≤1GiB、5 parts全部exact；任一不符不得擴7 datasets。 |

### Safe probe 的額外 acceptance

1. `verification_summary.json`: production 7/7 `pass=true`、truth flags absent、VCF keys missing/extra=0、unknown HP=0、duplicate exact conflicts=0。
2. `sha256sum -c code.sha256` 與每 sample `input.sha256/output.sha256` 全部 OK；canonical input stat digest無漂移。
3. 新 layered manifest必須只讀 immutable production run root，不讀 mutable working output。
4. COLO829 probe前 host available RAM≥128 GiB、big7 free≥500 GiB、iowait 5-minute average<20%、無另一個 layered/LongPhase job。
5. COLO829通過後，7-dataset full run最多 parallel=2；若第一對 sample wall或RSS超 probe projection 2×，降回PROBE，不宣稱full validation。

## 8. 觀察命令 ledger

所有命令均唯讀；沒有執行 `kill`、`renice`、signal、寫入run root或啟動 pipeline。

| 目的 | 輸入 | 命令摘要 | 實際結果／輸出 |
|---|---|---|---|
| Host/time | bip7 | `date; hostname; uptime` | cutoff 04:11:58 +0800；load 70.58/67.89/50.20 |
| CPU/RAM | host | `lscpu`; `free -h` | 48 logical CPU；503 GiB total、471 GiB available |
| Disk/inode | big7/big8 | `df -hT`; `df -ih` | 932/845 GiB free；inode充足 |
| Process exclusivity | host process table | `ps -eo ... | rg <exact pipeline patterns>` | wrapper 656363；LongPhase 656465/656470；無layered process |
| Per-process resources | PIDs 656465/656470 | `pidstat -dru -p ... 1 3` | each ~57 MB/s read、6.6–7.0 GiB RSS、~1 core |
| I/O pressure | host/NFS | `mpstat`; `vmstat`; `nfsiostat` | iowait 57–60%；big8 ~116 MB/s |
| Input writers | manifest exact paths | `jq ... | sort -u | xargs lsof --` | only `r` FDs；no writer observed |
| Input identity | same 29 paths | `stat -L ... | sha256sum` twice | equal；digest `681514...33b8` |
| Input volume | 14 BAMs | `stat -L` sum | 2.485 TiB raw BAM input |
| Alignment rows | 7 tumor BAI | `samtools idxstats` | 170,496,417 mapped alignments |
| Historical clean-tag baseline | 7 exact canonical logs | parse recorded `... Ns` fields | 18.75 sample-hours；per-sample table見§5.1 |
| Layered baseline | exact prior run root | `run_status.tsv`, `du -sh`, `/usr/bin/time` logs | 3h32m41s；484 MiB；max split RSS~0.97 GiB |
| Active source identity | exact production run | `sha256sum -c code.sha256` | 4/4 OK at cutoff |

## 9. 最終紅隊結論

**現在不是「資源不足所以永遠不能跑」，而是「相同job已經在跑、I/O已飽和、handoff artifacts尚未完成，所以不能再launch」。** FIFO sidecar把1.674 TiB persisted-BAM風險降到約20–40 GiB，RAM也有充足餘裕；資源設計方向可行。但截至截點，唯一安全決策是：

```text
NO-GO additional launch
  → existing production run 7/7 immutable validation
  → schema-2.1 manifest + COLO829 bounded consumer PROBE
  → only then conditional GO for 7-dataset sidecar-aware layered full run
```

此 verdict 的最高優先 falsifier不是再跑更多計算，而是等現有 run 產生第一份與最終 7/7 可讀 evidence，並以 §7 threshold readback。
