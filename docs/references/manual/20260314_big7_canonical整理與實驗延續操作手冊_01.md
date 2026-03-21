<!--
建立時間: 2026-03-14
目標: 提供 big7 canonical materialization、paired/TO 執行入口與檢查方式
-->

# Big7 Canonical 整理與實驗延續操作手冊

## 1. 先整理既有 archive

```bash
cd /bip7_disk/liaoyoyo2001/InterSubMod
./.venv/bin/python scripts/analysis/materialize_big7_canonical.py --max-workers "$(nproc)"
```

主要輸出：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/master_run_manifest.tsv`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/master_experiment_matrix.tsv`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/master_report.md`

## 2. 跑新的 paired benchmark

```bash
cd /bip7_disk/liaoyoyo2001/InterSubMod
./scripts/pipeline/run_benchmark.sh --sample HCC1395 --mode s-pure --threads "$(nproc)" --min-free-gb 800
./scripts/pipeline/run_benchmark.sh --sample HCC1395 --mode s-pure --vcf-source pileup --threads "$(nproc)" --min-free-gb 800
```

新 run 預設會寫到：

- `canonical/HCC1395/paired_full/...`
- `canonical/HCC1395/paired_pileup/...`

## 3. 跑新的 pure research round

```bash
cd /bip7_disk/liaoyoyo2001/InterSubMod
./scripts/analysis/run_pure_research_round.sh \
  --sample-set all_pure \
  --mode s-pure \
  --run-pipeline \
  --threads "$(nproc)"
```

研究 round 預設會寫到：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/`

## 4. 跑 tumor-only pilot

```bash
cd /bip7_disk/liaoyoyo2001/InterSubMod
./scripts/analysis/run_longphase_to_intersubmod_pilot.sh --threads "$(nproc)" --min-free-gb 800
```

TO pilot 預設會寫到：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/`

## 5. 檢查是否齊全

先看總表：

```bash
column -ts $'\t' /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/master_experiment_matrix.tsv | less -S
```

再看單一 run：

```bash
column -ts $'\t' /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260307_HCC1395_paired_full_full/metrics/completeness.tsv
```

## 6. CPU / 空間原則

1. 可用滿全機 CPU；預設 `THREADS=$(nproc)`。
2. `materialize_big7_canonical.py` 也可直接用 `--max-workers "$(nproc)"`；預設已改為全機 CPU。
3. `run_benchmark.sh` 與 `run_longphase_to_intersubmod_pilot.sh` 會檢查輸出卷最低剩餘空間；預設門檻 `800G`，可用 `--min-free-gb` 覆寫。
4. 大型 BAM/VCF 不重複複製，只在 canonical bundle 寫 manifest 與 summary。
5. 若 `completeness_state=partial` 且 `blocking_reason=missing_tagged_bam`，代表舊 run 的 tagged BAM 沒留在 archive，需要正式重跑補齊，不是統計表遺失。
6. 若 `expected=false`，代表目前 source tree 沒有該 caller/model 路線，不能把它當成 runnable gap。
7. 目前 `TO full-model` 在 source tree 中仍是 unavailable；`TO` 主線只應先排 `to_pileup`。
