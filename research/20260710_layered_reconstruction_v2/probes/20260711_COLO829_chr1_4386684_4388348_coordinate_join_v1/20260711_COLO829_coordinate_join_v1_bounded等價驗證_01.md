<!--
建立時間: 2026-07-11 04:40 +08:00
目標: 以真實 COLO829 alignment 驗證 direct BAM HP/PS 與 raw BAM + coordinate_join_v1 sidecar 的 bounded consumer 等價性
處理範圍: ⚠ PARTIAL SCOPE — COLO829 × chr1:4,386,684-4,388,348 × 2 sSNV × 35 primary alignment exposures
完整 scope 未驗證: production no-truth sidecar 7 datasets × chr1-22；本 probe 未涵蓋真實 supplementary alignment
補強計劃: synthetic supplementary/conflict/missing tests + production 7/7 global uniqueness/validation + 35/35 split runtime conservation
關聯檔案:
  - InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/probe_read_tag_sidecar_equivalence.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_read_tag_sidecar_contract.py
task_type: A bounded probe supporting a Type-B comprehensive gate
goals: G4,G5
cycle_id: 20260710-2244-layered-reconstruction-v2
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
證據等級: L2 ⭐⭐⭐⭐（真實 bounded data + synthetic edge cases；尚非 production 7/7）
-->

# COLO829 coordinate_join_v1 bounded 等價驗證

> ⚠ **PARTIAL SCOPE**：本結果只證明指定 COLO829 區段的 join 方法等價；historical tagged BAM 的 truth 設定不可作 clean production tag 證據，也不可直接外推至 7 datasets × chr1-22。

## 結論

在 `chr1:4,386,684-4,388,348` 的兩個 ClairS PASS sSNV，direct historical tagged BAM 與 raw BAM + 同源 HP/PS sidecar 產生完全相同的 calls、HP、PS 與 family payload。35 個被分析 alignment 全部 exact matched，missing=0、conflict=0。

| Claim | 結果 | Evidence |
|---|---:|---|
| alignment exposure 守恆 | 35 direct / 35 joined | ⭐⭐⭐⭐ L2 bounded real data |
| exact join | 35/35 | ⭐⭐⭐⭐ L2 bounded real data |
| missing / conflict | 0 / 0 | ⭐⭐⭐⭐ L2 bounded real data |
| combined payload digest | `b466e584ab820448769bd1aede31434021f3c6ce1762788afc7363c0fe2d71d7`（兩路相同） | ⭐⭐⭐⭐ L2 bounded real data |
| sidecar identity uniqueness | 35 rows、0 duplicate | ⭐⭐⭐⭐ L2 bounded real data |

## 輸入

| Role | Absolute path |
|---|---|
| historical tagged BAM | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/COLO829_tagged.bam` |
| raw alignment payload BAM | `/big8_disk/data/COLO829/ONT_PAO/PAO29420.bam` |
| ClairS PASS VCF | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz` |
| sSNV positions | `chr1:4386684,4388348` |

## 執行命令

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/probe_read_tag_sidecar_equivalence.py \
  --tagged-bam /big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/COLO829_tagged.bam \
  --raw-bam /big8_disk/data/COLO829/ONT_PAO/PAO29420.bam \
  --somatic-vcf /big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz \
  --chrom chr1 --positions 4386684,4388348 \
  --mapq-min 20 --baseq-min 0 \
  --output-dir InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1
```

Exit code：`0`。

實際輸出：

```text
REAL-DATA SIDECAR PROBE: PASS -> .../equivalence_probe.json
{"calls_sha256_equal": true, "hp_sha256_equal": true, "payload_equal": true,
 "payload_sha256_equal": true, "ps_sha256_equal": true,
 "sidecar_conflicts_zero": true, "sidecar_exact_matches_all_exposures": true,
 "sidecar_missing_zero": true}
```

## 輸出與 readback

| Artifact | Path / digest |
|---|---|
| machine-readable result | `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json` |
| compressed sidecar | `COLO829.coordinate_join_v1.read_tags.tsv.gz` / SHA-256 `489a8a04390e68484b1f05d3cda27142e6aa11b8fc8f87ee87b29728ef4128ec` |
| tabix index | `COLO829.coordinate_join_v1.read_tags.tsv.gz.tbi` / SHA-256 `80ee5038fc787f1778c11568feac16d3a9178d156d7250d19183ac4f26be38a6` |
| direct payload | SHA-256 `b466e584ab820448769bd1aede31434021f3c6ce1762788afc7363c0fe2d71d7` |
| joined payload | SHA-256 `b466e584ab820448769bd1aede31434021f3c6ce1762788afc7363c0fe2d71d7` |

HP counts：`1-1=11`、`2-1=11`、`3=11`、`.`=2；兩條路徑完全相同。PS-bearing alignments：`1-1=11`、`2-1=8`。

## 判斷細節

1. Sidecar key 為 `QNAME + RNAME + START0 + END0 + FLAG + BLAKE2b8(CIGAR)`；本區段 35 rows 無 duplicate，因此 actual-data key 是一對一。
2. HP/PS 只由 sidecar提供；base、quality、MAPQ與 pileup 仍由 frozen raw BAM提供。
3. `coordinate_join_v1` 不含 SEQ/QUAL digest，因此只能在 BAM identity另行凍結、global duplicate=0、runtime exact conservation的條件下使用；不可稱為 payload identity v2。
4. Synthetic contract另覆蓋同 qname primary/supplementary、缺列與衝突列；3/3 tests PASS。

## Scope Limitation

- 已驗證：COLO829、chr1兩個 sSNV、35 個 primary alignment，direct vs sidecar payload逐鍵等價。
- 未驗證：production no-truth tag本身、其他區段、其他 6 datasets、真實 supplementary/secondary、全基因組 global payload identity。
- 推論限制：本結果只解鎖 join 方法的 bounded gate；不能單獨解鎖新的 7-dataset full run。
- 補強動作：等待 production 7/7 `duplicate_exact_alignment_rows=0`、truth flags absent、hash binding；再由 v3 frozen validator與 35/35 split runtime gate重驗。
