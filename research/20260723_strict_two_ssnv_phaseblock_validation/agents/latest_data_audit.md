<!--
建立時間: 2026-07-23 13:18 +08:00
最後重跑: 2026-07-23 13:24 +08:00
目標: 獨立稽核截至當下最新 strict exact-PS/HP production 輸出、freshness 與 graph invariants
處理範圍: canonical 7 technical datasets readiness；對已完成 summary 的 5 datasets 全量重算 chr1-22
關聯檔案: InterSubMod/research/20260723_strict_two_ssnv_phaseblock_validation/agents/audit_latest_strict_outputs.py
-->

# 最新 strict production 資料與 freshness 獨立稽核

> **PARTIAL：截至 2026-07-23 13:18:04（Asia/Taipei），已有 5/7 technical datasets 具 chr1–22 strict receipts 與 all-PASS summary；這 5 組共 110 個 chromosome receipts 的獨立 graph／hash／mass 稽核為 0 violations。不可寫成 7/7 已完成。**  
> 任務類型：B Comprehensive validation；服務 G1／G4／G5。

## 1. Freshness 與完成狀態

| technical dataset | extraction | strict chromosome receipts | genome summary | 本輪納入 |
|---|---:|---:|---|---|
| HCC1395 | 22/22 | 22/22 | all-PASS | 是 |
| COLO829 | 22/22 | 22/22 | all-PASS | 是 |
| H1437 | 22/22 | 22/22 | all-PASS | 是 |
| HCC1395_DORADO | 22/22 | 22/22 | all-PASS | 是 |
| H2009 | 22/22 | 22/22 | all-PASS | 是 |
| HCC1937 | **5/22** | 0/22 | 尚無 | 否，extraction 正在進行 |
| HCC1954 | 0/22 | 0/22 | 尚無 | 否 |

狀態快照時間為 `2026-07-23T13:18:04+08:00`。HCC1937 的 extraction 輸出仍在進行／重啟狀態，因此沒有把未完成檔案混入 invariant 分母。

`audit_latest_strict_outputs.py` 於 `2026-07-23T13:24:18+08:00` 再次完整執行；五組 summary、110 個 chromosome receipts 與 772 個 unique identity files 的結果不變，`violations_total=0`。13:24:30 觀察 HCC1937 為 5/22 completed extraction receipts、0/22 strict receipts；HCC1954 為 0/22，因此仍不納入完整稽核分母。

完成資料的最新 strict artifact 時間：

| dataset | latest artifact（UTC） | 台北時間 |
|---|---|---|
| HCC1395 | 2026-07-23 02:54:29 | 10:54:29 |
| HCC1395_DORADO | 2026-07-23 03:44:41 | 11:44:41 |
| COLO829 | 2026-07-23 03:47:00 | 11:47:00 |
| H1437 | 2026-07-23 04:02:22 | 12:02:22 |
| H2009 | 2026-07-23 05:13:15 | 13:13:15 |

## 2. 稽核資料粒度與定義

- candidate catalog grain：`dataset × chromosome × site_index`。
- primary edge grain：同一 `dataset × chromosome × exact nonmissing PS × HP1/HP2` 內的一對不同 sSNV endpoints，且 `support_total>=3`。
- `W` grain：primary threshold=3 endpoint graph 的 connected component；只有 `k>=2` 是 tree-eligible，`k=1` 是 singleton abstention。
- 「無 missing PS/HP」只表示 missing PS、nonprimary HP 沒有進 primary edge／W；不表示 raw molecule 資料本身沒有這些列。

## 3. 五組完成資料的重算結果

| dataset | candidate S | primary edges | all components | singleton abstain | tree-eligible W | k>12 W | max k |
|---|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 79,687 | 76,202 | 39,846 | 28,384 | 11,462 | 90 | 153 |
| COLO829 | 37,788 | 41,736 | 42,189 | 28,256 | 13,933 | 30 | 40 |
| H1437 | 77,080 | 140,566 | 31,588 | 15,262 | 16,326 | 904 | 138 |
| HCC1395_DORADO | 79,739 | 19,453 | 43,310 | 36,482 | 6,828 | 34 | 106 |
| H2009 | 154,465 | 855,116 | 49,392 | 25,199 | 24,193 | 4,206 | 567 |
| **合計（dataset-scoped）** | **428,759** | **1,133,073** | **206,325** | **133,583** | **72,742** | **5,264** | — |

守恆關係：

- `206,325 all components = 133,583 singleton + 72,742 tree-eligible W`。
- primary threshold=3 active node memberships 共 542,454；逐 chromosome 驗證每個 membership 恰屬一個 component，且 component 的 `k` 等於實際 membership 數。
- 全部 observed endpoint-pair rows 共 1,564,789；state mass 共 20,410,212，逐列驗證 `support_total = RR + RA + AR + AA`。

## 4. Invariant 結果

| invariant | 稽核數量 | violations |
|---|---:|---:|
| primary edge `support_total>=3` 且 threshold flag 一致 | 1,133,073 edges | 0 |
| edge 兩端 site index／position 不同 | 1,133,073 edges | 0 |
| edge 兩端均精確映射 candidate sSNV catalog | 1,133,073 edges | 0 |
| edge `linkage_basis == PS_HP{hp}`、PS nonmissing、HP∈{1,2} | 1,133,073 edges | 0 |
| edge 兩端在同 exact PS／HP membership 與同 component | 1,133,073 edges | 0 |
| tree-eligible `W` 的 `k>=2` | 72,742 W | 0 |
| `W` members 全在同 chromosome／exact PS／HP | 72,742 W | 0 |
| `W` 由 retained primary edges 獨立重建後 connected | 72,742 W | 0 |
| singleton 均為 `k=1` 且不標 tree-eligible | 133,583 components | 0 |
| chromosome receipt `all_pass` 與所有 checks=true | 110 receipts | 0 |
| receipt `.sha256` sidecar | 110 receipts | 0 |
| receipt input/output identity（size＋SHA-256） | 880 identities | 0 |
| genome summary → chromosome receipt identity | 110 identities | 0 |
| genome summary aggregate 與 raw TSV 重算一致 | 5 summaries | 0 |

實際共 hash 772 個 unique files；builder／graph core 等共享檔案只 hash 一次，但每一份 receipt identity 都有逐項比對。

## 5. `phase_set` 特殊缺失 token 稽核

額外直接掃描 5 組完成資料的 authoritative raw molecule rows，以及所有 strict membership／edge rows。特殊 token 集合為：空字串、`.`、`NA`、`N/A`、`NaN`、`None`、`null`、`UNKNOWN`、`UNK`、`MISSING`、`NOT_AVAILABLE`、`NOT APPLICABLE`（非空 token 不分大小寫）。

| dataset | raw rows | raw 空字串 | raw `.` | raw 其他 sentinel | strict membership 空／sentinel | strict edge 空／sentinel |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 3,201,802 | 1,069,698 | 0 | 0 | 0 / 0 | 0 / 0 |
| COLO829 | 988,029 | 138,859 | 0 | 0 | 0 / 0 | 0 / 0 |
| H1437 | 2,899,543 | 453,907 | 0 | 0 | 0 / 0 | 0 / 0 |
| HCC1395_DORADO | 3,507,573 | 1,504,313 | 0 | 0 | 0 / 0 | 0 / 0 |
| H2009 | 7,544,561 | 1,222,875 | 0 | 0 | 0 / 0 | 0 / 0 |
| **合計** | **18,141,508** | **4,389,652** | **0** | **0** | **0 / 0** | **0 / 0** |

Strict membership 共 2,169,816 rows（4 個 thresholds 的 registry），strict edge 共 1,564,789 rows；兩層所有 `phase_set` 都是非缺失、非 sentinel。

再以 `^[0-9]+$` 且數值 `>0` 檢查所有非空、非 sentinel token：raw molecule、strict membership、strict edge 的「非正整數 PS」皆為 **0**。

**判定：目前五組資料不受 builder 只以 `if not phase_set` 排除空字串的影響**，因為 raw data 的缺失值實際全是空字串，其餘 PS 皆為正整數；未發現會被誤當合法 PS 的非空 sentinel 或異常 token。不過這只是目前資料觀察，不是 schema 保證；production gate 仍建議明確拒絕上述 sentinel，並固定 positive-integer validation。

### 5.1 Candidate catalog 的 SNV allele contract

另對 5 組完成資料的 428,759 筆 candidate catalog 全量檢查 `ref`／`alt`：

| dataset | catalog rows | `len(ref)!=1` | `len(alt)!=1` | ref 非 A/C/G/T | alt 非 A/C/G/T | ref=alt |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 79,687 | 0 | 0 | 0 | 0 | 0 |
| COLO829 | 37,788 | 0 | 0 | 0 | 0 | 0 |
| H1437 | 77,080 | 0 | 0 | 0 | 0 | 0 |
| HCC1395_DORADO | 79,739 | 0 | 0 | 0 | 0 | 0 |
| H2009 | 154,465 | 0 | 0 | 0 | 0 | 0 |
| **合計** | **428,759** | **0** | **0** | **0** | **0** | **0** |

因此本輪 strict graph 的 candidate universe 全部是 REF／ALT 皆為單一 canonical DNA base 且互異的 biallelic SNV；未發現 indel、MNV、symbolic allele 或 `N`。

Catalog 檢核命令型態：

```bash
for f in <dataset>/chromosomes/chr*/extraction/*.site_catalog.tsv.gz; do
  gzip -dc "$f" |
    awk -F '\t' 'NR>1 {
      if (length($4)!=1 || $4!~/^[ACGT]$/) bad_ref++;
      if (length($5)!=1 || $5!~/^[ACGT]$/) bad_alt++;
      if ($4==$5) same++
    } END {print NR-1, bad_ref+0, bad_alt+0, same+0}'
done
```

重算命令使用：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  research/20260723_strict_two_ssnv_phaseblock_validation/agents/audit_phase_set_tokens.py \
  --sample-root HCC1395=.../hcc1395_strict_regions_v2 \
  --sample-root COLO829=.../COLO829/strict_regions_v1 \
  --sample-root H1437=.../H1437/strict_regions_v1 \
  --sample-root HCC1395_DORADO=.../HCC1395_DORADO/strict_regions_v1 \
  --sample-root H2009=.../H2009/strict_regions_v1 \
  --workers 8
```

## 6. 輸入、命令與實際輸出

主要輸入：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/samples/{COLO829,H1437,HCC1395_DORADO,H2009}/strict_regions_v1`
- 各 chromosome receipt 所綁定的 `site_catalog.tsv.gz`、`endpoint_edges.tsv.gz`、`components.tsv.gz`、`site_component_membership.tsv.gz`。

執行命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  research/20260723_strict_two_ssnv_phaseblock_validation/agents/audit_latest_strict_outputs.py \
  --sample-root HCC1395=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2 \
  --sample-root COLO829=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/samples/COLO829/strict_regions_v1 \
  --sample-root H1437=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/samples/H1437/strict_regions_v1 \
  --sample-root HCC1395_DORADO=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/samples/HCC1395_DORADO/strict_regions_v1 \
  --sample-root H2009=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/samples/H2009/strict_regions_v1
```

H2009 增量實際輸出為 `22` chromosome receipts、`156` unique identity files、`1/1` sample all-PASS、`violations_total=0`。與先前 4 組合併後的精確累計：

```json
{
  "chromosome_receipts_audited": 110,
  "identity_files_hashed_unique": 772,
  "sample_roots_all_pass": 5,
  "sample_roots_audited": 5,
  "violations_total": 0
}
```

## 7. 判定與下一道 gate

- **可用範圍：** 五個已完成 technical datasets 的 strict region engineering evidence 可引用。
- **不可用範圍：** 不可寫「七資料集皆通過」；HCC1937、HCC1954 尚未完成同一層級的 receipt chain。
- **完成 gate：** 每個剩餘 dataset 必須先有 22/22 chromosome receipts、genome summary `all_pass=true`，再用本腳本重跑；應達 `violations_total=0` 才能升級到 7/7。
- **freshness 風險：** production 正在寫入；後續報告引用時應重新快照，不應把本文件的 5/7 當永久狀態。

## Claim ceiling

本稽核證明五組完成資料的 primary strict graph、receipt identity 與集合守恆正確；不證明 unique mutation-state topology、真實 clone/subclone 數量、A–C 直接共現，亦不驗證 k>12 block partition 的 downstream tree eligibility。
