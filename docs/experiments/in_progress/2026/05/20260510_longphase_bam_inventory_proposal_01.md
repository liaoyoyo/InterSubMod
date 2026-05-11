---
title: longphase-to BAM 輸出盤點 — Cleanup Proposal
date: 2026-05-10
status: pending_user_review
type: longphase_bam_inventory
classification: housekeeping_external_dirs
generated_by: AI scan (no deletion executed)
related:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260510_external_dirs_inventory_proposal_01.md
  - InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md
---

# longphase-to BAM 輸出盤點

> **Bottom line**：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/` 含 33 個 work 目錄（baseline → v2 → v2b → v3_fixed → v4_alt_guard → v5_somatic_fallback 多版本歷史），每個 active work 內含 ~270 GB tumor_tagged.bam。
> /tmp 已清乾淨（800 GB 災情已解決）。
> **可清候選總共 ~320-350 GB**（11 個 0-active-script-引用 work dir）：
> 主要：`pononly_v4_alt_guard/` (1 BAM = **267.6 GB**) + 10 個 ism_*/ 衍生目錄（~50-80 GB；含 50 萬-58 萬個 region-level small files）。
> **未執行任何刪除**；下列 bash command 用戶 review 後手動執行。

---

## §1 確認事項

### 1.1 /tmp 清空 ✅

`ls /tmp/` grep `longphase|.bam|.cram` = 0 hit。
2026-05-08 longphase-to 把 BAM 寫到 /tmp 寫滿 root volume 的 800 GB 災情已恢復。
v1.8 T1-2 `tmpdir_guard.sh` 已防再犯（pipeline launcher 強制 export TMPDIR=/big7_disk/...）。

### 1.2 longphase-to-mod/output/ 結構

33 個 work 目錄，多版本演進：

| 階段 | Work dir | mtime |
|---|---|---|
| 1. 基準 | `baseline/` | 2026-04-03 |
| 2. v2 | `pononly_v2/` (在 find 看到目錄但未 ls) | 推 4 月初 |
| 3. v2b | `pononly_v2b/` | 推 4 月初 |
| 4. v3_fixed | `pononly_v3_fixed/`、`v3f_no_pononly/`、`clairsto_v3fixed_work/` | 4-10 |
| 5. v4_alt_guard | `pononly_v4_alt_guard/` ⚠️ 0 active ref | 4-11 |
| 6. v5_somatic_fb | `pononly_v5_somatic_fallback/` | 4-12 |
| ISM 衍生 | `ism_*` × 多個 | 4-3 ~ 4-12 |
| 純度模擬 | `purity_06_simulation/` | 4 月 |

**每個 active work 主目錄結構**（用 baseline 為例）：
```
baseline/
├── tumor_tagged.bam       278.3 GB    ← 主 BAM（90%+ size）
├── tumor_tagged.bam.bai   118 MB
├── tumor_phased.vcf       655 MB
├── tumor_phased_GE.bed    840 KB
├── tumor_phased_LOH.bed   74 KB
├── haplotag.log           1.4 KB
└── run.log                2 KB
```

**ism_* 衍生結構**（用 ism_v3_fixed_tp 為例）：
```
ism_v3_fixed_tp/
├── debug/
├── filtered_snv_tp/{chr1..chr22,chrX,chrY}/  ← 24 染色體 region-level outputs
├── run.log                45 KB
├── significance_statistics.txt
└── significance_summary.csv  15 MB
```
(無 BAM，但檔案數量極大 — 30K-580K 個 region-level small files)

---

## §2 引用驗證 — Active vs 0-Reference

### 2.1 Active 主 work dirs（不可動）

| Work dir | scripts/ refs | 用途 |
|---|---|---|
| `baseline/` | **4** | v5_imbalance / v5_whole_genome / v3f_ablation 都引用 |
| `pononly_v2/` | 9 | v3f_ablation 對照組 |
| `pononly_v2b/` | 9 | v3f_ablation A5 |
| `pononly_v3_fixed/` | 1+ | run_vcf_all_snv 的 TO_TUMOR_BAM |
| `clairsto_v3fixed_work/` | 2+ | TO_VCF_TP / TO_VCF_FP |
| `v3f_no_pononly/` | 1 | v3f_ablation A3 |
| `pononly_v5_somatic_fallback/` | 2 | v5_imbalance / v5_whole_genome V5_BAM |
| `purity_06_simulation/{baseline_06,v3f_no_pononly_06,v5_06}/` | 1+ | v3f_ablation 純度模擬 B1/B3/B5 |
| `ism_baseline_*` (4 dirs) | 4 | active ISM 衍生 |
| `ism_pononly{,_v2b,_v2b_clairsto}_*` | 4 | active ISM 衍生 |

### 2.2 0-Reference candidates（可清，14 個 — du 完整 size 已拿到）

| Work dir | scripts/ refs | docs/ refs | Size | mtime |
|---|---|---|---|---|
| **`pononly/` (v1 無版本標號)** | **0** | 1 | **260 GB** ⭐⭐ | 2026-04-03 |
| **`pononly_v4_alt_guard/`** | **0** | 1 | **273 GB** ⭐⭐⭐ | 2026-04-11 |
| `v3f_ablation_ism_06/` | **0** | 1 | 1.8 MB | 4 月 |
| `v5_93_purity_fix/` | **0** | 1 | 627 MB | 4 月 |
| `ism_v3_fixed_tp/` | **0** | 0 | 8.1 GB | 4-10 |
| `ism_v3_fixed_fp/` | **0** | 0 | 1.2 GB | 4-10 |
| `ism_v3fixed_clairsto_tp/` | **0** | 0 | 7.5 GB | 4-10 |
| `ism_v3fixed_clairsto_fp/` | **0** | 0 | 3.1 GB | 4-10 |
| `ism_v4altguard_clairsto_tp/` | **0** | 0 | 7.5 GB | 4-11 |
| `ism_v4altguard_clairsto_fp/` | **0** | 0 | 3.1 GB | 4-11 |
| `ism_v5_somatic_fb_clairsto_tp/` | **0** | 0 | 7.5 GB | 4-12 |
| `ism_v5_somatic_fb_clairsto_fp/` | **0** | 0 | 3.1 GB | 4-12 |
| `ism_longphase_to_pass_tp/` | **0** | 0 | 7.5 GB | 4-10 |
| `ism_longphase_to_pass_fp/` | **0** | 0 | 5.2 GB | 4-10 |

**子小計**：260 + 273 + 0.6 + (~54 GB ism_*) = **~588 GB 可清** ⭐⭐⭐

### 2.3 為何這些 0-reference 仍存在？

- **`pononly_v4_alt_guard`**：v4_alt_guard 是 ClairS pipeline 中間版本（v3_fixed → v4 → v5），已被 `pononly_v5_somatic_fallback/` 取代
- **`ism_v3_fixed_*` / `ism_v3fixed_clairsto_*`**：ISM 跑 v3_fixed 階段的 region-level outputs，已被新 cycle 取代但 dir 留著
- **`ism_v4altguard_clairsto_*`**：v4_alt_guard 的 ISM 衍生 — 同上 obsolete
- **`ism_v5_somatic_fb_clairsto_*`**：v5_somatic_fallback 的 ISM 衍生 — 主目錄 `pononly_v5_somatic_fallback/` 仍 active，但 ism_v5 衍生已被新 ISM run 取代
- **`ism_longphase_to_pass_*`**：LongPhase pass test 的 ISM 衍生 — 已過時

---

## §3 推理與證據鏈

### 3.1 為何 active 的 work 不能搬

| Active work | scripts 來源 | 影響 |
|---|---|---|
| `baseline/` | `scripts/analysis/v5_imbalance_improvement.py` line 39 = BL_BAM | 搬 = v5 imbalance test 跑不動 |
| `clairsto_v3fixed_work/` | `scripts/run_vcf_all_snv.sh` line 41-42 | 搬 = run_vcf_all_snv 跑不動 |
| `pononly_v3_fixed/` | `scripts/run_vcf_all_snv.sh` line 43 | 搬 = TO_TUMOR_BAM 失效 |
| `pononly_v5_somatic_fallback/` | `scripts/analysis/v5_imbalance_improvement.py` line 40 | 搬 = V5_BAM 失效 |
| `pononly_v2`, `v2b`, `v3f_no_pononly` | `scripts/analysis/v3f_ablation_phased_vcf_f1.sh` 對照組 | 搬 = ablation 對照失效 |

### 3.2 為何 0-reference 安全可清

- `grep -rli` 對 InterSubMod/scripts/ 全掃 → 0 命中
- mtime 1 個月舊（4 月初）— 1 個月無修改
- 路徑名稱本身（v4_alt_guard / v3_fixed / longphase_to_pass）顯示是中間版本或測試 pass
- 主版本 `pononly_v5_somatic_fallback/` (active) 已取代 v3_fixed / v4_alt_guard 的角色

---

## §4 建議行動（用戶手動，不要直接刪）

### 🟢 Tier 1 — 主大頭（267.6 GB 立即釋放）

```bash
# pononly_v4_alt_guard 含 1 個 268 GB tumor_tagged.bam，0 active ref
mkdir -p /big7_disk/liaoyoyo2001/_archive_2026-05/longphase_obsolete
mv /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v4_alt_guard \
   /big7_disk/liaoyoyo2001/_archive_2026-05/longphase_obsolete/
```

### 🟡 Tier 2 — ism_* 衍生（~50-80 GB，但檔案 ~150 萬個 inode 釋放）

```bash
# 注意：每個 ism_*/ 含 30K-580K small files；mv 也會慢，但釋放大量 inode
for d in ism_v3_fixed_tp ism_v3_fixed_fp \
         ism_v3fixed_clairsto_tp ism_v3fixed_clairsto_fp \
         ism_v4altguard_clairsto_tp ism_v4altguard_clairsto_fp \
         ism_v5_somatic_fb_clairsto_tp ism_v5_somatic_fb_clairsto_fp \
         ism_longphase_to_pass_tp ism_longphase_to_pass_fp; do
  mv "/big7_disk/liaoyoyo2001/longphase-to-mod/output/$d" \
     /big7_disk/liaoyoyo2001/_archive_2026-05/longphase_obsolete/
done
```

### Final delete（你 review 後）

```bash
# 確認封存區內容無問題後，硬刪釋放空間
ls -lah /big7_disk/liaoyoyo2001/_archive_2026-05/longphase_obsolete/
# 若無問題（你看一眼確認後）：
rm -rf /big7_disk/liaoyoyo2001/_archive_2026-05/longphase_obsolete
```

---

## §5 還沒拿到 size 的目錄（不在 cleanup 列表）

部分 active work du timeout，但**它們有 active script 引用，本來就不動**：
- `baseline/` (~278 GB)
- `pononly_v2{,b}/`、`pononly_v3_fixed/`、`pononly_v5_somatic_fallback/`、`v3f_no_pononly/`
- `clairsto_v3fixed_work/`（已知小，純 VCF）
- `purity_06_simulation/{...}/`
- `ism_baseline*`、`ism_pononly*`

預估每 active work ~250-280 GB；7 個 active = **~1.7-2.0 TB**（不可動）

---

## §6 累計 cleanup 提案彙總（與其他 inventory 報告整合）

| 報告 | 提案 |
|---|---|
| InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md | InterSubMod 內 .agents/ + .codex + small synthesis 目錄 |
| InterSubMod/docs/experiments/in_progress/2026/05/20260510_external_dirs_inventory_proposal_01.md | data/bam 空 + InterSubMod_big7_runbook 7.5GB + bip8/big8 archive |
| **本報告**（longphase BAM）| **pononly_v4_alt_guard 267.6 GB + 10 ism_*/ ~50-80 GB** |

| 已驗證可清總計（不含尚未確認的 bip8/big8 archive）| Size |
|---|---|
| `cmp_ClairS_DeepSomatic/` (前述 inventory) | 886 GB |
| `pononly_v4_alt_guard/` （本報告主大頭）| 267.6 GB |
| ism_*/ × 10 (本報告) | ~50-80 GB |
| `InterSubMod_big7_runbook/` 7.5 GB | 7.5 GB |
| LongPhase test cram + GRCh38 重複 | 23 GB |
| 其他小目錄 | <1 GB |
| **總計** | **~1.24 TB+** |

---

## §7 References

- v1.8 retro：`InterSubMod/docs/experiments/in_progress/2026/05/20260510_v1.8_implementation_retro_01.md`
- 800 GB /tmp 災情教訓：memory `feedback_tmp_disk_full_pipeline_pitfall.md`
- 引用驗證：grep 全掃 `InterSubMod/scripts/ .claude/ docs/ research/`
- BAM size 樣本：`baseline/tumor_tagged.bam` = 278 GB（每 active work 同等規模）
