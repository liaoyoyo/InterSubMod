<!--
建立時間: 2026-06-12
報告類型: 大規模搜尋與驗證 handoff（給另一 session 執行；論文 claim 確認用）
任務類型: D handoff — 把「論文主軸 proof-of-concept → validated」所需的大規模驗證項目規格化
狀態: handoff spec / 待用戶確認完整性
data_sources:
  - docs/paper_focus/00_共識證據台帳_20260612_01.md (4 校正 + 13 矛盾 + open gap)
  - docs/paper_focus/00_論文章節材料對應總表_20260612_01.md (§D gap)
  - docs/thesis/20260612_thesis_alignment_table_01.md (修訂主軸 + R 脊柱)
  - docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md
provenance_note: 本檔列「待大規模驗證」項目，不產新數字；每項標現狀 tier + 方法腳本 + 成功判準。執行 session 跑完須回填真值（§13：先寫檔→讀回→才寫報告；勿填預期值）。
-->

# 大規模搜尋與驗證 Handoff — Subclonal reconstruction 論文 claim 確認

> **給誰**：執行大規模搜尋/驗證的另一 session（建議用 `parallel-benchmark` / `headless-research` agent；長 compute 走背景 Bash 落檔，見 CLAUDE.md §8）。
> **為什麼**：本論文主軸（個案驅動 subclonal reconstruction）目前是 **proof-of-concept + 誠實護欄**（假想樹/甲基 corroborate/大規模未驗證）。用戶要把這些升級為 validated → 需以下大規模驗證。**跑完回填真值 → 論文 R5 climax + de-confound 稀有性 claim 才能去掉 placeholder。**
> **鐵律（§13）**：每項跑完先落檔（.json/.tsv）→ Read 回真值 → 才寫進論文/報告；**勿把未測當已驗陰性、勿填預期值**。

---

## §0 一句話：四個必驗，決定論文能多強

| 驗證 | 一句話問題 | 決定論文什麼 |
|---|---|---|
| **V-1** 🔴 | 乾淨 somatic cis 真稀有，還是只測了 816？ | de-confound「稀有」claim 能否成立（headline 脈絡）|
| **V-2** 🔴 | 能否系統化找到「甲基+多點突變共同支持 subclone」的位點 + 標準化判準？ | **主軸 R5 從個案升 systematic** |
| **V-3** | chr2 主 demo 圖的數字真值？ | Fig 1（主 demo）能否寫真值 |
| **V-4** | BRCA2 單點個案有無乾淨 cis 成分（去 copy）？ | BRCA2 當單點 narrative exemplar 的誠實口徑 |

外加既有 open gap（V-5~V-9）解鎖 tier ⭐3→⭐4 與誠實定案。

---

## §1 🔴 V-1 — 乾淨 somatic cis 的真實稀有性（用戶頭號關切）

- **問題**：論文目前寫「乾淨 somatic cis = 816 可測中僅 1（chr17）」。但 cis-test **只跑了 HCC1395 的 816 個 HP-axis loci**；catalog 的 332,705 位點**從未全 cis-tested**（稽核校正 #2）。→ **「稀有」是真的生物現象，還是只是測得不夠多？**
- **現狀 tier**：🟡（under-tested，非已證實稀有）。
- **要做**：把 normal-anchored cis-test **擴大到所有合格位點 × 全 6 樣本**（合格 = 有 germline het + somatic allele + matched normal 甲基 + 足夠 read）。
  - 起手：先擴到 catalog 中所有 reliable (TAG-C, 12,868) + latent (TAG-E, 28,254) 位點；再評估能否覆蓋更大子集。
  - 每位點算 d_somatic/d_cis/d_drift（cis-candidate 判準 `p_cis<0.05 AND |d_cis|>1.8|d_drift|+0.02`）+ copy-partition（`|d_within|<0.5|d_HP|` → copy 主導）。
- **腳本**：`research/tsg_promoter_asm_reviewer/scripts/34_control_loci_cohesion_cistest.py`（cis-test）+ `37_copy_partition_confirm.py`（copy-partition）+ `03_step4_ism_methylation_diff.py`（Δβ）。需改 `--loci` 吃大位點清單、`--sample` 跑 6 樣本。
- **資料**：6 樣本 tumor+normal BAM、somatic VCF、甲基（normal 甲基 5/6 ready，COLO829 缺）；`catalog_skeleton.tsv`（332,705×32 欄）。
- **預期輸出**：每樣本 + pooled 的「乾淨 somatic cis 位點數 / 可測位點數」真實比率 + 位點清單 + 標準化判準。
- **成功判準**：(a) 若擴測後乾淨 cis 仍 <1%、跨樣本一致 → **「稀有」成立**，de-confound headline 站得住、chr17 例外更有力；(b) 若大量新乾淨 cis 出現 → **「稀有」是 under-test artifact**，論文須改述（floor 改「在 N 可測位點中 M 個」、headline 降調）。**兩種結果論文都站得住，但措辭不同——必須測。**
- **影響 claim**：de-confound RULE（R3）、headline 脈絡、chr17 例外的意義。

---

## §2 🔴 V-2 — subclone 重建：從個案 → 系統化 + 標準化判準（主軸核心）

- **問題**：主軸是「ISM 找到甲基+多點突變共同輔助區分/驗證 subclone+演化」的位點。目前是**個案（chr2 區）**，無標準化判準、無位點計數、樹是假想。→ **能否系統化找出這類位點 + 給判準 + 量化頻率？**
- **現狀 tier**：⭐3 個案 / proof-of-concept。
- **要做**：
  1. **定義「重建支持位點」標準化判準**（pre-register）：例如某窗內 ≥N 個 somatic 變異（兩 caller 一致：ClairS+DeepVariant）+ 甲基 read-clustering 與突變分群一致（ARI/PERMANOVA）+ LOH 脈絡 + 可建假想子單倍型分支（如 chr2 的 2-1/2-1-1/2-2-x）。
  2. **全基因組 × 6 樣本掃描**，套判準 → 位點清單 + 頻率。
  3. **每候選位點**：輸出 read-level 甲基+突變矩陣 + 假想樹 + 與資料一致性指標（甲基分群 vs 突變分群 ARI）。
- **腳本/工具**：ISM C++ engine（read×CpG→距離→分群→PERMANOVA→HP/LOH tag）；甲基 vs 突變分群一致性（ARI，参 `blind-ARI+imprinting 正控` 機制）；somatic 兩-caller 一致（ClairS + DeepVariant intersect）。
- **資料**：6 樣本 tagged BAM、somatic VCF（含 ClairS + DeepVariant 兩套）、甲基、LOH bed。⚠ 需確認 6 樣本是否都有 DeepVariant somatic call（圖2 用了 google=DeepVariant；若只有 HCC1395 有，須補或標範圍）。
- **預期輸出**：標準化判準文件 + 重建支持位點 catalog（位點數、每位點假想樹 + ARI 一致性）+ chr2-type pattern 頻率。
- **成功判準**：(a) 找到 ≥K 個符合判準的位點 + 跨樣本可見 → 主軸從「個案」升「systematic（有判準）」，R5 可寫 validated；(b) 仍只有零星個案 → 維持 proof-of-concept，但有了判準與頻率上界（誠實）。
- **影響 claim**：**主軸 R5 climax**（個案 → 系統化）；標題 reconstruction 的説服力。

---

## §3 V-3 / V-4 — 兩個 narrative exemplar 的數字溯源與重驗

### V-3 chr2:18M 主 demo 數字溯源
- **問題**：圖1/圖2（chr2:18,068,480–18,100,683，HCC1395）的數字（read counts、somatic A/G @18,086,020 ClairS 79:39A/40G、DeepVariant 58:29A/29G、各位點 C→G/G→C、假想分類 2-1/2-1-1/2-2-1/2-2-2/2-2-3）**目前未經本研究流程溯源**。
- **要做**：重跑 ISM + pileup 該區，確認所有數字；確認兩 caller 一致；輸出可放 Fig 1 的真值表 + 圖。
- **成功判準**：每個進 Fig 1 的數字都能 grep 到輸出檔。
- **影響**：Fig 1（主 demo）能寫真值（否則維持 {{待填}}）。

### V-4 BRCA2 單點 narrative exemplar 重驗（用戶偏好此例）
- **問題**：用戶要 BRCA2 當主要單點敘述例（甲基+haplotype+somatic 同見、知名癌症基因）。但稽核：BRCA2 = **copy-confounded**（d_within −0.023 ≪ d_copy −0.11, TAG-B），**非乾淨 cis**。
- **要做**：copy-partition + per-read residual permutation（`37_copy_partition_confirm.py` BRCA2 路徑：HP1-1 內 alt-ref vs 2000× 隨機 split null）重驗，**確定 BRCA2 有無乾淨 somatic-cis 成分**，或確認純 copy-confounded。
- **成功判準**：BRCA2 cis 成分有/無的定論 + 可放圖的 read-level 甲基+haplotype+somatic 整合數據。
- **影響**：BRCA2 當單點 exemplar 的**誠實口徑**（「整合展示」vs「乾淨 cis」——若 copy-confounded，敘述為 illustrative-integration 並標 caveat，cis 錨仍用 chr17）。

---

## §4 既有 open gap（解鎖 tier + 誠實定案）— V-5~V-9

| # | gap | 問題 | 方法 | 成功 = 解鎖 | 優先 |
|---|---|---|---|---|---|
| **V-5** | **G-A 跨樣本** | V10(normal not-copy)+V11c(亞群存在性)+catalog 跨 6 樣本重現？| 把 A0_assets 腳本 `--bam` 跑全 6 樣本（normal 甲基 5/6 ready；COLO829 補或標缺）| R3/R7 + 主軸 ⭐3→⭐4 | P0 |
| **V-6** | **G-B within-hap null** 🔴 | subclone 甲基是否 somatic-specific（非 germline-allelic）？| **正確 null = within-haplotype somatic-vs-baseline**（非 germline-het null，後者跨-allele 錯對照）| subclone-甲基 undetermined → 定案（能否寫「甲基 corroborate subclone」變 validated）| P0 |
| **V-7** | **FDR #24** 🔴 | reliable 12,868 + chr17 perm p=0.001 跨 816 的 null/FDR 校準（含 n_reads/coverage 校正）| label-shuffle null 重算 reliable rate + BH-FDR across loci + chr17 Bonferroni(×816)/genome-wide perm；**必含 n_reads 校正**（O11 教訓 epipoly AUC 0.845→0.530）| Methods 嚴謹度（投稿必補）；reliable/chr17 數字穩固 | P0 投稿前 |
| **V-8** | **HD-1 R-SELFREF** | phasing 脊柱 by-construction 循環（自我參照）？| 跑 R-SELFREF（~25-50hr C++）破循環 | R6 能否寫 positive（否則維持 ⭐3 支撐 + caveat）| P1（不押也可）|
| **V-9** | **G-C/G-D/G-E** | cis vs 突變足跡分離(±1-2kb)；亞克隆結構 demo 輸出；非 longphase-S 第二 pipeline | normal-anchored cis 寬空間控制 / 重建 demo pipeline / 第二定相工具 | R5 精確度 / Fig1 説服力 / 破 single-pipeline ⭐3 | P1-P2 |

---

## §5 資料與工具清單（執行 session 用）

- **6 樣本路徑**：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{sample}/paired_full/{date}_*_complete_matrix/longphase_s/{sample}_tagged.bam` + 同目錄 `somatic_pass.vcf.gz`；{date}：HCC1395=`20260314`、其餘=`20260315`。ISM 輸出在 `intersubmod_tp/` `intersubmod_fp/`。
- **normal 甲基**：5/6 ready（HCC1395 5mC+5hmC；HCC1937/1954/H1437/H2009 5mC；COLO829 缺 MM tag）；契約 `docs/data_specs/20260612_external_data_dependencies_01.md`；⚠ 6 normal 全 `zhenyu112` 帳號 = SPOF。
- **腳本**：cis-test/copy-partition/Δβ = `research/tsg_promoter_asm_reviewer/scripts/{34,37,03}*.py` + `pipeline/lib/cis_asm_core.py`；catalog = scripts 98-100 + `catalog_skeleton.tsv`；phasing-assist/V1-V12 = `docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/*.py`（可 `--bam` 換樣本）；ISM C++ = `src/core/`。
- **真值 SoT**：`research/autoresearch/evidence_ledger.jsonl` + `VERIFIED_RESULTS.md` + 共識台帳 §1。
- **資源 preflight**：6 樣本 × cis-test 可能重；長 compute **背景 Bash 落檔**（勿放 workflow agent step），≥3 樣本實跑用 `parallel-benchmark` 但先查 CPU/mem（CLAUDE.md §8 長 job 鐵則）。

---

## §6 驗證完成後，論文如何上修（誠實護欄解除條件）

| 護欄（現） | 解除條件 |
|---|---|
| 假想樹標 hypothetical | V-2 標準化判準 + ARI 一致性達標 → 改「依判準重建（criteria-based）」，仍標單樣本/未獨立驗證者 caveat |
| 甲基寫 corroborate 非 validated | **V-6 within-hap null 通過** → 才可寫「甲基（去混淆後）支持 subclone」；否則維持 corroborate |
| 乾淨 cis「816 中 1」| **V-1 擴測** → 改真實分母（「N 可測中 M」）+ 跨樣本率 |
| 大規模驗證 = future | V-2 完成 → 改「systematic（N 位點，判準如下）」 |
| single-sample ⭐3 | **V-5 G-A 6 樣本** → ⭐4 |
| BRCA2 copy caveat | **V-4 重驗** → 依結果定 illustrative vs 部分 cis |

> **回填紀律**：每項回填時，更新 `VERIFIED_RESULTS.md` / 共識台帳 / 對齊表 + 論文對應節，並標 build_commit。**勿在未跑完前先改論文措辭**。
