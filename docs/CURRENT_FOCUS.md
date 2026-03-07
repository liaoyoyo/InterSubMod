<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-03-08 23:00
狀態: validated
資料來源:
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - scripts/analysis/check_ai_agent_readiness.sh
-->

# 當前目標

## 1. 目前狀態

1. docs 重整已完成核心落地：
   - 命名：`YYYYMMDD_主題_流水號.md`
   - 報告分層：`reports/validated|finalized`
   - 實驗分層：`experiments/in_progress|validated|finalized`
2. `output/` 入口已固定為軟連結：
   - `output -> /big8_disk/liaoyoyo2001/InterSubMod_runs/output`
3. Knowledge MCP 已接入：
   - `.mcp.json` 指向 `/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py`

## 2. AI Agent 主要入口

1. docs 導航：`docs/README.md`
2. 研究歷史索引：`docs/experiments/INDEX.md`（已試驗方向、成功/失敗結論、建議後續）
3. Agent 手冊：`docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md`
4. 健康檢查：`scripts/analysis/check_ai_agent_readiness.sh`
5. 文件規範：`docs/standards/README.md`

## 3. 當前進行中

1. 研究主線已重排為「純樣本優先」：
   - 第一主實驗：`HCC1395 5kHz` paired tumor/normal
   - 第一交叉驗證：`HCC1395 ONT_Dorado + HCC1395BL`
   - 後續順序：其他 pure samples → pure tumor-only → mixed purity/subsample
2. 當前最重要的研究問題：
   - `label-first` 與 `cluster-first` 的設計本身是否合理、可用、可解釋
   - 加入 HP/Allele 後，哪些 region 會從 `Weak` 升級為 `Strong` 或 `Subclone`
   - 哪些流程節點是真正影響結果的因子，哪些只是暫時工程預設
3. `2026-03-07` 新版完整欄位 paired pure 正式 rerun 已完成：
   - `HCC1395 5kHz`：`LongPhase-S F1=0.8522` → `InterSubMod F1=0.8532`（`+0.0010`）
   - `HCC1395_DORADO`：`LongPhase-S F1=0.8592` → `InterSubMod F1=0.8590`（`-0.0002`）
   - 跨 7 個 pure paired 樣本比較後，`低 VAF + 偏高 AlleleDelta + 不穩定全域結構` 只在 `HCC1395 5kHz` 顯示明顯正增益
   - `Noise->Strong` 也幾乎只在 `HCC1395` 系列明顯出現，尤其 `5kHz` 最值得警戒
   - 新增 refined label 分析後，`Strong + 低 VAF + 高 AlleleDelta` 在 `HCC1395 5kHz` 可移除 `68` 個 FP、僅誤傷 `1` 個 TP（`F1 +0.000806`），但 `DORADO` 與其他樣本不支持將其升級為全域規則
   - `Noise->Strong` 本身不適合直接刪除，因為在 `5kHz` 會移除 `1659` 個 TP、僅移除 `51` 個 FP
   - 反向救援分析顯示：在 `HCC1395 5kHz`，`LongPhase-S` 丟掉的 caller TP 中，有 `45` 筆可用單純 `GQ>=20` 救回且不重新引入精確對齊集合中的 FP（`F1 +0.000739`）
   - 目前的 InterSubMod `VerificationClass / PairwiseDist` 對 `LongPhase-S` 的 TP rescue 沒有超過 caller `GQ` 的表現，甲基訊號現階段更適合用於 FP triage 與 Strong 細分類
4. 研究流程與紀錄制度化：
   - 研究藍圖：`docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md`
   - 執行計畫：`docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md`
   - 研究手冊：`docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md`
   - 最新進度報告：`docs/experiments/in_progress/2026/03/20260307_純樣本round執行與進度報告_01.md`
   - 候選規則與 shift 跨樣本分析：`docs/experiments/in_progress/2026/03/20260307_低VAF高AlleleDelta與shift跨樣本分析_01.md`
   - Strong 細分類與 haplotagged samtools 驗證：`docs/experiments/in_progress/2026/03/20260307_Strong細分類與samtools驗證分析_01.md`
   - LongPhase 救回 TP 與甲基救援分析：`docs/experiments/in_progress/2026/03/20260308_LongPhase救回TP與甲基救援分析_01.md`
   - LongPhase-TO 空間需求與中間產物確認：`docs/experiments/in_progress/2026/03/20260307_LongPhaseTO_空間需求與中間產物確認_01.md`
   - HCC1395 5kHz TO pilot 啟動與執行紀錄：`docs/experiments/in_progress/2026/03/20260307_HCC1395_5kHz_TO_pilot啟動與執行紀錄_01.md`
5. `2026-03-07` 新確認：
   - 真正要做 `LongPhase-TO + InterSubMod` 驗證時，仍需 `haplotagged BAM`；目前 repo pipeline `02_intersubmod.sh` 直接吃 `--tagged-bam`
   - `benchmark_split_snv_vcf.sh` 已補上 plain `.vcf` 自動 `bgzip + tabix`；`run_longphase_to_intersubmod_pilot.sh` 已補上 `--tag-prefix`
   - `20260307_hcc1395_to_pilot_1` 已完成：
     - `ClairS-TO`
     - `LongPhase-TO phase`
     - `LongPhase-TO haplotag`
     - `InterSubMod`
     - benchmark / design validation / round summary
   - `HCC1395 5kHz tumor-only` 本輪結果：
     - `ClairS-TO F1=0.7127`
     - `LongPhase-TO F1=0.7127`
     - `InterSubMod F1=0.7130`
     - `InterSubMod` 移除 `36` 個 FP、誤傷 `3` 個 TP，`F1 delta=+0.0003`
   - 在這個 `5kHz TO` pilot 上，`LongPhase-TO phase` 本身沒有改變 benchmark call set；目前正增益來自後續甲基過濾層
   - TO 場景下 `label-first` 仍非空：
     - `label_upgrade=2756`
     - `conflict=3227`
     - `consistent_strong=1471`
   - `5kHz TO tagged BAM` 實測大小為 `260G`，`.bai` 約 `113M`
   - 截至本輪主流程完成，磁碟大約剩：
     - `/big8_disk` `228G`
     - `/bip8_disk` `207G`
   - 因此目前**不需要立刻手動清理**，但在進入下一個 `200G+` 級 tagged BAM 任務前，應先將 `5kHz TO tagged BAM` 搬到 `Archive_pending_delete/`
   - `HCC1395_DORADO` 已補齊 TO baseline benchmark（`F1=0.7226`），仍是下一個最合理的第二階段驗證樣本
   - `2026-03-08` TO 特徵回查已確認：
     - `label_upgrade=2756` 中 `TP=2526`、`FP=230`，`FP rate=8.35%`
     - `conflict=3227` 中 `TP=2465`、`FP=762`，`FP rate=23.61%`
     - `Weak->Strong` 與 `Noise->Strong` 在 TO 下都以 `TP` 為主，不適合直接當 artifact 刪除規則
     - TO 下真正仍有效的 `FP` 特徵是 `low VAF + high AlleleDelta`，而且這批幾乎不落在 `label_upgrade/conflict`
     - 該特徵在 TO 下主要落在 `cluster_only + Strong->Weak + 低 GQ/QUAL`
     - `low VAF + high AlleleDelta + low CramersV` 在 `5kHz TO` 可移除 `36 FP`、僅誤傷 `2 TP`，`F1 +0.000261`
   - `2026-03-08` TO diagnostics 與 scheme 調整新確認：
     - 本輪 `5kHz TO` 是 `ClairS-TO pileup` 路線，不是 full model；目前結論先只對 pileup 成立
     - top `Strong->Weak + lowVAF/highAD` 可疑 FP 中有 `4/5` 聚集在 `chr11` 的 `1402 bp` 區段
     - `HP skew` 本身不足以區分 FP/TP；TP 對照也可呈現極端單一 HP
     - 比較有辨識力的是 `高 AlleleDelta + 低 PairwiseMedianDist + 很低 alt fraction + 低 QUAL/GQ`
     - 現行 core 的 `Strong` 會混入「只有 allele sig、沒有 HP sig」的位點；這批不適合整包刪除（`5104 TP / 2042 FP`）
     - 新增 `cluster_plus_weak_label` 後，原先的 TO 可疑子集更精確地落在：
       - `cluster_plus_weak_label + Strong->Weak + lowVAF/highAD`
     - 分類語意調整目前主要提升的是可解釋性，尚未超過既有 `low VAF + high AlleleDelta + low CramersV` 的 F1 增益
   - `2026-03-08` ClairS 邊緣 TP rescue 與甲基輔助評估第一輪已完成（**只用 final VCF，不用 pileup candidate/tensor**）：
     - `HCC1395 5kHz paired`：
       - `caller_lost_tp=1052`（`candidate_eligible=920`）
       - `caller_removed_fp=12974`（`candidate_eligible=5266`）
       - 依規劃排序的最佳 safe rule：`candidate_gq_ge_15`，`106 TP / 75 FP`，`F1 +0.000825`
       - 若只看最佳 `delta F1`：`gq>=20`，`59 TP / 8 FP`，`F1 +0.000871`
       - 甲基 summary 目前只覆蓋 `111/1052` lost TP 與 `802/12974` removed FP；本輪 `methylation support/veto` 尚未超過 caller-only
     - `HCC1395_DORADO paired`：
       - `caller_lost_tp=1489`（`candidate_eligible=1122`）
       - `caller_removed_fp=2533`（`candidate_eligible=1658`）
       - 依規劃排序的最佳 safe rule：`candidate_gq_ge_15`，`97 TP / 88 FP`，`F1 +0.000502`
       - 若只看最佳 `delta F1`：`gq>=20`，`53 TP / 7 FP`，`F1 +0.000782`
       - 尚缺 `DORADO paired` 的 candidate-specific InterSubMod summary，因此本輪不能把甲基 rescue 解讀成有效 negative result
     - `HCC1395 5kHz TO`：
       - `caller_lost_tp=2108`（`candidate_eligible=773`）
       - `caller_removed_fp=2720008`（`candidate_eligible=298`，`pon_hit_fraction=0.999019`）
       - 最佳 caller-only safe rule：`qual>=10 or gq>=20`，`491 TP / 118 FP`，`F1 +0.006824`
       - 第二名：`gq>=15`，`365 TP / 71 FP`，`F1 +0.005233`
       - 這輪 `TO` 的 downstream lost/removed 與現有 baseline InterSubMod summary overlap 為 `0/0`，因此 **TO 甲基 rescue 尚未被正確測到**
     - 本輪正式報告：
       - `docs/experiments/in_progress/2026/03/20260308_ClairS邊緣TP_rescue與甲基輔助評估_01.md`
   - `2026-03-08` ClairS 邊緣 FN 探勘（Pool B）與甲基救援分析已完成（**實驗 25**）：
     - `Pool B FN`（ClairS non-PASS ∩ truth set）：**840 個**
     - 策略 B 補跑後甲基覆蓋率：**99.9%（839/840）**
     - **最佳 Caller-only**：`no_varcluster`（F1 delta=+0.003391）；最乾淨：`no_varcluster_and_gq15`（precision=71.9%，F1 delta=+0.001452）
     - **最佳 Combined（Caller+Methyl）**：`gq15_and_allele_delta_low`（F1 delta=+0.002702），僅比 Caller-only 好 +0.000386
     - **甲基唯一規則全數為負**：Pool B non-FN 含大量 germline 變異，Strong rate 比 FN 更高（44.7% vs 30.8%）
     - **AlleleDelta 具 somatic 特異性**：FN 中位 AD=0.010，non-FN 中位 AD=0.065，但 precision 仍只 ~32%
     - **結論**：甲基訊號不適合單獨用於 Pool B rescue；現行標準無需調整，問題在 Pool B non-FN 本身混入 germline
     - **不推薦升為正式規則**；`no_varcluster_and_gq15` 建議在 `HCC1395_DORADO TO` 做交叉驗證
     - 正式報告：`docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md`
6. docs 維護工作仍持續進行，但現已降為研究主線的次要維護事項：
   - 降低舊語彙殘留（例如舊輸出路徑與舊文件分層名稱）
   - 收斂腳本路徑設定（優先使用可配置 `OUTPUT_ROOT` 或 `output/` 入口）
   - 清理不必要空目錄，降低 Agent 探索噪音

## 4. 阻塞與風險

1. `archive/deep/` 為 immutable 快照，保留歷史失效連結（不回寫）。
2. 活躍腳本仍有部分硬編碼歷史路徑，可能影響跨專案重用。
3. `purity-aware` / subsample 結果容易被 tumor-normal 組織甲基差異混淆，暫不作主證據。
4. `label-first` 與 `cluster-first` 雖已有輸出，但方法學正當性與適用邊界仍需驗證。
5. 若只報 F1 而不附 `truth_total/TP/FP/FN`，容易高估方法改善幅度或誤判結論。
6. `pileup` 與 `full model` 的 TO 輸入若混用，會直接改變 caller VCF 與 downstream 結論；後續報告必須明確標示來源。
7. `TO` 的 borderline rescue 若直接拿 baseline kept-set 的 InterSubMod summary 來做 join，會得到錯誤結論；必須使用 candidate-specific `lost_tp / removed_fp` run。

## 5. 每次研究啟動必查

1. `docs/CURRENT_FOCUS.md`
2. `docs/experiments/INDEX.md`
3. `docs/README.md`
4. `docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md`
5. `docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md`
