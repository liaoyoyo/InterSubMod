<!--
建立時間: 2026-08-11 23:59 CST
目標: 依現行 InterSubMod 研究與機器可重算數據，完整稽核 CCU Bioinformatics Lab lab-tutorial 公開網站的科學敘述、數字口徑與證據上限
處理範圍: 25 個 live HTTP 200 頁面、1 個 source-only/HTTP 404 頁面、網站 commit 46b6f5b3016c187ad742fecbfa813f835b09e605、7 technical datasets/6 biological IDs/chr1-22 exact-PS authority、paired 與 tumour-only 研究證據、LongPhase-S/TO KB 與 runtime/code contract
關聯檔案: InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/00_INDEX.md; InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/route_status.tsv; InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/key_number_recheck.tsv; InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/page_verdict_matrix.tsv
branch: feat/lineage-tag-methylation-axes
commit: 73afaeac8e61c767241fa59c1ca6043a1c95290c
website_commit: 46b6f5b3016c187ad742fecbfa813f835b09e605
worktree: dirty; 未覆寫既有未追蹤檔案
status: VALIDATED_CLAIM_AUDIT
report_class: B_COMPREHENSIVE_VALIDATION
framework: SCQA + Claim-Evidence-Verdict
-->

# CCU lab-tutorial 網站科學敘述完整稽核

用 SCQA + Claim–Evidence–Verdict：**網站的基礎教學大致合理，但整體目前只能「有條件通過、需重大修訂」；首頁、M9、M12、M13、SR1、SR2b 的分子→細胞升格或舊結論會實質誤導，M6/M7 另有可由程式碼反駁的工具事實錯誤（影響：高，信心：高）。**

> **整體判定：CONDITIONAL PASS — MAJOR REVISION REQUIRED**  
> 基礎概念可以保留；在 P0 頁面修正前，不宜把網站當作 InterSubMod 目前科學結論的權威入口。

## 1. 30 秒結論

### Situation

網站把 bulk long-read、phasing、somatic haplotagging、purity、甲基化與 local topology 串成一條容易理解的教學線。M0、M1、M4、M8、M10 與 capstone 的核心框架多數合理，尤其「read 是 DNA 分子，不是細胞」及「tumour DNA fraction 不等於 cellular purity」兩個 guardrail 是正確的。

### Complication

同一網站的首頁、M9、M12、M13 與 SR1 卻把 read/haplotype/local molecular state 升格為 cell/cell group/clone lineage；M13 又宣稱特定 synthetic setup 中 DNA fraction 與 cellular purity 恰好相等，直接和 glossary 的公式衝突。現行 authority 明訂：CN/LOH 未整合、甲基僅 association-only、confirmed cellular subclone=0、linear biological ancestry=0。

### Question

網站能否依目前研究與數據被視為合理、哪些敘述可被反駁、哪些數字應補充？

### Answer

- **可保留**：long-read 分子連鎖、HP/PS 基礎、paired 與 tumour-only 的問題差異、VAF/purity/CN/CCF 的概念拆分、評估與 leakage guardrails。
- **需立即改**：任何「read 屬於某癌細胞群」「已重建 clone tree」「局部 state 就是 cellular lineage」「DNA fraction 就是癌細胞比例」的句子。
- **數字正確但解釋過強**：M9 的 51.2%–80.4%、9.5%、35.4% 可由 M1 screen 重算，不能解讀成 cell-group prevalence 或正式 subclone 證據。
- **有直接反證**：TO 特徵的跨樣本安全過濾沒有泛化；0.909 是技術再現性而不是準確度；64.89% 找不到原始出處；LongPhase-TO 不只支援 ONT，且其 haplotag 是 `HP:i`，不是全程 `HP:Z`。
- **可補充的正向結果**：paired-pure read classifier 有 ΔF1=+0.01116，但這是 73,230 sampled reads／537 regions 的 read-level PoC，不是 whole-genome variant-call F1。

## 2. 稽核契約與可重現範圍

| 項目 | 固定值 | 驗證結果 |
|---|---|---|
| Task type | B — Comprehensive validation | 不採 subset |
| 服務目標 | G2、G3、G4、G5 | paired/TO 分離、甲基 claim ceiling、多樣本重現、外部可驗證性 |
| 網站 snapshot | commit `46b6f5b3016c187ad742fecbfa813f835b09e605` | live 內容 hash 與該 commit 的 `site/*.html` 相符 |
| Live coverage | 25 routes | 25/25 HTTP 200 |
| Source-only | `sr6.html` | GitHub source 存在；live HTTP 404，且不在導覽 |
| Canonical research scope | 7 technical datasets、6 biological IDs、chr1–22 | HCC1395 與 HCC1395_DORADO 是同一 biological cell line |
| Authority integrity | 13 artifacts | SHA-256 13/13 MATCH |
| Primary analysis grain | dataset × chromosome × exact PS × primary HP × strict read-linked component × bounded local block | 不是 cell、clone 或 whole-sample phylogeny |

完整 route receipt：[`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/route_status.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/route_status.tsv)。網站 snapshot 可由 [GitHub commit 46b6f5b](https://github.com/CCU-Bioinformatics-Lab/lab-tutorial/commit/46b6f5b3016c187ad742fecbfa813f835b09e605) 重現。

### 證據等級

| Tier | 本報告定義 | 可支持的語氣 |
|---|---|---|
| L1 | machine artifact、程式碼/runtime、hash、原始 numerator/denominator | 可重算事實；仍受 scope 限制 |
| L2 | 跨資料集、LOSO、技術重複或獨立重算 | 可談再現／泛化，但不自動等於生物真值 |
| L3 | validated report、audit card、已審查整合 | 可作整合結論，不能凌駕 L1 |
| L4 | spec、toy model、方法假說 | 只能說「設計／預期／候選」 |
| L5 | 教學直覺、未附來源精確數字 | 不可當驗證證據 |

## 3. 現行 claim ceiling：網站最需要補上的一道牆

目前資料能支持的因果鏈只到「局部分子狀態／候選數學拓撲」：

```text
bulk BAM/VCF/MM/ML
        ↓ 直接觀測
read-level allele / exact-PS / HP / methylation pattern
        ↓ 有條件推論
bounded local molecular state + model-conditional candidate graph
        ╫ 仍缺 CN/LOH/purity/multiplicity/CCF、跨 block bridging、正交 biological truth
confirmed cellular clone / lineage / whole-sample phylogeny
```

Authority 明列 methyl role 是 `pattern-conditioned association-only sidecar`，禁止升格為 confirmed cellular clone/subclone、whole-sample cancer phylogeny、methylation-rescored topology 或 cross-HP cellular pairing。來源：[`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`](/big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json)。

因此，下列兩句看似相近，證據地位完全不同：

- 合理：**「這批 reads 在同一 exact-PS/HP/local block 中呈現不同的分子狀態。」**
- 目前不合理：**「這批 reads 分別來自兩群癌細胞，且兩群具有某一祖先關係。」**

## 4. 必須優先修正的 P0 敘述

| 頁面 | 網站目前傳達的意思 | 現行證據 | Verdict | 建議改法 |
|---|---|---|---|---|
| 首頁 | 軟體把每條字串還原到正常／癌細胞來源，並算癌細胞比例 | read 可指派到 haplotype/molecular carrier state；LongPhase-S 輸出欄位雖名 purity，實質 estimand 是 tumour DNA fraction | **被 claim ceiling 反駁** | 改為「指派到局部 haplotype 狀態、辨識 somatic-supporting molecules、估計 tumour DNA fraction」 |
| M9 | 同一群細胞甲基較像；ALT reads 全來自同一 cell group | 正式結果只有 3/811 robust regional associations；沒有 cellular truth，且 truncal ALT 可跨多 clone | **推論不成立** | 改為「固定 genetic pattern 後的 regional methylation association，可作 bounded auxiliary annotation」 |
| M12 | somatic haplotagging 可看出哪些 reads 屬於同一癌細胞群 | LongPhase-S 的 tag 是 germline/somatic haplotype-carrier 類別，不是 cellular clone ID | **層級錯置** | 改為「哪些 reads 支持同一局部 germline/somatic haplotype state」 |
| M13 | 真變異必然形成一條只差候選一格的新路徑；GT3 兩支代表不同 cell groups；synthetic 中 f=p | 路徑是理想化假設；local state≠cell；在模型 `f=pκ/[pκ+2(1−p)]` 中只有 κ=2 才有 f=p | **多項直接錯誤** | 加上模型條件與反例；移除「必然／細胞群」；用 glossary 的換算式 |
| SR1 | 目前工作已由 long reads 重建 clone tree，甲基可當細胞 clock | 現行輸出是 local candidate molecular topology；甲基 association-only | **研究願景被誤寫成已完成** | 分成「長期目標」與「目前 validated output」兩欄 |
| SR2b | 71,955 vs 85,941 尚待釐清；0.909 是最強真訊號 | 分母是已解釋的巢狀 funnel；0.909 是 technical reproducibility，不是 accuracy | **已被 SR2c 與 authority 更新** | 加 OBSOLETE banner/redirect 或同步重寫 |

來源頁面：[首頁](https://ccu-bioinformatics-lab.github.io/lab-tutorial/index.html)、[M9](https://ccu-bioinformatics-lab.github.io/lab-tutorial/m09.html)、[M12](https://ccu-bioinformatics-lab.github.io/lab-tutorial/m12.html)、[M13](https://ccu-bioinformatics-lab.github.io/lab-tutorial/m13.html)、[SR1](https://ccu-bioinformatics-lab.github.io/lab-tutorial/sr1.html)、[SR2b](https://ccu-bioinformatics-lab.github.io/lab-tutorial/sr2b.html)、[SR2c](https://ccu-bioinformatics-lab.github.io/lab-tutorial/sr2c.html)。

### 一個站內即可完成的反證：f 不等於 p

M13 的「coverage-mixed synthetic 樣本中兩者恰好相等」與同版 glossary 自相矛盾。若 tumour DNA fraction `f=0.2`、tumour ploidy `κ=3.2`：

\[
p=\frac{2f}{\kappa-f(\kappa-2)}
=\frac{0.4}{3.2-0.24}
=0.135135\ldots
\]

即 cellular purity 約 13.5%，不是 20%。只有在 κ=2 且模型其他假設成立時，才有 f=p。這不是措辭偏好，而是數學上的直接矛盾。

## 5. 數字完整重算：哪些支持、哪些反駁、哪些要補分母

完整 machine-readable receipt：[`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/key_number_recheck.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/key_number_recheck.tsv)。

### 5.1 Exact-PS topology：88.26% 有數據，但不能稱生物拓撲率

| Funnel | Numerator / denominator | 比例 | 正確解讀 |
|---|---:|---:|---|
| PASS biallelic autosomal sSNV dataset-records | 469,849 / 469,849 | 100% | 7 technical datasets／6 biological IDs／chr1–22 |
| strict components | 255,752 | — | endpoint read-linkage components |
| k=1 linkage abstain | 170,131 / 255,752 | 66.5219% | 無法做 multi-sSNV tree inference |
| mutation-bearing units | 85,941 / 98,955 | 86.8486% | 至少一個 active ALT；仍不是 clone |
| family complete | 75,224 / 85,941 | 87.5298% | search family 完整 |
| ranked complete | 71,955 / 75,224 | 95.6543% | read-AF usable 且 recurrence screen 通過 |
| unique best | 39,648 / 71,955 | 55.1011% | exact arithmetic unique representative |
| tied, same rooted-unlabeled shape | 23,858 / 71,955 | 33.1568% | 多棵 labeled tree，但同一 shape |
| tied across shapes | 8,449 / 71,955 | 11.7421% | shape 仍未決 |
| one rooted-unlabeled shape | 63,506 / 71,955 | **88.2579%** | ranked-only、model-conditional mathematical shape |

以上 funnel 由 [`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv) 重算。不同樣本的 one-shape rate 從 H2009 78.2082% 到 HCC1937 99.3875%，也不應只報 pooled 88.26%。

**直接修正 SR2c**：88.26% 可保留但必須附上述限定；「64.89%–88.26%」的 64.89% 在現行 repo 與 authority 中找不到原始出處，應刪除或標成 unverified，不能用兩個不同定義的數字組成範圍。

### 5.2 M9 methylation：數字可重現，生物解釋不可升格

網站的 51.2%–80.4%、COLO829 9.5%、H2009 35.4% 可在 M1 artifact 中重算：

| Dataset / quantity | 重算 | 結果 |
|---|---:|---:|
| HCC1395_DORADO residual / M1 | 7,569 / 14,789 | 51.1799% |
| H2009 residual / M1 | 43,937 / 54,644 | 80.4059% |
| COLO829 M1 / all sSNV | 3,579 / 37,788 | 9.4713% |
| H2009 M1 / all sSNV | 54,644 / 154,465 | 35.3763% |
| all-dataset M1 / all sSNV | 102,842 / 469,849 | 21.8883% |

但 M1 是 **operational stable-null multigroup screen**。`residual` 只表示已測軸沒有解釋該結構，不表示所有 confound 被排除，更不表示 subclone prevalence。

Formal multi-sSNV assay 問的是另一個、更嚴格的問題：

| Formal status | Count | 解讀 |
|---|---:|---|
| formal units | 1,045 | frozen pattern-conditioned units |
| evaluable | 811 | 77.6077% 通過 eligibility |
| robust association | 3 | **3/811=0.3699%**，association-only |
| evaluable, no robust | 627 | 沒有通過 robust claim |
| confounded gate failure | 181 | 至少一個 geometry/dispersion/robustness gate 失敗 |
| not evaluable | 234 | support/common CpG/exchangeability 不足 |
| full RR/RA/AR/AA robust | 0 | 沒有完整四狀態 robust unit |

六個 biological cell lines 都有 M1 locus screen；formal assay 只有五株進入，COLO829 formal units=0。這個 0 的意思是沒有單元通過完整 R/A、read 數與 state-support gate，不是「COLO829 沒有甲基資料」或「沒有甲基異質性」。來源：[`InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md)。

**因此 M9 正確寫法應是**：「6/6 cell lines 有 M1 operational screen；5/6 有 formal multi-sSNV unit，COLO829 formal units=0；3/811 evaluable units 有 robust regional association。這些結果不確認 cellular clone 或 ancestry。」

### 5.3 F1：必須區分 read classifier 與 variant caller

| 分析層級 | Baseline | 加入 InterSubMod/methyl context | ΔF1 | 可下結論 |
|---|---:|---:|---:|---|
| paired-pure sampled multi-bio read classifier，external | 0.861066 | 0.872226 | **+0.011159** | 73,230 reads、537 regions 的 PoC；95% bootstrap CI +0.004437 到 +0.018808 |
| canonical HCC1395 paired variant-level | 0.852208 | 0.853190 | **+0.000981** | 此 canonical run 的 variant-call 小幅變化 |
| canonical HCC1395 tumour-only variant-level | 0.712697 | 0.712971 | **+0.000274** | 此 canonical run 幾乎無變化 |

paired-pure +0.0112 的 machine source 是 [`incremental_comparison.tsv`](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1/incremental_comparison.tsv)。它不是 HCC1395 whole-genome variant F1，也不是普遍「甲基讓 caller 提升 1.1 個百分點」。

Tumour-only 的跨樣本反證更強：

- read-level LOSO AUC=0.721，但在 TP loss≤2% 的安全約束下，FP removal=0%；
- TO pileup LR 的 HCC1395 in-distribution ΔF1=+0.02236，held-out=-0.00012；
- 5-sample LOSO mean ΔF1=-0.00004，0 positive／4 negative／1 approximately zero，best threshold 退化為 keep-all；
- 全量觀察沒有單一 TO 特徵 AUC>0.58，germline-FP VCF feature 也沒有超過 0.64。

這些結果**不能直接推翻 LongPhase-TO 論文／特定 benchmark 的 paper-scope 結果**，但足以反駁教材若暗示「haplotype/methylation rescue 是可一般化且安全的 TO filtering rule」。

### 5.4 0.909：支持再現性，不支持正確性

HCC1395 與 HCC1395_DORADO technical pair：

- callset Jaccard=0.927646；
- VAF concordance correlation coefficient=0.933922；
- VAF Spearman=0.9093；
- 但獨立 PyClone fits 的 subclonal Jaccard=0.3810，state kappa=0.544，cluster structure 並未同等穩定。

所以 SR2b 所稱「0.909 代表訊號是真的」不成立。它只能支持 measurement/process reproducibility；若要談 accuracy，需要 single-cell、multi-region、synthetic truth 或其他獨立 ground truth。

### 5.5 見樹也見林：同時展示 aggregate、canonical、outlier 與 well-explained case

| 展示層 | Case | 數值／結果 | 為何需要 |
|---|---|---|---|
| Aggregate | pooled ranked topology units | 63,506/71,955=88.2579% | 呈現整體 model-conditional shape rate |
| Canonical | HCC1395 current variant-level | paired ΔF1=+0.000981；TO ΔF1=+0.000274 | 避免把 read-classifier +0.0112 誤植為 canonical variant F1 |
| Extreme outlier | sample one-shape range | H2009 78.2082%；HCC1937 99.3875% | 顯示 pooled 88.26% 掩蓋樣本異質性 |
| Well-explained | COLO829 formal methyl units | 0 | 由 eligibility gates 解釋；不是沒有 methylation data |

## 6. M6/M7 的可機械反駁工具錯誤

這部分依 P-14 先查工具 KB，再以 local code/runtime contract 核對。

| 網站敘述 | L1/L2 核對 | Verdict |
|---|---|---|
| LongPhase-S 有 `HP1-1-1` | 完整字串集合是 `. / 1 / 2 / 3 / 4 / 1-1 / 2-1 / 1-2 / 2-2` | **錯；移除不存在的 tag** |
| LongPhase-TO 只支援 ONT | `longphase-to phase --help` 要求 `--ont` 或 `--pb` 二擇一 | **錯；同時支援 ONT/PacBio route** |
| LongPhase-TO 全程使用 `HP:Z` | haplotagged BAM 使用 integer `HP:i:1/2`，sub-state 可能為 11/21/33 | **錯；不要與 LongPhase-S 的字串 tag 混用** |

KB 來源：[`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md)、[`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md)。外部工具精確行為仍應在頁面釘 release/commit；不能把這次 local revision 的 runtime 永久泛化到所有版本。

## 7. 逐頁完整判定

Machine-readable matrix：[`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/page_verdict_matrix.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/page_verdict_matrix.tsv)。

| 頁面 | 判定 | 優先級 | 核心理由 |
|---|---|---:|---|
| index | REVISE | P0 | 把分子來源與細胞來源、DNA fraction 與 cellular purity 混用 |
| glossary | CLARIFY | P1 | 多數正確；triplet/read-AF 定義需模型限定，且應統一約束 M13 |
| capstone | KEEP WITH PROVENANCE | P2 | 「read≠cell」guardrail 正確；精確案例需標真實或 DEMO 與來源 |
| print-all | INHERITS | P0–P2 | 是其他頁複本；修正來源頁後再重建 |
| M0 | KEEP | P2 | 分子觀測模型合理；lineage reconstruction 要標長期目標 |
| M1 | KEEP | P2 | bulk/within-read linkage 說明合理；保留 toy assumptions |
| M2 | CLARIFY | P1 | 0.909 應稱 reproducibility；truth-set provenance 太粗 |
| M3 | ADD PROVENANCE | P1 | F1 表缺 dataset/version/truth region/caller mode |
| M4 | KEEP WITH VERSION PIN | P2 | 格式概念合理；命令與 tag behavior 需釘版本 |
| M5 | CLARIFY | P1 | 是理想化 intuition；需補 systematics、LOH/CN、mapping bias、self-phasing |
| M6 | CORRECT TECHNICAL | P1 | `HP1-1-1` 不存在；平台支援說法錯；需 circularity guardrail |
| M7 | CORRECT TECHNICAL/SCOPE | P1 | TO 的 HP type 與 rescue 泛化都需修正 |
| M8 | KEEP WITH PAPER SCOPE | P2 | VAF/purity/CN/CCF 關係合理；效能只屬 cited preprint/version |
| M9 | REVISE | P0 | 數字可重算，但 screen/formal eligibility 混用且 cell-group 推論不成立 |
| M10 | KEEP WITH PROVENANCE | P2 | 方法論合理；精確百分比與模型數值缺來源 |
| M11 | VERSION PIN | P1 | 工具 internals 高度易過期，需 commit/release/runtime help |
| M12 | REVISE | P0 | haplotype-carrier classes 不能寫成同一癌細胞群 |
| M13 | REVISE | P0 | 「必然」路徑、cell group、f=p 均超出證據或直接錯誤 |
| SR1 | REVISE CLAIM CEILING | P0 | 文獻背景可留；目前實作不能稱已重建 cellular clone tree |
| SR2 | SPEC / UNVALIDATED | P1 | 方法方向可留；首屏需持續標規格而非結果 |
| SR2b | OBSOLETE / REVISE | P0 | 分母與 0.909 已被 SR2c/authority 修正 |
| SR2c | KEEP WITH CORRECTION | P1 | 最接近現行 local molecular ceiling；64.89% 無來源，仍有 cell wording |
| SR3 | KEEP WITH PAPER SCOPE | P2 | purity/CN 綜述大致嚴謹；外部常數/軟體需版本限定 |
| SR4 | SPEC / UNVALIDATED | P1 | 已承認不是 estimator；再補 self-phasing、CN/LOH independence |
| SR5 | SPEC / UNVALIDATED | P1 | 已明示未跑 real/simulated data，應持續顯示 |
| SR6 | SOURCE-ONLY SPEC | P1 | source 存在但 live 404；部署狀態與 CN/LOH 未整合需明示 |

## 8. 哪些外部／paper-scope 敘述本輪不能直接判錯

本 repo 的 negative result 不等於外部工具論文失效。以下項目缺同 scope、同版本或完整 numerator/denominator，因此本輪判「需 provenance／限縮」，不判 false：

- M3 的 ClairS／DeepSomatic 精確 F1 表；
- M8、M12 的 LongPhase-S MAE、R²、F1/recall 增益；
- M10 的 4%/19%、0.997/0.994、λ≈1.45、epoch 12/22；
- capstone 的 HCC1395 精確座標與 confusion counts；
- SR1 的「methylation clock 快五個數量級」與固定 20 kb／5 variants per Mb toy assumptions；
- LongPhase-TO HCC1395 LOH bp-level F1≈0.9625：本輪可重算，但只適用單一 cell line、固定 reference/parameters，不能外推六株。

建議每個精確數字至少附六欄：

```text
source/version | dataset | scope/region | numerator/denominator | metric definition | evidence tier
```

可作 paper/version 入口的 primary sources：[LongPhase-S preprint DOI](https://doi.org/10.1101/2025.11.20.689492)、[ClairS, Nature Methods 2026](https://www.nature.com/articles/s41592-026-03152-4)。引用 LongPhase-S 效能時應以 Results 的具體 metric 為準，不把摘要中的 recall 增幅誤寫成 F1 增幅。

## 9. 建議直接替換的文字

### 首頁

> 我們把 long reads 指派到局部 germline/somatic haplotype state，辨識支持 somatic variants 的 DNA molecules，並在明示模型假設下估計 tumour DNA fraction。這些輸出不等同單一細胞身分、cellular purity 或完整 clone lineage。

### M9

> ONT read 同時提供序列與 methylation observation。固定 exact PS、HP 與 genetic pattern 後，可檢驗區域 methylation 是否仍與分子 pattern 關聯；這是 association-only auxiliary evidence，不能單獨把 reads 指派為不同 cellular clones。

### M12

> Somatic haplotagging 把 reads 指派到局部 germline/somatic haplotype-carrier classes。相同 tag 表示在該模型與 phase block 下共享分子狀態，不代表已確認來自同一癌細胞群。

### M13

> 在特定圖模型與 error assumptions 下，真正候選可能呈現較一致的 local haplotype path；這是可檢驗假說，不是必然規則。若候選本身參與 phasing/haplotagging，HP concentration 不是獨立驗證。

### 全站共通 guardrail

> **Claim ceiling**：read 是分子觀測；local PS×HP state 是局部分子狀態。CN/LOH/purity/multiplicity/CCF 與正交 biological truth 未整合前，不得把 local topology、methylation association 或 HP concentration 稱為 confirmed cellular clone、ancestry 或 whole-sample phylogeny。

## 10. 修正順序與驗收條件

1. **P0：改首頁、M9、M12、M13、SR1**  
   → 驗證：全文搜尋不再用 read/HP/local state 直接推出 cell group/clone；所有 cellular claim 都附獨立 truth 條件。
2. **P0：處理 SR2b**  
   → 驗證：首屏有 obsolete banner/redirect；71,955/85,941 巢狀分母與 0.909 語義不再錯。
3. **P1：修 M6/M7 工具事實**  
   → 驗證：tag vocabulary、`HP:i`/`HP:Z`、`--ont`/`--pb` 與 pinned runtime/code 一致。
4. **P1：把 M9 兩套 denominator 分開**  
   → 驗證：M1=102,842/469,849 與 formal=3/811 各自命名；COLO829=0 附 eligibility 解釋。
5. **P1：加 self-phasing guardrail 與 TO negative evidence**  
   → 驗證：M5–M7、M12–M13 不再把同一候選建立的 phase evidence 當獨立確認；區分 paper benchmark 與 repo generalization test。
6. **P2：為所有精確數字補六欄 provenance**  
   → 驗證：頁面任一 performance/percentage 都可回到 dataset、scope、denominator、metric、version、tier。
7. **重建 print-all 與部署 SR6 決策**  
   → 驗證：25 live routes 無舊敘述；SR6 要嘛部署並標 SPEC，要嘛從 source/navigation contract 明確排除。

## 11. 限制與不應過度解讀之處

- 本報告是 **claim audit**，不是重新跑全套生物流程；數值重算使用已凍結 artifact 與 authority manifest。
- `confirmed cellular subclone=0` 是目前沒有案例跨過正交確認門檻的 claim ceiling，不是腫瘤內真實 subclone prevalence=0。
- 3 個 robust methyl associations 沒有 same-locus cross-sample replication；既不能升格為 clone，也不表示 methylation 完全無資訊。
- HCC1395 與 HCC1395_DORADO 是 technical pair；技術再現性不能冒充跨病人泛化。
- `CURRENT_FOCUS.md` 是 working state；精確 topology/methyl 終值以 20260801 authority bundle 與 machine artifacts 為主。
- 網站往後更新後，這份判定只對 commit `46b6f5b…e605` 有效。

## 12. 執行 IO 與實際輸出片段

### 網站

- **輸入**：`https://ccu-bioinformatics-lab.github.io/lab-tutorial/*.html`；GitHub commit `46b6f5b3016c187ad742fecbfa813f835b09e605`。
- **命令**：對 26 routes 執行 `curl -L -sS -o <temp> -w '%{http_code}' <url>`，並對 25 個 200 response 執行 `sha256sum`，與該 commit 的 `site/*.html` 比對。
- **輸出**：`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/route_status.tsv`。
- **實際片段**：`25 HTTP 200; sr6 HTTP 404; 25/25 live SHA256 match`。

### Authority 與數字

- **輸入**：`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`、`denominator_registry.tsv`、M1 JSON、formal methyl JSON、paired benchmark TSV、TO/technical-pair artifacts。
- **命令**：manifest artifacts 逐筆 `sha256sum`；以 `jq`/TSV numerator-denominator arithmetic 重算 topology、methyl、F1、reproducibility；所有最後採用命令 exit code 0。
- **輸出**：`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/key_number_recheck.tsv`。
- **實際片段**：`13/13 MATCH`; `63506/71955=88.25793899%`; `3/811=0.36991369%`; `7569/14789=51.17993035%`; `43937/54644=80.40590001%`; `ΔF1=0.011159235`。

## 13. Source map

| Claim family | Primary/authority source |
|---|---|
| Website content | [live lab-tutorial](https://ccu-bioinformatics-lab.github.io/lab-tutorial/index.html)；[commit 46b6f5b](https://github.com/CCU-Bioinformatics-Lab/lab-tutorial/commit/46b6f5b3016c187ad742fecbfa813f835b09e605) |
| Exact-PS scope/ceiling/hash | [`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`](/big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json) |
| Funnel denominators | [`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`](/big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv) |
| Formal methyl assay | [`InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md`](/big7_disk/liaoyoyo2001/InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md) |
| Paired-pure read classifier | [`incremental_comparison.tsv`](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1/incremental_comparison.tsv) |
| Current negative/caveat register | [`InterSubMod/docs/CURRENT_FOCUS.md`](/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)；[`InterSubMod/research/autoresearch/evidence_ledger.jsonl`](/big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/evidence_ledger.jsonl) |
| LongPhase-S/TO code contracts | [`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md)；[`/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md`](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md) |

---

**Validated scope footer**：完整檢查 25/25 live pages；另記錄 1 個 source-only/HTTP 404 route。結論適用網站 commit `46b6f5b3016c187ad742fecbfa813f835b09e605` 與 InterSubMod 20260801 authority。這不是 subset、demo 或 cellular truth validation。
