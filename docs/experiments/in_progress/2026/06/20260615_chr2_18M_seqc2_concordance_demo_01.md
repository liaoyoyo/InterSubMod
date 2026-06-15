<!--
建立: 2026-06-15
任務類型: F demo / illustration（scope = chr2:18M 6 位點，非全基因組 → 標 PARTIAL）
data_sources:
  - scripts/seqc2_subclone_concordance.py（比對工具）
  - data/demo_loci_chr2_18M.tsv（輸入 6 位點）
  - data/chr2_18M_seqc2_concordance.tsv / .md（本次輸出真值）
  - /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf（TVAF 來源）
  - /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed（LOH/CN 來源）
  - /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed（可評估遮罩）
provenance_note: 表內每個數字由 scripts/seqc2_subclone_concordance.py 於 2026-06-15 直接讀 SEQC2 磁碟 truth 產出（exit 0），非手填。CCF=TVAF×CN/(purity×mult) 近似值。
-->

# 🟡 DEMO（PARTIAL — 6 位點，非全基因組）— chr2:18M × SEQC2 外部 subclone 比對

> **scope**：HCC1395 chr2:18M 旗艦區 **6 個 subclone-defining SNV**（你方候選樹）vs SEQC2 外部 truth。**非全基因組、非全樣本**；為比對方法的可跑 demo。工具可重用於任意 loci（見文末）。

## 一、做什麼（3 層 concordance，非 1:1 金標準）
- **L1 CN/LOH 一致性**：每位點是否落在 SEQC2 gain/loss/LOH 段。
- **L2 SNV clonal/subclonal**：查 SEQC2 sSNV truth 的共識 `TVAF` → 近似 `CCF = TVAF×CN/(purity×mult)`。
- **HC 可評估性**：是否在 SEQC2 High-Confidence 區（落空隙 = truth-unevaluable）。

## 二、本次驗證結果（grep 真值，purity=0.99, mult=1, LOH 視 CN=2）

| id | 位點 | 你方角色 | SEQC2 狀態 | inHC | CN | TVAF | CCF≈ | 分類 |
|----|------|---------|-----------|:--:|:--:|:--:|:--:|------|
| 1 | chr2:18,068,480 C>G | beta-left | **SEQC2-HighConf** | Y | loh(2) | 0.403 | 0.814 | subclonal-major |
| 2 | chr2:18,072,546 G>C | beta-left | **SEQC2-HighConf** | Y | loh(2) | 0.389 | 0.786 | subclonal-major |
| 3 | chr2:18,086,020 G>A | **alpha-branch（樹的 pivot）** | **out-of-HC（truth-unevaluable）** | **N** | loh(2) | . | . | n/a |
| 4 | chr2:18,096,269 C>G | beta-right（homopolymer）| SEQC2-MedConf | Y | loh(2) | 0.242 | 0.489 | subclonal-major |
| 5 | chr2:18,099,697 G>C | alpha-1（nested, low-VAF）| **SEQC2-HighConf** | Y | loh(2) | **0.048** | **0.097** | **subclonal-minor** |
| 6 | chr2:18,108,828 C>G | beta-right | **SEQC2-HighConf** | Y | loh(2) | 0.382 | 0.772 | subclonal-major |

**region summary**：chr2:18,068,480–18,108,828，全落 SEQC2 **LOH**（`chr2:16,146,119–22,100,000=loh`）；in_truth=**5/6**、evaluable=**5/6**、subclonal(CCF<0.85)=**5/6**。

## 三、解讀（可防守）
1. **L1 LOH 100% 一致** ✅：6 位點全在 SEQC2 確認的 LOH 區 → 你方「此區 LOH」前提外部成立。
2. **L2 出現明確 CCF 梯度** ✅：同一區內 CCF 從 **0.10（18,099,697，minor subclone）→ 0.49 → 0.77–0.81（3 個 major）**。同區位點不是單一 clonal block，而是**有亞群 fraction 差異** → 外部 TVAF 獨立佐證「此區有 subclonal 結構」。其中 18,099,697 的 CCF≈0.10 **吻合 Fang 2021 報告的 subclonal VAF peak（0.15/0.08）**。
3. **🔴 樹的 pivot（18,086,020）外部無真值**：落 SEQC2 HC 空隙 → 無 TVAF/可評估性。你方演化樹最關鍵的 alpha-branch SNV **只能靠內部證據（HKU/DORADO read linkage）**，外部無法確認。
4. **5/6 SNV 是 SEQC2 truth-confirmed somatic**（4 HighConf + 1 MedConf）→ 你方候選 SNV 多數外部真實，非 caller artifact。

## 四、§13 誠實註記（重要）
- 較早一輪手動 `grep TVAF= | head -1` 取到的是**單一 aligner 的 bwaTVAF**（0.396/0.297/0.050/0.389），非共識值。本 script 用 anchored regex `(?:^|;)TVAF=` 取**共識 `TVAF`**（0.389/0.242/0.048/0.382）→ **以 script 值為準**。結論方向（梯度、minor subclone、5/6 confirmed）不變。

## 五、限制（比對時必守）
- **非金標準對拍**：SEQC2 樹是單一 PhyloWGS pipeline（DREAM：19–35% 演算法變異）；本比對是「一致性 sanity check」。
- **CCF 為近似**：mult=1、LOH 視 CN=2；真實 multiplicity/CN 會移動 CCF（±）；CCF≈0.8 可能是「major clone 在誤差內接近 clonal」。
- **跨平台**：SEQC2 TVAF=短讀 bulk；你方=ONT。CN/purity 校正後可比 fraction 帶，非逐 read。
- **passage drift**：SEQC2 stock ≠ 你 HKU/DORADO stock，fraction 不可假設相同。
- **無甲基**：SEQC2 任何檔皆無甲基 → 你方「甲基 coherence」層**無外部檔可比**（只能內部 cross-basecaller）。
- **18,086,020 unevaluable**：樹 pivot 無外部真值（見上）。

## 六、如何重用（全基因組 / 其他樣本）
```bash
# 1) 準備 loci TSV：id  chrom  pos1  ref  alt  our_label
# 2) 跑：
python3 InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/seqc2_subclone_concordance.py \
  --loci <your_loci.tsv> --out-prefix <out_prefix> [--purity 0.99 --mult 1]
# 輸出 <out_prefix>.tsv + .md；預設讀 /big8_disk/data/HCC1395/SEQC2/ 全套 truth
```
> COLO829 可改 `--ssnv-vcf/--cnv-bed` 指向其 truth（NYGC/Craig）後同法跑；其 subclone 另有 Velazquez-Villarreal 2020 單細胞 4-clone 可作 prior。
