<!--
建立時間: 2026-06-15
任務: 獨立審查 20260615_chr2_18M_subclone_verification_HCC1395_01.md 的完整性、可信度、錯誤、矛盾與可防守結論
驗證範圍:
  - 指定 report 01 與其 assets
  - 同日 independent verdict 02 / independent_subclone_audit.py / independent_audit.json
  - 原始 HKU/DORADO tumor-normal BAM、SEQC2 truth/HC/LOH、GRCh38
  - thesis flagship ISM output 與 run.log
  - external_validation 中 subclonal reconstruction 與 ASM 文獻卡
狀態: critical_audit_complete
-->

# chr2:18M report 01 獨立可信度稽核

## 最終裁決

指定報告 `20260615_chr2_18M_subclone_verification_HCC1395_01.md` 是一份**有價值但不可直接當最終真值**的中間報告。

它正確修正了原始五群模型的多個過度宣稱，也抓到真實的 read-level linkage、LOH、局部甲基關聯與 pos4 homopolymer 問題；但仍存在數個會改變主結論的錯誤或未揭露矛盾：

1. **不能稱為乾淨的 somatic-subclone methylation 案例**：matched normal 在 `3.3/3.4/3.5/corrected-4.1` 已有顯著 haplotype ASM。
2. **不能稱樹拓樸被 parsimony「強制」**：beta-left `(1)(2)` 與 beta-right `(4)(6)` 缺直接 genetic bridge，完整 root/tree 未被 spanning reads 觀察。
3. **pos3 不能寫成已證實真 somatic**：它是強 tumor-specific somatic candidate，但位於 SEQC2 HC gap，且未有正交 truth。
4. **DORADO 甲基複製的敘述混合了不同比較軸**：若比較 `pos3 A vs pos3 G`，FDR 顯著的是 `3.1/3.2/3.3`，不是 report 01 所寫的 `2.1/2.2/3.1/3.2/3.3`。
5. **report 01 未直接驗證 thesis flagship pipeline**：旗艦 ISM 實際使用另一個 `longphase-to-mod v6` BAM 與人工建立、註明 `germline_het_anchor_not_somatic` 的 anchor VCF。

因此，最合理總結是：

> **此區域強力支持一組 LOH-constrained、由局部 somatic-allele linkage 定義、並帶 allele-associated methylation structure 的 regional molecular states；它適合作為 ISM 去混淆與 subclonal-structure characterization 案例，但不是已確認的完整 subclone tree，也不是乾淨的 somatic methylation 或演化順序證明。**

## 重大發現

### Critical 1：report 01 與 thesis flagship 不是同一 pipeline

report 01 審核的是：

- HKU LongPhase-S tagged tumor BAM
- HKU matched normal
- DORADO raw/tagged tumor-normal
- 六個候選 SNV 的手動 read-linkage reconstruction

但 thesis 旗艦 ISM 輸出的 `run.log` 顯示實際輸入是：

```text
Tumor BAM: /big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam
Normal BAM: None
Somatic VCF: research/flagship_chr2_18086020_20260612/anchor_18086020.vcf.gz
```

該 anchor VCF 又明確寫：

```text
NOTE=germline_het_anchor_not_somatic
```

同一批 57 個 read name 在旗艦輸出與目前 HKU LongPhase-S BAM 的 HP tag 也不一致：

| Flagship output | Current HKU LongPhase-S BAM | Reads |
|---|---|---:|
| `1-1` | `2-1` | 36 |
| `1` | `2` | 16 |
| `1` | `2-1` | 3 |
| `2` | `1` | 2 |

這表示 thesis 的 `Normal→HP1→HP1-1` 是特定 pipeline 的標籤敘事，不可直接由 report 01 的 `HP2/HP2-1` 結果驗證。HP1/HP2 本來可任意翻轉，但此處不只整體翻轉，另有 5 reads 分類不同，必須先完成 pipeline provenance 對帳。

### Critical 2：「乾淨 somatic-genotype 內甲基分群」不成立

report 01 將案例描述為非 germline ASM 的乾淨 operational-subclonal methylation case；獨立重算顯示 matched normal 已有顯著 ASM：

| CpG | Normal HP1 mean | Normal HP2 mean | MW BH-FDR |
|---|---:|---:|---:|
| 3.3 | 0.814 | 0.046 | 0.0070 |
| 3.4 | 0.033 | 0.850 | 0.0070 |
| 3.5 | 0.583 | 0.012 | 0.0496 |
| corrected 4.1 | 0.887 | 0.078 | 0.0070 |

所以 distal 區段至少部分可由 retained-haplotype / pre-existing ASM 解釋。可防守的說法是：

> 此區同時包含 normal haplotype ASM 與 tumor-associated allele-methylation structure；ISM 的價值是拆解兩者，而不是把全部甲基差異歸因於 subclone。

此外，normal HP1/HP2 與 tumor retained haplotype 的 phase orientation 尚需用共同 germline markers 明確對齊；目前 normal ASM 證明「confound 存在」，但尚未完整量化它解釋 tumor 差異的比例。

### Critical 3：CpG 座標錯誤未在 report 01 揭露

report 01 使用的原始 extraction script 將 `4.1` 寫成 `chr2:18,096,041`；此位置在 GRCh38 是 `A`，不是 CpG。正確圖示位置為 G coordinate `18,096,341`，對應 CpG C coordinate：

```text
chr2:18,096,340
```

原腳本還用 `site / +1 / -1` 搜尋 modification call，適合處理 CpG 兩股座標，但不能修復相差 299 bp 的 `4.1` 誤植。report 01 的原始 4.1 結論因此不可用；應採 independent audit 的 exact-reference normalization 重算結果。

## 高重要性發現

### High 1：pos3 是強候選，不是 truth-confirmed somatic

成立的觀察：

- pos3 `G>A` 位於 SEQC2 HC gap，不能被 SEQC2 判為 TP 或 FP。
- BQ20 tumor：HKU `A/G=29/30`；DORADO `21/27`。
- BQ20 matched normal：HKU `A=0/18`；DORADO `A=0/27`。
- ALT 支持跨 strands、跨 read starts，且兩個資料來源重現。

report 01 的問題：

- 寫成「真 somatic」超過 truth 能力；正確用語是 **strong tumor-specific somatic candidate**。
- 「30/28 分子、0 低品質」不成立：BQ10→BQ20 時 HKU A `30→29`，DORADO A `28→21`。
- 「truth-gap artifact」用語不準確；應寫 **benchmark-unevaluable HC gap**，不是 truth artifact。

### High 2：「0 違反」只在狹義 strict-alpha vs strict-beta 成立

若 beta 嚴格定義為 `(1)G or (2)C or (6)G`，兩資料確實未觀察到 `pos3 A` 同時帶 strict-beta event。

但 report 01 同時把 `pos4 non-ref` 納入 beta-like state；DORADO 的 `(3)A vs (4)non-ref` 有：

```text
00/10/01/11 = 0/7/5/1
```

因此廣義分叉不是「0 違反」。另外，DORADO 對 `(1) vs (3)` 與 `(3) vs (6)` 沒有足夠 spanning reads，不能把「0 可觀測違反」寫成跨資料強複製。

### High 3：beta-left 與 beta-right 未直接連接

直接 read linkage 支持：

- `(1)G + (2)C` 同支。
- `(4)non-ref + (6)G` 在 HKU 同 state。
- `(3)A` 與上述事件主要互斥。
- `(5)C` 巢狀於 `(3)A`，但 direct support 僅 HKU 4 reads、DORADO 2 reads。

未直接支持：

- `(1)(2)` 與 `(4)(6)` 相距約 36-40 kb，缺足夠直接 genetic bridge。

所以不能寫成一個已證實 beta lineage，也不能稱 all-REF root 與 alpha/beta siblings 被 parsimony「強制」。可畫候選樹，但 beta 左右區塊間應用虛線，root 應標 `inferred/unsampled`。

### High 4：DORADO 甲基複製敘述混合不同 null

report 01 的 methylation script 實際把：

```python
pos3 A -> alpha
pos3 G -> beta
```

也就是把所有 pos3 reference-G reads，包括 ancestral/reference reads，都稱為 beta。這與 script docstring 所宣稱的 strict beta 定義不一致。

以 corrected independent audit 的 `pos3 A vs pos3 G` 比較，DORADO 通過 MW BH-FDR 的 CpG 只有：

```text
3.1, 3.2, 3.3
```

`2.1` 只有 `n=1 vs 1`，不可檢定；`2.2 q=0.135`。report 01 把 `2.1/2.2` 寫成 DORADO alpha-vs-beta 顯著，實際上混入了其他 variant-anchor ALT-vs-REF 比較。

跨資料可保留的結論是：

- 多個 **local allele-CpG associations** 可複製。
- `pos3 A vs G` 的直接跨資料穩健核心是 `3.1/3.2/3.3`。
- 不能把所有複製的 allele-CpG association 合併成單一 alpha-vs-beta lineage methylation claim。

### High 5：pos4 應撤回三 subclone，但不能說 G/T/DEL 全是同一 artifact

SEQC2 truth 明確支持 pos4 `C>G` MedConf；因此真實 biological `C>G` event 存在。poly-T context 與兩資料中的 G/T/DEL 混合則表示 read-level allele assignment 容易受 homopolymer/basecalling/alignment jitter 影響。

正確結論：

> `G/T/DEL` 不足以定義三個 biological subclones；應視為 truth-supported `C>G` 位點周圍的 homopolymer-uncertain observed state。

不應寫：

> `G/T/DEL` 三者全部都是 artifact。

## 中重要性與統計限制

1. **normal=REF 必須帶過濾條件**：raw normal 中仍有 DEL/T/C 等非 REF call；正確說法是「unique primary、MAPQ/BQ 過濾後，未觀察到預期 somatic ALT」。
2. **all-REF read 數字應統一**：independent audit 的 ≥2 點 all-REF 為 HKU `17`、DORADO `15`，不是 report 01 的 `22/17`。
3. **HKU normal dedup 敏感性未完成**：390 records→195 unique names，且 16 個 duplicate names 的 records 不一致；normal ASM 結論應補 dedup-selection sensitivity。
4. **候選位點與 CpG 是事後選擇**：目前 BH-FDR 只在每次 10 CpG 比較內校正，未涵蓋多 variant anchors、多 grouping 定義與手動挑選 locus 的 selection burden。p/q 值應視為 exploratory。
5. **「0 methyl-lineage contradiction read」缺可重現定義**：甲基是連續機率且多 reads 有混合 pattern；沒有預先定義的 per-read contradiction rule，不能宣稱 0。
6. **HKU 與 DORADO 不是單純同 reads 換 basecaller**：此區 read name overlap=0，BAM header 也顯示不同資料來源。它是有價值的 cross-dataset/technical replication，但仍不是獨立 biological replicate。
7. **Flagship ISM run 本身 `Total Significant=0`**：雖 cluster PERMANOVA `F=10.612, p=0.01`，其最終 significance gate 顯示 `Total Significant: 0`；論文必須區分 cluster structure、HP association、allele association 與 final significant flag。

## Claim-by-claim 修正版

| Report 01 claim | 稽核裁決 | 可防守修正版 |
|---|---|---|
| 5 TP + 1 FP | 部分錯誤 | 5 個 SEQC2 truth-confirmed + 1 個 out-of-HC 強 tumor-specific candidate |
| 分叉乾淨、0 違反 | 過寬 | strict alpha-vs-beta 未見違反；廣義 pos4-beta 有 1 DORADO discordant，且多遠距 pair 無 coverage |
| 約 3 lineage | 過度確定 | 至少 3 個 distinguishable regional states；beta-left/right linkage 未證 |
| 甲基把兩 lineage 分開 | 部分成立 | 多個 local allele-CpG association 跨資料重現；單一 alpha-vs-beta methylation claim 不成立 |
| clean subclonal methylation | 不成立 | mixed normal-ASM + tumor-associated allele methylation case |
| LOH | 成立 | SEQC2 LOH BED + HP imbalance 支持 |
| normal 六點全 REF | 需改寫 | 過濾後未見預期 ALT；raw calls 並非全 REF |
| pos3 真 somatic | 過度宣稱 | strong tumor-specific somatic candidate；需正交驗證 |
| all-REF root→alpha/beta siblings 被強制 | 不成立 | 可防守的 parsimonious candidate topology，非唯一或已確認 tree |
| pos3→pos5 | 有界支持 | perfect nesting among observed spanning reads；僅 4+2 direct reads，provisional |
| pos4 G/T/DEL 三 subclone | 否定 | truth-supported C>G locus with homopolymer-uncertain read calls |

## 與 ISM 研究目標的邏輯結論

此案例**支持** ISM 的核心價值：

1. 在長讀 read 層次聯合讀取 mutation、haplotag 與 native methylation。
2. 找出跨資料重現的 allele-CpG association。
3. 用 matched normal 揭露 retained-haplotype ASM confound。
4. 在單樣本中提出有明確不確定性的 regional molecular-state / subclone hypothesis。
5. 顯示單一 fine HP tag 不能完整代表 subclone，甲基可提供附加 characterization。

此案例**不支持**：

1. 甲基造成突變。
2. 甲基獨立重建 subclone。
3. 完整 clone number、clone fraction 或 phylogeny。
4. 已確認完整演化順序。
5. 跨病人或跨 biological replicate 泛化。

這與外部方法學邊界一致：

- [Tarabichi et al., Nature Methods 2021](https://doi.org/10.1038/s41592-020-01013-2)：正式 subclonal reconstruction 需要 genetic/CN/purity/clustering/phylogeny 與不確定性處理。
- [Onuchic et al., Science 2018](https://doi.org/10.1126/science.aar3146)：sequence/haplotype-dependent ASM 廣泛存在，支持 matched-normal baseline 是必要 confound control。
- [Do et al., Genome Biology 2020](https://doi.org/10.1186/s13059-020-02059-3)：cancer ASM 真實存在，但其跨個體/DMR 設計不能直接證明本案例的 somatic-subclone specificity。

## 升級為論文旗艦前必做

1. **凍結 canonical pipeline**：明確選定 LongPhase-S、longphase-to-mod v6 或其他版本，對帳 57 reads 的 HP tag 差異，禁止混寫 `HP1/HP1-1` 與 `HP2/HP2-1`。
2. **重跑同一 canonical BAM 的完整證據鏈**：ISM cluster、六 SNV linkage、corrected exact CpG、normal-anchor、DORADO replication 必須指向同一分析定義。
3. **正交驗證 pos3**：使用獨立 short-read/PacBio/duplex 或 targeted assay；完成前保持 strong candidate。
4. **normal retained-haplotype 對齊**：用共同 germline markers 對齊 normal/tumor phase orientation，再做 within-retained-haplotype somatic-vs-baseline null。
5. **pos4 局部重比對/原始訊號驗證**：區分 truth `C>G` 與 T/DEL jitter，禁止以 G/T/DEL 定義三 clone。
6. **預先註冊統計 family**：固定 grouping、CpG 範圍、FDR family、最低 n、discordance rule，並在第二 LOH locus 或第二樣本複製。
7. **限制樹主張**：beta-left/right 未有 direct genetic bridge 前，候選樹用虛線並標 uncertainty；不使用 confirmed reconstruction。

## 可重現性確認

- `independent_subclone_audit.py` 已重新執行成功。
- 重跑 JSON 與既有 JSON SHA256 完全一致。
- 重跑 Markdown 與既有 Markdown SHA256 完全一致。
- 腳本已通過 `python3 -m py_compile`。

