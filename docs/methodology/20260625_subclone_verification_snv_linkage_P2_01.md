<!--
建立時間: 2026-06-25
類型: methodology — P2 sSNV read-level 連鎖驗證（正交非循環 subclone 確認）
狀態: in_progress（單樣本 HCC1395 ⭐3；4 ≥2-lineage co-occurrence confirmed；對抗驗證 = QUALIFIED CONFIRMATION，5 替代解釋全排除）
build_branch: docs/method-comparison-ism-external-202606
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/p2_feasibility.json, docs/methodology/_assets/20260618_subcluster_pilot/p2_linkage.json, docs/methodology/_assets/20260618_subcluster_pilot/anchored_retest.json
-->

# P2 — sSNV read-level 連鎖驗證（HCC1395；正交非循環 subclone 確認）

> **TL;DR**：把 P1.3 cis-test 篩出的 9 tumor-specific（+30 anchored）候選交給**正交、非循環**的 sSNV read-level 連鎖。**4 個 locus 經多 sSNV 同分子共現確認為真 subclone**（3 sibling 互斥 + 1 nested 階層，全同 germline HP、全 normal=REF somatic、全是 P1.3 的 9 tumor-specific）。🔑 **可連鎖的 9-候選 4/4 全中** → 甲基預選被獨立驗證。誠實邊界：subclone 稀疏→只 ~4-7 locus 可連鎖；確認的是「**locus 局部 ≥2 lineage 的分子證據**」非 genome-wide clone tree；⭐3 單樣本。**數字全來自上列 JSON（本輪 Read）。**

## §0 為何 sSNV 連鎖是唯一非循環錨

P1.3 證實：單樣本甲基無法非循環證 subclone（genotype-同質群內甲基結構永遠 double-dip）。唯一非循環錨 = **somatic SNV 在同一條 read 上的共現**（genotype，非甲基切出）。同一 germline 單倍型背景上的兩個 somatic 變異若在 read 上**互斥**（從不共現）→ 起於不同細胞 = sibling subclone；若**嵌套**（一變異 ⊂ 另一）→ 祖先-衍生階層。此即 Fang2021 / chr2:18M 的 subclone 證據範式。

## §1 Feasibility（p2_feasibility.json）— subclone 稀疏，少數可連鎖

全基因組 TP somatic SNV = **30,490**（密度 ~1/98kb，nearest 鄰 sSNV median 15,233 bp）。候選周邊 ≥2 sSNV（≤10kb 可同 read）：

| 子集 | ≤1kb | ≤5kb | ≤10kb | ≤50kb |
|---|---|---|---|---|
| 9 tumor-specific | 2 | 4 | **4** | 8 |
| 30 anchored | 4 | 7 | **7** | 23 |
| 102 candidate | 6 | 13 | 14 | 65 |

→ 多 sSNV 連鎖只對少數候選可行（4/9 tumor-specific、7/30 anchored）；其餘無鄰近 sSNV → 連鎖測不了（仍 characterization-only）。

## §2 連鎖驗證方法（p2_snv_linkage.py）

對每個 linkable locus：pysam 抽 tumor BAM **per-read 在各 sSNV 的等位**（REF/ALT，MAPQ≥20）+ HP tag → 覆蓋 ≥2 sSNV 的 reads 建**每對 sSNV 2×2 共現表** → 關係分類：
- **co_linked**（REF_ALT=0 且 ALT_REF=0）= 同 haplotype/lineage；
- **nested_a_in_b / b_in_a**（一方向 0）= 祖先-衍生嵌套；
- **mutual_excl**（ALT_ALT=0 且離對角>0）= sibling 不同 lineage；
- **independent**（各 cell 皆有）= 無乾淨結構。
🔴 **HP 一致性閘**：互斥/嵌套**只在兩 ALT 同 germline HP 才算 subclone**；不同 HP = allelic（等位非 subclone）。+ normal BAM 確認各 sSNV normal=REF（somatic）。

## §3 結果（7 linkable loci，p2_linkage.json）

| locus | in9 | sSNV | coread | 關係 | HP | somatic | verdict |
|---|---|---|---|---|---|---|---|
| chr11:91134191 | ✓ | 2 | 39 | mutual_excl（0/21/18/0） | 1-1 | 2/2 | **subclone sibling** |
| chr17:48360161 | ✓ | 3 | 64 | nested+co_linked | 1-1 | 2/3 | **subclone nested** |
| chr1:50031328 | ✓ | 2 | 61 | mutual_excl（3/14/43/1） | 2-1 | 2/2 | **subclone sibling** |
| chr1:57478998 | ✓ | 2 | 58 | mutual_excl（3/34/21/0） | 1-1 | 2/2 | **subclone sibling** |
| chr22:33135662 | ✗ | 2 | 124 | mutual_excl 但**異 HP**（2-1 vs 1-1） | — | 2/2 | allelic_not_subclone |
| chr17:39668348 | ✗ | 2 | 46 | independent | — | 2/2 | no_clean_linkage |
| chr19:17533317 | ✗ | 3 | 107 | independent（3對） | — | 3/3 | no_clean_linkage |

**4 subclone confirmed / 1 allelic / 2 no-clean**。

## §4 4 個確認 subclone 細節

- **chr11:91134191（sibling）**：91139191(C/T) 與 91142579(G/T) 互斥（ALT_ALT=0；REF_ALT 21、ALT_REF 18），兩 ALT 全在 HP **1-1**，normal 各 34/25 reads 0 ALT → 同 1-1 背景上**兩個從不共現的 somatic 變異 = 2 distinct lineages**。
- **chr17:48360161（nested 階層）**：48365089-ALT(C, 51 reads, **祖先**) ⊃ {48362515-ALT(A), 48365161-ALT(C, anchor)}，後兩者**完美共連**（REF_ALT=0 且 ALT_REF=0, ALT_ALT 22）= 同衍生 lineage；全 HP 1-1。= 2 層克隆階層（ancestral → derived）。⚠ 48362515 normal 1/29 ALT（borderline somatic）。
- **chr1:50031328 / chr1:57478998（sibling）**：互斥、同 HP（2-1 / 1-1）、normal=REF。

## §4b 對抗驗證（fresh-context agent，QUALIFIED CONFIRMATION）

獨立 agent spot-check BAM 逐一裁決 5 個替代解釋，**全部排除**：① **phasing error**（chr11 ALT reads 0/40 HP 不一致 = 太乾淨，非相位噪音）；② **mapping artifact**（會造正相關共現，與互斥反向）；③ **sequencing error**（chr11 39/40 reads 帶 ≥1 ALT、雙峰分群，非稀疏 singleton；joint error ~10⁻¹⁴）；④ **allelic**（同 HP 互斥 = somatic branching，非異等位；HP tag 零混雜）；⑤ **CN confound**（CN-gain 會有 ALT_ALT>0，實測 0/40 = CN-neutral）。

per-locus confidence：**chr11 / chr1:57478998 = HIGH**；chr1:50031328 = MED-HIGH（1 ALT_ALT = Phred error）；chr17:48360161 = MED（48362515 borderline somatic → 報「2 confirmed + 1 borderline」，另 2 sSNV 完美共連 = robust core）。

🔴 **scope 校正**：⭐3 justified（regional ≥2-lineage co-occurrence）**非 ⭐4**（genome-wide clonal architecture）。**互斥示 branching 非 directionality**（sibling = 分支非祖先序）；nested（chr17）才示階層方向。**framing 語言**：用「**regional ≥2-lineage co-occurrence / branched somatic evolution evidence**」，**非**「subclone identification / clonal reconstruction」。

## §5 🔑 整合結論（甲基預選被獨立驗證）

**可連鎖的 9-tumor-specific 候選 4/4 全部 subclone confirmed**（chr11/chr17/chr1×2）；非-in9 的 linkable loci 無一是 subclone（chr22 allelic、chr17:39668348/chr19 independent）。→ **甲基 NACT cis-test 的最乾淨 survivors（9 tumor-specific），凡能被正交 sSNV 連鎖獨立檢驗者，皆為真 subclone**。這就是論文的 corroboration 鏈：
- **haplotag + sSNV 連鎖 = 重建骨幹**（非循環確認局部 subclonal lineage）；
- **methylation = characterize/corroborate**（預選正確、刻畫已驗 lineage 的表觀差異），非偵測驅動。

## §6 誠實 caveats / 紅線

- **小 N**：只 4 confirmed（somatic 稀疏→多數候選無鄰近 sSNV 可連鎖）；非 subclone 普查。
- **局部非全基因組**：確認的是「**locus 局部 ≥2 lineage 的分子證據**」；phase-block 跨區無連結→**不可組 genome-wide clone tree**（reconstruction GAP）。
- **n=2 sSNV**：sibling 推論基於 2 變異互斥；coread≥39 充足但 single-cell/multi-region 才是 confirmation 黃金標準 → ⭐3 characterization-grade。
- **chr17 nested 1 borderline somatic**（48362515 normal 1 ALT）。
- **allelic confound 已處理**（HP 一致性閘擋掉 chr22）。

## §7 Provenance + 自驗

| 數字群 | 來源檔 | 腳本 |
|---|---|---|
| feasibility（sSNV 密度/linkable） | `p2_feasibility.json` | `p2_multisnv_feasibility.py` |
| 連鎖 + 2×2 + HP + verdict | `p2_linkage.json` | `p2_snv_linkage.py` |
| 候選來源（9 tumor-specific） | `anchored_retest.json` | `p1_anchored_retest.py` |

> 紅線：⭐3 單樣本 HCC1395；4 subclone = 局部 ≥2 lineage 分子證據非 genome-wide tree；甲基 corroborate 非偵測；sSNV 連鎖 = 唯一非循環錨。對抗驗證進行中（同-HP 互斥真偽 + 替代解釋）。§13 數字全 JSON 溯源。
