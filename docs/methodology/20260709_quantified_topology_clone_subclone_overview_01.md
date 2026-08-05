<!--
建立: 2026-07-09 | 類型: 7樣本量化總覽(拓撲+clone/subclone) | build: build_quantified_overview.py 一鍵重算
資料: 新骨幹(ClairS PASS) region-view + mlhp col_coverage_by_hp
-->

# 7 樣本分層重建 — 完整量化總覽（拓撲關係 × clone/subclone 頻率）

> 新骨幹（ClairS PASS/LongPhase-S；is_somatic 已移除）。分母=germline(1/2) lineage。一鍵重算：`build_quantified_overview.py`。


## §1 單位階層 + 骨幹 sSNV（每樣本）

| 樣本 | 骨幹 sSNV | region | germline lineage | 形狀確定% |
|---|---|---|---|---|
| HCC1395 | 113,997 | 7,927 | 12,093 | 66% |
| HCC1395_DORADO | 112,387 | 7,958 | 11,568 | 56% |
| COLO829 | 38,196 | 7,774 | 13,699 | 69% |
| H1437 | 75,578 | 8,865 | 14,375 | 74% |
| H2009 | 157,405 | 9,717 | 15,543 | 69% |
| HCC1937 | 49,548 | 2,695 | 3,960 | 84% |
| HCC1954 | 20,969 | 4,023 | 7,635 | 75% |

## §2 拓撲軸（樹形狀 = clone 關係）

> single=單clone · **linear=巢狀subclone後代** · **branched/star=姊妹subclone分歧** · mixed=形狀未定 · capped=太密

| 樣本 | single | linear(巢狀) | branched+star(姊妹) | mixed | capped |
|---|---|---|---|---|---|
| HCC1395 | 38% | 19% | 5% | 33% | 5% |
| HCC1395_DORADO | 43% | 9% | 2% | 42% | 4% |
| COLO829 | 43% | 22% | 0% | 29% | 5% |
| H1437 | 33% | 31% | 2% | 23% | 11% |
| H2009 | 22% | 28% | 3% | 23% | 24% |
| HCC1937 | 39% | 32% | 11% | 16% | 3% |
| HCC1954 | 56% | 12% | 6% | 24% | 2% |

## §3 clone/subclone 頻率軸（VAF；🔴 有 CN confound）

> VAF≥0.75=clonal(主幹) · [0.1,0.75)=subclonal(亞群) · founding=全clonal · has_subclone=≥1亞群突變 · lowcov=覆蓋<6

> 🔴 **CN confound**：CN-gain 稀釋 VAF→誤判 subclone；LOH 抬高→誤判 clonal。**has_subclone% 是 CN-未控上界**；只 HCC1395 有 SEQC2 CN 可控。

| 樣本 | founding主幹 | has_subclone(上界) | weak | lowcov |
|---|---|---|---|---|
| HCC1395 | 10% | 54% | 34% | 2% |
| HCC1395_DORADO | 15% | 52% | 32% | 1% |
| COLO829 | 23% | 11% | 64% | 1% |
| H1437 | 16% | 32% | 52% | 0% |
| H2009 | 12% | 40% | 45% | 3% |
| HCC1937 | 13% | 56% | 31% | 0% |
| HCC1954 | 4% | 58% | 38% | 0% |

## §4 樹選擇（read-likelihood）— 實證結論

> 對 mixed_shape(形狀未定)區用 read 數排形狀：**全 7 樣本僅 5% 可排**，95% 是純隱藏祖先(中間節點 0-read 分不出)。
> → **read 數不適合當 tree-selector**（多數歧義是隱藏祖先）；read 數的正確用途 = §3 clone/subclone 頻率軸。


## §5 HCC1395(5khz) vs DORADO 重現性（同細胞株）

- 共同形狀確定區 3,944；**拓撲一致 80%**、**生物類別(單/巢狀/姊妹)一致 80%**。
- clone/subclone(§3)：has_subclone 54% vs 52%（一致）。→ 不一致主要 depth-driven（5khz 深→linear，DORADO 淺→single）。


## §6 誠實限制

1. 🔴 **clone/subclone VAF 有 CN confound**：has_subclone% 是上界，需 CN 控制（只 HCC1395 有 SEQC2 CN）。
2. **樹選擇 read-lik 僅 5%**：多數樹/形狀歧義是隱藏祖先，read 分不出。
3. **COLO829 65% lowcov**：其低-coread artifact（多 lineage 覆蓋不足）。
4. **region-local scale**：每區送 solver 的位點 ≤8 sSNV；region 由「相鄰 sSNV gap ≤50kb」形成 connected component，**總 span 可 >50kb**。每個 mutation-bearing HP1/HP2 family 分開建局部結構，非全腫瘤整樹。
5. 跨樣本只比比例（絕對受深度混淆）；p 值 pseudoreplication（真 n=7）不報。
