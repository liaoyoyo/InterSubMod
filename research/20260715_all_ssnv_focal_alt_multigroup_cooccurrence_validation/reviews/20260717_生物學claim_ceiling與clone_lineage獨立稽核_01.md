<!--
建立時間: 2026-07-17
目標: 獨立檢核 focal-ALT 甲基多群、sSNV 共現與 clone/lineage claim ceiling
處理範圍: claim-contract-v5、final dataset/report builders、active execution state
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# 生物學 claim ceiling 與 clone/lineage 獨立稽核

## 初審判定

**NO-GO for cellular clone/lineage claim；GO only for bounded molecular claims after report wording fixes.**

bulk long reads 每次只觀察一條分子，不能把 HP1 與 HP2 分子配回同一細胞。因此即使資料符合
RR、AR、AA 且 RA 接近固定錯誤上限，也只能說與 focal-before-partner 的簡化模型相容；不能確認
clone count=2、兩個 cellular subclones、直系祖先或唯一 linear tree。

## Claim ladder 判讀

| 層級 | 可支持 | 不可支持 |
|---|---|---|
| M1 | focal-ALT read-level epigenetic heterogeneity | subclone prevalence、clone count |
| M2 | 排除八個已量測軸後的 residual robust epigenetic partition | latent cellular state 已被識別 |
| G1 | same-read partner allele 的 local molecular co-segregation candidate | cellular clone |
| G2/R1 | multi-marker molecular-haplotype candidate 與 partition robustness | cellular identity 或 lineage |
| B1 | 與預設 fixed-error relation model 相容或不相容 | 唯一 mutation order、一般 branching 已排除 |
| Tumor-REF/normal | 削弱 locus-wide/background/normal 替代解 | clonality 證據 |
| C1 | 需 prespecified focal-partner joint CN/CCF model | 現行 focal-only annotation 不足以 PASS |
| L1/L2 | 需 single-cell、colony、spatial、multi-region/longitudinal 證據 | bulk marginal 直接升級 |

## 初審 findings

1. G1 family 文字誤把 999-permutation 可執行性寫成 global family membership 條件；實作是先納入
   全部 M2 exact-testable pairs 做 global BY，再套 effect、conditional permutation 與 callability。
2. G2=0 fallback 未區分 G1 是否也為 0，可能在 G1=0 時仍宣稱局部共分離。
3. C1=0 應寫成缺 joint model 的結構性 `NOT_EVALUABLE`，不能寫成 CN/CCF 否證 0 個候選。
4. 「共同 ancestral ALT 被觀察」應改成「資料與共同 ancestral ALT 模型的分子預測相容」。
5. background superset 只成立於同一 background payload 的 lenient/ARI-qualified predicates，不能拿
   ALT 與 REF 的實際 flag 集合互比。

## 必須保留的結論

使用者提出的 clone1/clone2 情境是合理且可產生可檢驗預測的生成模型。甲基群若與同一 read 上的
partner R/A 共分離，可把證據升級到 latent molecular-state 或 molecular-haplotype candidate；它仍是
提出候選與輔助分離，而不是 cellular subclone 的確認。L1 需要正交 cellular identity；L2 另需可識別
mutation order。最終數值必須等 active producer 與 source-attested release 完成後再審。

## Reviewer

- Independent agent: `019f6e3d-1027-7cc1-a755-b06922ff2e04`
- 查核方式: 唯讀 source、contract、CURRENT_FOCUS 與 active process inspection
- 初審範圍不包含尚未產生的 final G1/G2/R1/B1/C1 數值

## Report-layer closure

上述 5 項文字／fallback findings 已在 report builder 與回歸測試封閉：G1 family 先全域 BY；G1=0
不再宣稱局部共分離；C1=0 明示結構性 `NOT_EVALUABLE`；ancestral ALT 與 background superset 用語均限定
在可觀察模型預測與同一 payload predicates。Python compile exit=0，定向測試 `47 passed`。此 closure
只涵蓋敘述與報告邏輯；最終數值仍須在 source-attested release 後另做 result-level review。
