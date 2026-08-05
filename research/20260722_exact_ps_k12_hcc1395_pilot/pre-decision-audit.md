<!--
建立時間: 2026-07-22
目標: 在 HCC1395 上驗證 exact PS × HP fail-closed 分區與 k<=12 read-supported segmentation
處理範圍: PARTIAL / exploratory pilot / HCC1395 chr1-22 only
關聯檔案: InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/implementation-notes.md
-->

# HCC1395 exact PS × HP、k≤12 實作前稽核

> **PARTIAL / EXPLORATORY PILOT**：本輪只處理 HCC1395 chr1–22；不可作七樣本、production 或論文終版證據。

## 1. Verdict

**PROBE**。先以 HCC1395 驗證資料契約、Python/C++ parity 與舊版差異；通過後才由使用者決定是否推廣至其餘樣本。

服務目標：G1、G2、G4、G5；LongLineage 對應 LL-G1、LL-G2、LL-G3、LL-G4、LL-G5。

## 2. 問題與決策

舊描述性流程以 HP family 聚合不同 PS block，可能把沒有全域方向關係的 local HP1／HP2 當成同一條 haplotype。HCC1395 舊資料中 865 個 multi-PS 建樹區域，有 802 個（92.7168%）曾出現 HP family 跨 PS 聚合，因此不能默認安全。

本輪固定契約：

1. primary evidence unit 是 `(dataset, chromosome, HP family, exact non-missing PS, read-linked component)`；
2. 每條 molecule 只屬於 sidecar 指定的單一 exact PS × HP unit；
3. 不同 PS 的 constraint 永不默認合併，沒有 signed SAME／FLIP anchor 就維持分開；
4. HP3、missing HP、missing PS 不進入 primary topology constraint，另列 ABSTAIN／排除統計；
5. k＞12 只在同一 exact PS × HP component 內，以 read-support objective 分成連續、非重疊的 k≤12 blocks；
6. 分割後保留與切斷的 constraint 分開守恆；不宣稱小塊可無損拼回原 component 的 global tree；
7. max k=12 是目前 exact solver 的工程上限，不是生物學上的 clone 上限。

## 3. 關鍵假設與反證條件

| 假設 | 風險 | 反證／處置 |
|---|---|---|
| authoritative sidecar 的 HP、PS 可與 BAM alignment exact join | identity mismatch 會錯配 read | 任一 exposed alignment 缺 sidecar row或 identity 衝突即 fail closed |
| scalar PS 足以隔離 local phase block | 它不能推導跨 block SAME／FLIP | 本輪禁止跨 PS 合併；未建立 signed bridge claim |
| fixed R/A calls 可定義 linkage hyperedge | partial X 不能當 fixed link | X 保留給後續 likelihood，但 segmentation hyperedge 只使用 fixed R/A positions |
| 每個 exact unit 的連續分割可作局部運算 | 跨切點 evidence 會遺失 | 每個 constraint 僅一次分類 retained/cut/unavoidable，並量化 retention |
| k≤12 block 可交給 exact solver | read support 仍可能不足 | 無足夠支持者必須 ABSTAIN，不輸出 winner |

Hard NO-GO 條件：跨 PS block 出現在任何 primary block；unit-site 未完整且唯一指派；constraint 質量不守恆；Python/C++ semantic parity 不一致；輸入不是已鎖定的 big7 HCC1395 BAM。

## 4. Step → Verify

1. 鎖定 big7 HCC1395 BAM、BAI、VCF、HP/PS sidecar。  
   → 驗證：輸入絕對路徑、size、抽樣 SHA-256、index quickcheck 與 frozen authority 相符。
2. 建立 exact PS × HP evidence units。  
   → 驗證：每個 block 恰有一個 HP、一個 non-missing PS；`cross_ps_blocks=0`、`cross_hp_blocks=0`。
3. 同 unit 內對 k＞12 執行 read-supported DP。  
   → 驗證：所有 unit-site memberships 恰好指派一次、`max(block_k)<=12`、constraint total=retained+cut+unavoidable。
4. 建立獨立 C++ parity kernel。  
   → 驗證：synthetic negative fixtures 先通過；Python/C++ block membership、disposition、summary semantic SHA 完全一致。
5. 先跑小型 chromosome smoke，再跑 HCC1395 chr1–22。  
   → 驗證：完整命令退出碼 0、產出 receipt、舊版／新版比較表、PARTIAL claim boundary。

## 5. 範圍邊界

本輪回答「是否能避免亂跨 PS、如何合理分解 k＞12、Python/C++ 是否一致、HCC1395 數據改變多少」。本輪不回答跨 PS orientation、不完成七樣本、也不把 local mutation-state tree 直接稱為真實 clone／subclone lineage。

## 6. 2026-07-23 production adapter gate

**Verdict: GO（HCC1395-only）**。使用者已確認 exact PS 應為 fail-closed evidence boundary；本次只把已通過 22/22 chromosome 的 extraction、Python partition 與 C++ parity 接到 topology adapter，不啟動七樣本。

### 研究問題與假說

- 問題：舊 production 把 `PS` 只當 QC，同 HP family 的不同 PS read 會共用一棵樹；能否以 exact-PS blocks 取代該輸入，並完成 HCC1395 topology 重算？
- 假說：若 tree input 單位固定為 `exact PS × HP × read-linked component × bounded block`，則 primary tree 不再出現 cross-PS constraint，且 non-capped units 可通過 V1–V7。
- 成功：`cross_ps=0`、`cross_hp=0`、Python/C++ mismatch=0、stable region key 無覆寫、constraint disposition/molecule weight 守恆、所有 non-capped tree units V1–V7 PASS。
- 失敗：任一 cross-PS/collision/parity 違規；或 adapter 無法把 retained sparse pattern 無重複地投影到 block。

### 權威邊界

1. Production C++ `inter_sub_mod` 目前是單 anchor-sSNV 甲基化 region analyzer，沒有 multi-sSNV topology schema；不在 dirty `RegionProcessor` 中假裝加一個 PS 欄位就完成建樹。
2. 新路徑將 standalone exact-PS C++ kernel 納入可重現 build/validation gate，Python adapter 才產生 layered-tree schema。
3. `k≤12` 是 segmentation bound；若 exact enumeration 觸及 budget 仍須標 capped，不可報為 complete topology。
4. Read weight 本輪用於 pattern MINREAD 與 segmentation 保留量；現有 Steiner solver 仍以觀測 pattern presence 枚舉，尚未是 read-likelihood/VAF ranking。

### 修改前 SHA baseline

```text
run_layered_v3.py                 3f036b653c4580e1589f6d4ca60f2862f29575a2a0d3a69c6a6de52bd4583d5f
sm_linkage_genomewide.py          83ef74971566459dc7bfeeeb7842f2a6cf0056d033171d8c49326c862d574677
sm_multilocus_combinations.py     12ae8e9b79fecc7e66266cf39e13b10a81723dd954d1ff98c3ad22e0434e10bc
tree_enumeration_solver.py        36727f4e1d8d7ce8abf869606211c93d8c1a0506dd7d142e855863c412ca0d61
layered_tree_reconstruction.py    70d28fa4dfe3e69a0ea610e346457a894e7e024ffdb437895343fc66516639d0
build_region_view.py              cd39bb7799f6b62190f9dca2a3afad321d0b384e254ab3fb4a102e594d02d872
exact_ps_k12_partition.py         b71215a8f52cbf5eac0755f93cc28e4ee8490919cfb5c276f1e989dba450e39d
exact_ps_k12_partition.cpp        5ce6e8170dd97c5e10b2df4487400fbef227ac6144dcf5dc5e9de38970829dcf
```

### Step → Verify

1. 建立 exact-partition → mlhp adapter。  
   → 驗證：every retained constraint 只投影到一個 same-PS/same-HP block；stable `region_id` 無重複。
2. Layered driver 傳遞 PS/component/block identity。  
   → 驗證：trace/JSON 可從 tree 回溯 exact PS 與 source block。
3. C++ kernel 成為正式 build target，並作為 adapter precondition。  
   → 驗證：CMake build exit 0；cross-PS negative fixture 被拒絕；Python/C++ semantic parity=0 mismatch。
4. chr22 smoke 後再跑 HCC1395 chr1–22。  
   → 驗證：舊/新同分母比較表、machine-readable receipt、PARTIAL claim boundary。
