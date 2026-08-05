<!--
建立時間: 2026-06-21
報告類型: 任務規格 + 背景知識（handoff，給其他 AI / session 接手 V4 單分子 multi-sSNV co-occurrence）
任務類型: D handoff — 任務拆解 + 背景 + 暴力先算→後驗證流程
build_branch: research/subclonal-reconstruction-202606
status: task-spec（尚未實作；chr2:18M 旗艦為手動先例，本任務=系統化它）
data_sources:
  - memory project_chr2_18m_subclone_locus_verification（手動先例：6 sSNV，α/β 互斥 0 違反）
  - external_validation/axis2 foltz-somatichaplotype-linkedread-2024（genetic-only 先驗，PMID 39149342）
  - InterSubMod/docs/method_comparison/20260621_subclone_confirmation_construction_and_ONT_transferability_01.md（V4 定位、移植裁決）
relurl: 紅線同主軸（reconstruction 由 somatic haplotagging 驅動、甲基佐證、regional proof-of-concept）
-->

# 任務規格 + 背景：單分子 multi-sSNV co-occurrence（V4）

> **給接手的 AI**：這是把 ISM **chr2:18M 旗艦的手動分析**系統化成「**暴力先算 → 後驗證**」流程的任務。先讀 §1 背景理解原理，再看 §2 怎麼做、§3 困難、§4 兩階段流程、§6 紅線。**它是 ISM 在單樣本上最有「直接證據」味道的一塊，但只 LOCAL（phase-block 內）。**

## §0 一句話 + 目標
**一句話**：一條 long read 同時跨越多個 somatic SNV 時，「哪些 read 帶哪一**組合**的 ALT」就是**局部譜系的直接證據** —— 不需要 VAF→CCF 推論，直接從同一條分子上讀出「這幾個突變是否在同一個 lineage」。
**目標**：把 chr2:18M（手動做過 6 sSNV）系統化成一個工具 + 驗證協議，能對任意 region 暴力枚舉 co-occurrence、分類 lineage 關係、再過驗證。
**定位**：proof-of-concept / regional characterization，**非** genome-wide tree（紅線見 §6）。

---

## §1 背景知識（原理 — 先理解這個）

### 核心：read 跨 2 個 somatic SNV → 4 種組合 → lineage 關係
一條 read 跨越 somatic SNV **A** 與 **B**，在這條分子上只有 4 種可能：
`A+B（兩個都 ALT）` / `A-only` / `B-only` / `neither（REF-REF）`。

**觀察到哪些組合 → 推哪種 lineage 關係**：

| 觀察到的組合 | lineage 關係 | 意義 | 證據強度 🔴 |
|---|---|---|---|
| 只見 {**A-only, B-only**}，無 A+B | **互斥 (mutually exclusive)** | A、B 在**不同 subclone**（不同分支）| ✅ **可靠**（標準推論，Tarabichi：同 copy 互斥 SNV→branching）|
| 只見 {neither, A-only, **A+B**}，無 B-only | **巢狀 (nested) 候選** | 與「B 在 A-clone 內後發生」**一致** | ⚠ **必要非充分** —— 同 haplotype/molecule co-occurrence **不能單獨推同 subclone**，須加 VAF/fraction（見下 🔴）|
| 4 種組合都見 | **獨立 / germline / artifact** | 需 normal 排 germline；或真獨立事件；或定序錯誤 | 須 normal + QV |

> 🔴 **2026-06-21 外部驗證更正（Foltz NRAS 反例，PMID 39149342）**：NRAS **G13R + Q61K phase 到同一條 haplotype（H2）**，但 **VAF 不同（35.7% vs 22.2%，Q61K relapse→0）→ 它們是不同 subclone 的獨立事件**。→ **「co-occurrence/同-haplotype → 同 lineage」是 *necessary-not-sufficient***；**互斥→不同 subclone 可靠，但 共現→同譜系必須再加 cellular fraction（VAF/CCF）才能定 lineage**。先前「不需任何模型推論的直接讀出」**錯** —— 至少需要 fraction。⚠ Foltz 細節為 search-snippet（full-text 403），投稿前過 /citation-verification 核 PMC11326269。

### 這是什麼的轉移
- = **scSNV-tree 的 co-occurrence primitive**（SCITE / COMPASS：哪些突變在同一細胞）**轉移到 read/block span** → 但**只 LOCAL**（read 只跨 ~kb-Mb，兩 SNV 必須同 read 或同 phase block）。
- **先例 1（genetic-only）**：**Foltz et al. 2024 "SomaticHaplotype"**（PMID 39149342，linked-read 骨髓瘤）：NRAS **Q61K vs G12 無 ALT-ALT barcode → 推獨立 subclone、Q61K 後消失**。已 formalize 此邏輯（但 linked-read、零甲基、無結構檢定）。
- **先例 2（ISM 自己，手動）**：**chr2:18M 旗艦**（HCC1395，6 sSNV）：α/β **互斥 0 違反**、LOH HP1 1.4% → 手動證了局部 lineage。**本任務 = 把這個手動分析變成系統化工具。**

---

## §2 如何實現（暴力先算 — 可先做的部分）

**輸入**：tumor BAM（已 somatic-haplotag）+ somatic VCF（ClairS-TO 等）+ normal BAM + LOH/CN BED。

**步驟（brute-force，無推論）**：
1. **找可跨對/集**：列出彼此距離 < read 長度（或同一 phase block 內、經 haplotag 可連結）的 somatic SNV 對 / 集。
2. **逐 read 取 base**：對每對，掃過該 region 的 read，取兩位點的 base call（ALT/REF/其他）；過 **base-quality + mapping-quality** 濾（ONT error 防呆，見 §3）。
3. **建 co-occurrence 計數表**：2×2（A×B）或 multi-way（≥3 SNV 的組合計數）。同時記每組合的 read 數、HP-tag、strand。
4. **分類 + 統計**：依 §1 表分 nested / mutually-exclusive / independent / germline；配統計量（Fisher exact on 2×2，或 LD-style D′/r²；⚠ 計數小須慎，見 §3）。
5. **建局部 micro-lineage**：從「哪些組合存在」推 phase-block 內的 lineage 拓撲（不跨 block）。
6. **落檔**：每對 → TSV/JSON（位點、組合計數、分類、統計、HP/strand/coverage），供 §4 Phase 2 驗證。

> **「暴力」可行性回答用戶**：✅ **Phase 1 可以暴力先算** —— 它是 cheap、exhaustive、**無推論**的枚舉（只是數同一 read 上的 base 組合），不需要任何模型。難的不是算，是**解讀**（§3）與**驗證**（§4）。

---

## §3 困難與限制（必懂，否則會 over-interpret）

| 困難 | 說明 | 緩解 |
|---|---|---|
| **read span 瓶頸（最大）** | 多數 somatic SNV 對**不在同一 read 上**（somatic 稀疏、相隔常 >read 長）→ 跨兩位點的 read 覆蓋是限制因子 | 只對「可跨對」做；用 haplotag 延伸（但有代價↓）|
| **haplotag 延伸的代價** | 可用 phasing 把連結延伸到整個 phase block，但就從「**同分子**」變「**同 haplotype**」→ 弱化，依賴 phasing 正確 + **跨 phase-set 無連結**（物理天花板）| 明標「same-molecule」vs「same-haplotype」兩級證據強度 |
| **低 VAF / carrier reads 少** | subclonal SNV read 本就少 → 跨點 read 更少 → 低 power（= 既有 carrier-limited 414 問題）| 報每對的有效跨點覆蓋；低於門檻標 under-powered |
| **ONT base error** | 定序錯誤造假 ALT/REF → **假 co-occurrence pattern** | 嚴格 QV 濾；考慮 error model；strand 一致性檢查 |
| **CN / multiplicity** | amplification 區一條 read 的解讀不同；multiplicity 影響 | 標 CN/LOH context（held-const）|
| **只 LOCAL** | read/block span 內才有「同分子」證據；**跨 phase-set 不可連結** | 紅線：禁宣稱 genome-wide lineage（§6）|
| **需 normal** | 分「獨立 subclone」vs「germline/artifact」需 normal 對照 | Phase 2 normal cis-control |
| **統計 underpowered** | co-occurrence 計數小 → Fisher/LD 不穩 | 慎報；多重檢定校正；不對單對 over-claim |

---

## §4 暴力先算 → 後驗證（兩階段流程，直接回答用戶）

### Phase 1 — 暴力枚舉（cheap / exhaustive / 無推論）✅ 可先做
列所有可跨對 → co-occurrence 表 → 分類 → 輸出「**co-occurrence event catalog**」（TSV/JSON）。這階段**只描述**，不下 subclone 結論。

### Phase 2 — 驗證與觀察（對 Phase 1 的候選事件）
對每個候選「同分子 subclone」事件（尤其互斥 / 巢狀）：
1. **normal cis-control**：normal 該位點有沒有此 ALT / 此甲基 pattern？→ 排 germline / artifact。
2. **甲基（兩步，🔴 不可跳第一步）**：
   - **Step 2a 先扣 cis**：「甲基隨 SNV 組合分開」**不能直接當佐證** —— 它可能只是 **cis-ASM**（SNV 在局部直接改甲基 = 同一個突變事件被看兩次，**非獨立證據**）。先過 **normal cis-control + 看離 SNV 較遠的獨立位點**。
   - **Step 2b 才算佐證**：**只有扣掉 cis 後甲基仍與 lineage 共分離**，才算 SNV 之外的**獨立 lineage 佐證**（corroborate genetic 骨幹，非甲基驅動）；否則標 **cis-ASM（同事件，不計入獨立證據）**。
3. **CCF gradient / 跨樣本**：是否有 CCF 階梯（chr2:18M 法）；6 cell-line 是否復現。
4. **誠實 framing**：每事件標證據級（same-molecule > same-haplotype）+ tier；catalog = characterization 非 confirmation。

> **核心節奏（回答用戶「暴力先算後驗證」）**：✅ **對** —— Phase 1 暴力先建 catalog（無推論、可重現），Phase 2 才對候選逐一過 normal + 甲基 + tier。**先算後驗，不要在枚舉階段就下 subclone 結論。**

---

## §5 與 ISM 既有的關係（接手前先對齊）
- **chr2:18M 旗艦 = 手動先例**（6 sSNV，α/β 互斥 0 違反，LOH HP1 1.4%）→ 本任務 = 系統化它；**Phase 1 工具須能重現 chr2:18M 手動結果作 sanity check**。
- **Foltz 2024 = genetic-only 先驗**（linked-read）→ ISM 增量 = **原生 ONT 單分子 + 甲基共分離軸 + normal cis-test**。引用時對照（must-distinguish）。
- **甲基角色** = **佐證**：同分子 SNV co-occurrence 是骨幹，甲基看「是否也跟著 SNV 組合分」→ corroborate，**非 driver**（紅線）。

---

## §6 🔴 紅線 + 給其他 AI 的處理指引
1. **只 LOCAL（read/block span）**：禁宣稱 genome-wide lineage / 完整 tree；跨 phase-set 物理不可連結。
2. **co-occurrence ≠ confirmed subclone**：需 normal 排 germline/artifact + 慎防 ONT error；單樣本只能 characterize（gold standard 仍是 single-cell/multi-region）。
2b. **🔴🔴 甲基共分離 ≠ 佐證（baseline-dependence 鐵則）**：**「甲基分群/共分離」永遠只是中間產物，不是結果；其意義 = f(分群, 你拿它對照的 baseline)**。甲基隨 somatic SNV 共分離在 **cis 基準**下 = cis-ASM（同事件、非獨立證據），須先過 **normal cis-control + 遠端獨立位點** 扣掉 cis，**只有殘餘共分離才算獨立 lineage 佐證**。**禁把「甲基跟著 X 分開」直接寫成佐證/結果而不講對照誰**（cis vs lineage vs germline-allelic 是不同 baseline → 不同結論）。同理：對齊 a-priori 標籤 ≠ subclone（85% = cis-ASM）；PERMANOVA 有結構 ≠ subclone；切 k 群依 null 而定。
3. **same-molecule vs same-haplotype 兩級**：read 直接跨 = 強；haplotag 延伸 = 弱（依賴 phasing）。明標。
4. **§13 數據誠信**：所有 co-occurrence 數字**先落檔（TSV/JSON）→ Read 讀回 → 才寫報告**；產數據與寫報告**不同批**。
5. **紅線一致**：reconstruction 由 somatic haplotagging 驅動、甲基佐證、regional proof-of-concept、⭐3 單樣本。
6. **error/統計**：低跨點覆蓋標 under-powered；多重檢定校正；不對單對 over-claim。

---

## §7 驗收標準（其他 AI 完成的判準）
- **Phase 1**：對給定 region 輸出每對 co-occurrence 表 + 分類 + 統計 + HP/strand/coverage，落 TSV/JSON；**重現 chr2:18M 手動結果**（α/β 互斥）作 sanity check。
- **Phase 2**：候選事件過 normal cis-control + 甲基共分離 + tier 標註；誠實 framing（LOCAL、characterization 非 confirmation）。
- **文件**：方法 + 困難 + 限制 + 紅線寫入報告；數字溯源（§13）。

---

*關聯 memory：`project_chr2_18m_subclone_locus_verification`（手動先例）· `project_subclone_confirmation_construction_ont_transferability`（V4 定位 + 移植裁決）。*
*關聯卡：`external_validation/axis2/foltz-somatichaplotype-linkedread-2024`（genetic-only 先驗）。*
*紅線權威：`InterSubMod/docs/method_comparison/20260621_subclone_confirmation_construction_and_ONT_transferability_01.md`。*
