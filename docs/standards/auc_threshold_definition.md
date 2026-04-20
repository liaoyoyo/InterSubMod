# AUC Threshold Definition（AUC 0.58 門檻統一標準）

> **建立日期**: 2026-04-19（P2-B CHECKLIST 選 A）
>
> **狀態**: Canonical（本文件為全專案 AUC 門檻唯一定義）
>
> **前身**：分散於 `CURRENT_FOCUS.md`、`03_ISM分析價值界定.md`、`06_結論穩定性審查.md`、多個 audit cards 的 0.58 / 0.60 / 0.58-0.60 引用

---

## 1. 門檻定義

本專案所有 feature-level AUC 分析共享以下門檻系統：

| 門檻 | 數值 | 意義 | 使用場景 |
|------|------|------|----------|
| **極低信號** | AUC < 0.55 | 實質隨機，不具區分力 | TO QS=0.497、HP_Ratio 跨模式 r=0 等 |
| **弱信號** | 0.55 ≤ AUC < 0.58 | 邊緣信號，不可用於 filter | 多數 TO 甲基化特徵 |
| **實用信號（無 HP）** | 0.58 ≤ AUC < 0.65 | HP-free 特徵可接受下限；需 within-group 驗證 | HPFineNGroups residualized 0.617、Self-phasing read-level 特徵 |
| **中強信號** | 0.65 ≤ AUC < 0.80 | 有實用價值；需 cross-sample + BH-FDR 驗證 | LOSO AUC=0.721、LOH zone 內某些 CNV 特徵 |
| **強信號** | AUC ≥ 0.80 | filter 候選；需確認無 caller_af 洩漏 | caller_af=0.654（但非 ISM 特徵） |

### 1.1 為什麼是 0.58

- **生物學合理性**：1000× bootstrap 下 AUC 0.58 約對應 95% CI 下界 >0.54（7 樣本 cross-validation）
- **歷史一致性**：O11/O12/O13 residualized 驗證 threshold、R1-R5 評估、Option C 雙路 NEGATIVE 判定均以 0.58 為「是否值得深究」的門檻
- **confound control**：pooled OLS residualized AUC <0.58 視同「confound 解釋所有原始信號」

### 1.2 為什麼 0.60 仍偶爾出現

歷史報告（2026-03 前）使用 0.60 作為經驗 threshold。2026-04 後統一降為 0.58（因 within-group 驗證下許多邊緣信號落在 0.58-0.60 區間）。**新文件一律使用 0.58**。

---

## 2. 使用規範

### 2.1 報告撰寫

引用 AUC 門檻時：
- ✅ `AUC < 0.58` / `AUC ≥ 0.58` — canonical
- ✅ `AUC 0.58 實用信號下限`
- ⚠️ `AUC < 0.60` — 僅限歷史報告；新文件改 0.58
- ❌ `AUC 0.5-0.6 之間` — 模糊；應明示 0.55 / 0.58 邊界

### 2.2 audit card 引用格式

```markdown
**Effect size**: feature X AUC=0.617 (>0.58 實用下限 ✅；但 <0.65 中強門檻)
```

### 2.3 統計 pipeline 規範

```python
AUC_THRESHOLD_USEFUL = 0.58   # 全專案 canonical；見 docs/standards/auc_threshold_definition.md
AUC_THRESHOLD_STRONG = 0.65   # 中強信號下限
```

### 2.4 信號分類矩陣

| AUC | within-group 驗證 | 建議下一步 |
|-----|------------------|-----------|
| <0.55 | 不需 | 關閉方向 |
| 0.55-0.58 | 選擇性 | 評估是否跨樣本穩定 |
| 0.58-0.65 | **必做** | within-group + bootstrap CI + BH-FDR |
| 0.65-0.80 | **必做** + cross-validation | LOSO + external validation |
| ≥0.80 | **必做** + caller_af 洩漏檢查 | 確認後評估 filter 可行性 |

---

## 3. 例外與注意

### 3.1 不適用本門檻的情境

- **效應量度量**：Pearson r / Cohen's d 有各自門檻（Cohen 1988 convention）
- **分類器集成 AUC**：Voting / stacking 需單獨評估 overfitting（如 O12 Voting AUC=0.577 雖 <0.58 但說明 5 特徵全耗盡）
- **per-region AUC vs pooled AUC**：pooled AUC 高不代表 per-region AUC 高（Simpson's Paradox，如 cnLOH 表面 0.587 → per-sample mean 0.50）

### 3.2 門檻不可作為**唯一**判定依據

AUC 0.58 僅為「是否值得深究」的門檻，不等於「實用性」或「生物學意義」。例如：
- HPFineNGroups residualized AUC=0.617 >0.58 但**非 FP filter**（因 AUC<0.85），僅作 somatic heterogeneity marker
- Phase 1A 甲基化+context F1 增益雖小，但 CI 下界正值 → 結論為「proof-of-concept 方向」

---

## 4. 對照歷史結論的門檻使用

| 結論 | AUC | 門檻判定 | 意義 |
|------|------|----------|------|
| C03 TO AUC ceiling（60+ 特徵）| <0.64 | 未達中強 | NO-GO filter |
| C04 O11 residualized | 0.509-0.578 | <0.58 | confound 解釋所有信號 |
| C05 O12 L3 | <0.59 | 邊緣 | NEGATIVE |
| C07 G1-G7 | <0.64 | 未達中強 | NO-GO filter |
| C08 Read-level LOSO | 0.721 | 中強 | 但 FP removal=0% 失敗 |
| C14 TO QS | 0.497 | 極低 | 隨機 |
| C15 Option C HP-free combo | 0.564 | <0.58 | 全 HP 信號來源確認 |
| C16 HPFineNGroups residualized | 0.617 | >0.58 | 實用信號下限 ✅ |
| caller_af | 0.654 | 中強 | 但 ISM 無法超越 |

---

## 5. 維護責任

- 新結論 audit card 引用 AUC 門檻 → 鏈結本文件
- 發現 0.60 新引用 → 改 0.58 + 說明 canonical 來源
- 每季回顧本門檻是否需調整（如 within-group 發現 0.58 過寬，可提至 0.60）

---

## 關聯文件

- `docs/reports/audit/decisions/03_P2_precision_decisions.md#P2-B`
- `docs/reports/audit/decisions/CHECKLIST.md`
- `docs/reports/audit/cross_cutting/Pooled_OLS_Audit.md`
- `docs/CURRENT_FOCUS.md`
