---
id: ism-kb-02-samples-hcc1395-dorado
name: "HCC1395_DORADO Sample"
description: "HCC1395 的 Dorado basecall 變體；與 HCC1395 共用 truth set；TO ΔF1 僅 +0.000476（低信噪比）。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "Dorado variant paths and F1"
related_ids:
  - ism-kb-02-samples-index
  - ism-kb-02-samples-hcc1395
tags: [sample, hcc1395, dorado, basecall, variant]
canonical_paths: [02_samples/02_hcc1395_dorado.md]
alias_paths: []
---

# HCC1395_DORADO Sample

- 一句結論：HCC1395 的 Dorado basecall 變體；與 HCC1395 共用 truth；TO 模式信噪比低
- 適用對象：Basecall 版本比較、Dorado 相容性測試
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/
  ```

---

## 基本資訊

| 項目 | 內容 |
|------|------|
| 與 HCC1395 關係 | 同細胞株，不同 basecall 方法 |
| Platform | ONT Dorado basecaller |
| Truth Source | 共用 HCC1395 SEQC2 truth |
| TRUTH_TOTAL | 39,447 |

---

## F1 現況

| Pipeline | Δ F1 | 狀態 |
|----------|-------|------|
| paired_full | ~類似 HCC1395 | ✅ 完成 |
| paired_pileup | ~類似 HCC1395 | ✅ 完成 |
| TO | **+0.000476**（僅正增益但極低） | ⚠ 低信噪比 |

---

## Canonical Run

- `output/canonical/HCC1395_DORADO/paired_full/20260315_*/`
- `output/canonical/HCC1395_DORADO/paired_pileup/<timestamp>/`
- `output/canonical/HCC1395_DORADO/TO/<timestamp>/`

---

## 使用建議

- **Basecall 版本比較**：HCC1395（5kHz）vs HCC1395_DORADO 對比
- **勿與其他樣本混用**：因與 HCC1395 非獨立（同細胞株）

---

## 相關

- HCC1395 主樣本：[01_hcc1395.md](01_hcc1395.md)
- 索引：[00_index.md](00_index.md)
