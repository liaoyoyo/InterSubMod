# 修正過程紀錄（changelog）資料契約

使用者常要「結果頁附上**修正的過程紀錄**」：哪些修正已落地、哪些待做、對結果影響、權威 vs 過時。本 skill 把它做成 **data-driven 的 `spec.changelog`**，由 generator 渲染成 §⓪b section + 頂部 banner。**禁止把 changelog 寫死在模板**。

## 結構

```json
"changelog": {
  "data_status":   {"summary": "本頁數據的真實狀態（如 pre-correction / ±1000bp / 凍結日期）", "src": "run.log 路徑"},
  "binary_status": {"summary": "程式/版本狀態（如 修正後 binary 已存在但未重跑）", "src": "binary mtime / git HEAD"},
  "phase":         {"current": "...", "next": "...", "decided": "日期+決策"},
  "corrections": [
    {"id": "SHORT-ID", "what": "修正內容", "status": "in-HEAD|planned-phase2|not-done|authoritative|superseded",
     "effect": "對結果的影響（標 🔴 若 meaning-changing）", "src": "verified provenance"}
  ],
  "audit_conclusion": {"summary": "方法稽核對這些修正的總結（是否翻轉方向）", "src": "audit doc 路徑"}
}
```

## status 語意

| status | 徽章 | 意思 |
|---|---|---|
| `in-HEAD` | ✓ 已落地(綠) | 修正已在當前 binary/碼，但**資料可能還沒反映**（`data_reflected:false`）|
| `planned-phase2` | ◔ 計畫中(琥珀) | 已決定要做、尚未跑（如 ±5000 重跑）|
| `not-done` | ✗ 尚未做(紅) | 已知問題、尚未修（如 Fisher→beta-binomial 需 C++ Hard Gate）|
| `authoritative` | ✓ 權威(綠) | 用於 authoritative-vs-draft：此為定本來源 |
| `superseded` | ✗ 已過時(紅) | 用於 authoritative-vs-draft：此版已被取代（如 chr2 報告A 的某數字）|

## src 驗證規則（§13.7：宣稱前先驗證）

寫每條 `src` 前**親驗存在**，否則降級或不寫：

- commit：`git log -1 --format='%h %ad %s' <c>` 確認存在 + `git merge-base --is-ancestor <c> HEAD && echo ancestor` 確認在當前分支。
- 文件：`ls <path>` 確認存在；引行號則開檔確認。
- 數字效應（如「AUC 0.64→0.835」）：來自 commit message / 文件，引用時標出處，不憑記憶。

## 兩種用法

1. **corrected-data refresh**（如 HCC1395 ASM）：data_status=pre/post-correction、corrections=程式/參數修正清單、phase=兩階段計畫。落地：`research/tsg_promoter_asm_reviewer/.../display_v2/correction_changelog.json` + script 102 §⓪b。
2. **authoritative-vs-draft**（如 chr2:18M）：corrections 用 `authoritative/superseded` 標「報告A 某數字 → 權威值」的校正表。落地：chr2 generator 的 `CONFLICTS` 區塊。

兩者都讓讀者一眼看到「這頁數字是哪個版本、修了什麼、還缺什麼」——避免把過時/未修的結果當定論。
