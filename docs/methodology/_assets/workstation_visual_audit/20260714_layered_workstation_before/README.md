<!--
建立時間: 2026-07-14
目標: 保存 layered workstation 改版前的完整桌面/手機視覺與互動基線
處理範圍: index + 7 sample pages；desktop 1440x1000 + mobile 390x844；非抽樣
關聯檔案: capture_layered_workstation_before.py, metrics.json, 20260714_layered_workstation_before_visual_audit_01.md
-->

# Layered workstation · before baseline

這是 Task B（comprehensive validation）的唯讀基線：8/8 HTML、16 個 page×viewport 狀態、32 張已逐張目視的 PNG。來源 HTML、generator 與 data 均未修改；前後 SHA-256 一致。

本任務主要服務 G4（跨樣本一致性與 reproducibility）與 G5（可被外部驗證的展示品質）。

## 輸入與輸出

- 輸入：`/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/layered_workstation/*.html`
- 擷取器：`InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py`
- 自動量測：`InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/metrics.json`
- 截圖：`InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/screenshots/`
- 逐張稽核：[20260714_layered_workstation_before_visual_audit_01.md](20260714_layered_workstation_before_visual_audit_01.md)

## 可重跑命令

```bash
python3 docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py --help

python3 docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py

# 只重算 DOM / interactions / hashes，不重拍 PNG
python3 docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_before/capture_layered_workstation_before.py --reuse-screenshots
```

最終實際輸出：

```json
{
  "status": "fail",
  "exit_code": 1,
  "screenshots": 32,
  "pages": 8,
  "page_viewports": 16,
  "checks_total": 289,
  "checks_failed": 16,
  "source_hashes_unchanged": true
}
```

`exit 1` 表示截圖與檢查已完成，但 before baseline 留有真實展示問題；不是 runner failure。`exit 2` 才代表執行器或設定失敗。

## 最終基線摘要

- 通過 273/289 checks；0 console error、0 page error、0 request failure、0 HTTP response error。
- 156 個 link / anchor / sample-switch target 全部存在；所有主要控制項通過。
- 76 次 sample details toggle/restore 全通過；16 個狀態的 heading hierarchy 無跳級。
- 鍵盤遍歷每頁 16–26 個目前可見 focus targets，缺少 focus indicator = 0；skip link 可見並可到達目標。
- 7/7 sample 的手機候選詳情表寬 551.4–693.6px，卻位於 390px viewport 中且沒有局部水平捲動，右欄被裁切。
- 7/7 sample 的手機樹圖被壓縮到 293px；估計節點標籤只剩 3.9–6.9px。
- index 已明示 region denominator，但 desktop/mobile 都沒有明確 `全基因組 / whole-genome` scope 標籤。

完整證據、32 steps 與優先改善順序見稽核文件。
