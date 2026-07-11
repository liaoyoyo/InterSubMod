<!--
建立時間: 2026-07-11
目標: 固定分層重建全景報告的 HTML renderer dependencies
處理範圍: Data Analytics v0.2.6 shell/runtime/embed/chart-contract source bundle
關聯檔案: ../scripts/build_layered_observation_report.py
-->

# Layered panorama report runtime bundle

此目錄讓 `build_layered_observation_report.py` 不再依賴會隨 plugin 更新而消失的 cache 路徑。2026-07-11 16:01，Data Analytics cache 由 v0.2.6 輪替到 v0.2.8；原 shell 與 embed helper 路徑因此失效。以下檔案已固定在本 cycle，並以原產物 provenance 與 deterministic rebuild 驗證。

| 檔案 | 角色 | SHA-256 |
|---|---|---|
| `html-report-shell.v0.2.6.html` | authored report shell | `1a524c1bce1b6e21715f9a1890119a1c705a2cf1d516f3c45650db3073085a56` |
| `embed_html_report_runtime.v0.2.6.py` | deterministic embed helper；改為讀本目錄 bundle | `fb76020ad214d34d094695db7584a90e01a9c7e9e32f54d353a3359dc34dcd5d` |
| `html-report-runtime-source.v0.2.6.txt` | gzip+base64 module runtime source | `5fe85fa510047937d131120e0d5061498fb4fdadb6543c361df441892eab259b` |
| `html-report-runtime-styles.v0.2.6.css` | runtime styles | `ee2a210983813cff554a7c60cc2fbe282d8aaaf28840f220551836fce4a4f9a0` |
| `html-report-source-tooltips.v0.2.6.js` | embedded source helper；本報告正文 tooltip 已停用 | `78532275cf179ae62f97088d542c42aa722db6e776410f4d4841451b1aa6c97f` |
| `chart_contract_v0_2_6.py` | official chart contract | `7d2f7cfdf3726bd5930970547c878d3809e455dffcadb57ac36a47637140fb86` |
| `package_utils.py` | chart contract dependency | `f419da66e0eb28dad2a40213faa6c9d09287c77708ea67bce940678ae34bb1e1` |

重生時使用：

```bash
python3 research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py \
  ... \
  --shell-template research/20260710_layered_reconstruction_v2/report_runtime/html-report-shell.v0.2.6.html \
  --embed-helper research/20260710_layered_reconstruction_v2/report_runtime/embed_html_report_runtime.v0.2.6.py
```

驗證要求：builder exit `0`；companion validation `artifact_validation_status=PASS`；data regression `70/70 PASS`；Chromium desktop/mobile `37/37 PASS`。此 bundle 只解決 renderer custody，不提高研究數據的 scientific evidence tier。
