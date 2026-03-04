<!--
建立時間: 2026-03-01 03:12
目標: 固化 Knowledge MCP 在 InterSubMod 與其他專案的接入與驗證規範
處理範圍: .mcp.json 設定、連線驗證、提問時知識庫提醒
關聯檔案:
  - /big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py
  - .mcp.json
  - .claude/settings.local.json
  - scripts/hooks/knowledge_check.sh
-->

# Knowledge MCP 跨專案接入規範

- 生效日期：2026-03-01
- 適用對象：InterSubMod 與其他需要使用實驗室知識庫的專案

## 1. 其他專案直接接入範本

在目標專案根目錄的 `.mcp.json` 加入：

```json
{
  "mcpServers": {
    "knowledge": {
      "type": "stdio",
      "command": "python3",
      "args": ["/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py"]
    }
  }
}
```

## 2. InterSubMod 既有狀態

1. `.mcp.json` 已接入 `knowledge` server。
2. `knowledge://stats` 可讀，索引文件數為 35（active）。
3. `.claude/settings.local.json` 已啟用 `UserPromptSubmit` hook：
   - `bash scripts/hooks/knowledge_check.sh`

## 3. 最小驗證清單

1. MCP 資源可列出：可見 `knowledge://catalog`、`knowledge://stats`。
2. 文件可讀：`03_file_formats/vcf_clairs_to.md` 可成功讀取。
3. Prompt 自動提醒可觸發：輸入包含 `VCF`、`HCC1395`、`資料路徑` 等關鍵字時，hook 會提示應查閱 Knowledge 子目錄。

## 4. 使用原則

1. 涉及資料路徑、樣本規格、VCF/BAM 欄位定義、工具參數時，必須先查 Knowledge 再執行。
2. 回覆中需標註來源，例如：`根據 Knowledge/03_file_formats/vcf_clairs_to.md`。
3. 若發現 Knowledge 與專案文件不一致，先以 Knowledge 為準並回報文件差異。
