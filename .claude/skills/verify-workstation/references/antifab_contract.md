# §13-A 反捏造契約（verify-workstation）

承襲 CLAUDE.md §13 / §13.0 / §13.7。判讀工作站給人**肉眼判真偽**，所以它自己的數字**絕不能是捏造的**——否則人會基於假數字做判讀。

## 三條構造級保證

1. **注入，不手打**：每個顯示的 metric 都在 spec 的 `items[].metrics[].v`，由 generator 注入。generator 內**沒有任何寫死的 metric**。數字來自 W0 落檔的 `*.json/*.csv`，撰 spec 時 Read 讀回真值填入，並記 `src`。
2. **缺值即 refuse（不畫 dash）**：`item_config.required_metrics` 列必填欄；任一 item 缺 → `build_workstation.py` **exit 3 中止**。這正是修掉 page-04 Fig4「DORADO·β 該有值卻顯示『—』還說有複製」那類 v1/v2 漂移的構造手段。合法的 NA（如 n=1 無法算 FDR）要在 W0 明確標成字串 `"NA"`（真實狀態），不是「缺」。
3. **分類由資料算出，不手填**：CLEAN/CONFOUNDED、PASS/可能漏掉/篩掉、TierA/B… 一律在 W0 資料層依**可 grep 的判準**算出（例：normal HP1-vs-HP2 MW-FDR<0.05 ⇒ confounded；PERMANOVA gate ⇒ PASS）。spec 收到的是已解析結果。判準要寫進 `src` 或 changelog。

## 撰 spec 前自我檢查（全 ✓ 才跑 generator）

- [ ] 每個 metric **現在**都能在某個檔案 grep 到？
- [ ] 都是**這一輪 Read 回來**的真值（非記憶/預期）？
- [ ] 分析回傳 success（非 error/INCONCLUSIVE/檔案不存在）？
- [ ] 分類是**程式算的**（附判準 src），非肉眼/記憶指定？
- [ ] 跑 generator 的指令與產數字的分析**不在同一 tool-call batch**？
- [ ] changelog 每條 `src`（commit/doc）已親驗存在（commit 用 `git merge-base --is-ancestor <c> HEAD`）？

## 機械後盾（generator 之外）

- `scripts/number_provenance.py audit <html> [--sources ...]` — 抽報告內 metric 形數字逐一去 bounded 來源 grep（§13-C 溯源表）。
- `scripts/fill_report.py <template> <data.json>` — 另一條 template+data 注入 refuse-on-missing 管線（與本 generator 同精神）。
- hook `number_provenance_check.sh` — validated/pi 路徑寫入含無來源 metric → exit 2。

## 反例（禁止）

- 在 generator 或 spec 寫「預期值/合理範圍」。未齊 → 該 item 不收，或該 metric 標 `NA` 並在卡上可見。
- 把 AI 推測的分類當真值填 spec。推測 = L3，需實測升 L1 才寫定論（見 `feedback_researcher_claim_needs_empirical_verification`）。
- 跑分析的 Bash 與寫 spec/HTML 的 Write 同批送出（平行批讓 Write 拿不到當批數字 → 用記憶補 = 捏造根因）。
