<!--
建立時間: 2026-07-12
目標: 提供碩士論文整合初稿交付物索引、狀態與閱讀順序
處理範圍: 本目錄 4 份主文件、1 份 HTML、6 張 SVG 與圖組索引
關聯檔案: 20260712_碩士論文重點流程計劃書_01.md；20260712_碩士論文_01.md
-->

# 碩士論文整合初稿索引

> **狀態：可審閱整合初稿／PARTIAL EVIDENCE。** Raw-all producer 已 7/7 完成；clean layered-v3 的 7 個 workers 均 exit 0，但 final verifier 於 2026-07-12 12:27:13 +08:00 exit 7，terminal state=`FAILED`。因此本交付不將未通過驗證的 topology outputs 或 historical rates 寫成正式 Results。

## 執行摘要

| 問題 | 目前結論 |
|---|---|
| 論文核心是什麼？ | 以 somatic haplotagging 與同一 DNA 分子上的 sSNV 共現，重建**可稽核的區域 mutation-state candidate-tree set**。 |
| 已完成什麼？ | 題名、六章論證、演算法適用域、資料關係、口試主線、六張圖與離線 HTML 均已建立；producer evidence 可正式引用，solver 只保留為有明示限制的 historical audit。 |
| 還缺什麼？ | 全樣本 downstream verifier 尚未通過，沒有可 release 的 funnel、determinacy 或 topology rates。 |
| 現在能否正式口試／送審？ | **科學結果凍結仍為 NO-GO。** 本版可供指導教授審稿與口試演練；需先修復 verifier、取得 `_SUCCESS` 與 U0–U7 PASS，再判定 submission readiness。 |
| 最高可說到哪裡？ | Regional hypothesis generator；不能說 confirmed biological clone genealogy、patient generalization 或 methylation-confirmed lineage。 |

建議正式題名為：**奈米孔長讀長資料之可稽核區域突變狀態候選樹重建：整合體細胞單倍型標記與單分子突變共現**。甲基化保留在相關研究、負向結果與有界輔助證據，不放入主標題。

## 建議閱讀順序

1. [重點流程計劃書](20260712_碩士論文重點流程計劃書_01.md)：先看研究主軸、pre-decision audit、scope 與完成 gates。
2. [論文架構與重點](20260712_碩士論文架構與重點_01.md)：確認六章 assertion-evidence、claim registry 與圖表配置。
3. [口試重點敘述講稿](20260712_口試重點敘述講稿_01.md)：12 張主投影片、22–25 分鐘逐頁講稿與常見問答。
4. [碩士論文 Markdown](20260712_碩士論文_01.md)：可 diff、可繼續修改的 source draft。
5. [碩士論文 HTML 摘要檢視版](20260712_碩士論文_01.html)：離線閱讀、sticky TOC 與列印版面；完整章節與附錄以 Markdown 為準。
6. [圖組索引](figures/README.md)：6 張原創 SVG 的用途、來源與更新條件。

## 本次新增交付

| 類型 | 數量 | 備註 |
|---|---:|---|
| Markdown 主文件 | 4 | 計劃、架構、口試、全文 |
| HTML | 1 | 濃縮摘要檢視 companion；零 CDN、零外部 JavaScript |
| SVG | 6 | 全部位於 `figures/`，相對連結 |
| 圖組索引 | 1 | 明示資料來源與更新 gate |

## 下一個必要動作

修正 verifier 報告的 `E_SCHEMA_INVALID` 與 `E_POST_INPUT_IDENTITY`，重新取得 clean layered-v3 `_SUCCESS` 且 U0–U7 全部通過後，才解除正式口試／送審的科學結果 NO-GO，並更新：

- 論文第四章 4.6 與摘要／結論；
- `figures/04_current_evidence_snapshot.svg`；
- 全樣本 funnel、canonical/sensitivity 與 bounded auxiliary-result 終版圖；
- 口試投影片 8–12 的數字；
- number provenance 與 immutable hash table。
