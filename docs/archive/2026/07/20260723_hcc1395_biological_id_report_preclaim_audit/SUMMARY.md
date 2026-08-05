<!--
建立時間: 2026-07-23 17:35 Asia/Taipei
目標: 保存 HCC1395 biological-ID HTML 報告在 claim audit 前的版本與兩次 renderer QA 失敗截圖
處理範圍: 僅封存本輪新產生、隨後被主張稽核取代的 HTML 與 QA diagnostics
關聯檔案:
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/20260723_HCC1395_biological_ID_specificity與融合樹驗證_01.html
-->

# HCC1395 biological-ID report pre-claim-audit archive

## 封存原因

初版 HTML 的數值與圖表通過 canonical artifact validation，但文字將 bulk raw-VAF backbone 與局部無向 read-linkage 過度延伸為「主要 trunk」及「方法沒有重大問題」。主張稽核後，正式版改為：

- 支持 strong biological-ID specificity。
- 支持 candidate-site 與 bulk raw-VAF profile 高度相似。
- 支持 shared-opportunity 中可辨識局部無向連結高度非衝突。
- 不宣稱 DORADO strict directed topology、clone 數或 parent–child tree 已重現。

## 封存內容

- `20260723_HCC1395_biological_ID_specificity與融合樹驗證_01.html`：claim audit 前的可攜式 HTML。
- 兩張 `verification-failure.png`：shared portable reader top bar 因 `100vw` 超出 scrollbar-reduced client width 8 px 的 QA diagnostics。

正式修正版在原研究目錄使用相同 canonical artifact 流程重建，因此原使用路徑不變。
