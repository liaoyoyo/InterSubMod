<!--
建立時間: 2026-07-27
目標: 修正 layered_workstation chr1–22 等寬顯示，建立共用 GRCh38 bp 尺度與桌機/手機回歸證據
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/ 7 個 exact-PS 樣本頁
關聯檔案:
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_layered_workstation.py
  - InterSubMod/docs/methodology/_assets/layered_workstation/README.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10/receipt.json
-->

# layered_workstation GRCh38 共同尺度修正與驗證

## 任務分類與研究目標

- Task type：E（UI semantic bugfix）。
- Scope：7 datasets × chr1–22 × desktop/mobile；不是 subset validation。
- 服務 G4/G5：讓跨樣本頁面使用一致、可自動驗證且可被外部讀者正確解讀的基因組尺度。
- 資料結論不變：仍是 98,955 final groups、71,955 ranked units；本次只修正座標到畫布的映射與說明。

## 根因

GRCh38 染色體長度資料本身正確；錯誤在 Canvas renderer。舊公式把每個點除以
自己的染色體長度，且每條骨架都畫滿相同寬度：

```text
old point x = left + plot_width × midpoint / L_chr
old track width = plot_width
```

因此單一染色體內的相對位置正確，但 chr1–22 沒有共用 bp→px 比例；chr22
被放大到和 chr1 一樣寬。

## 修正契約

新的 `shared-grch38-bp-v1` 契約以 chr1 的 248,956,422 bp 作共同上限：

```text
track width(chr) = plot_width × L_chr / 248,956,422
point x(chr, midpoint) = left + plot_width × midpoint / 248,956,422
```

selected point、一般 region point 與染色體骨架全部使用同一公式。桌面顯示
0/50/100/150/200/249 Mb；手機顯示 0/125/249 Mb，避免窄畫面刻度互相擠壓。

![HCC1395 桌面完整 chr1–22 共同尺度](../20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10/HCC1395_desktop_genome_canvas_all.png)

![HCC1395 手機完整 chr1–22 共同尺度](../20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10/HCC1395_mobile_genome_canvas_all.png)

## Step → Verify

1. 修正 builder 與 UI contract
   → 驗證：7 個樣本頁均含
   `layered-workstation-exact-ps-v5` 與
   `shared-grch38-bp-v1`。
2. 加入跨面板 regression
   → 驗證：Python unittest 8/8 PASS；禁止舊
   `row.mid / data.chromLengths[row.chrom]` 公式。
3. 重建 7 個 standalone HTML
   → 驗證：builder exit 0，實際輸出
   `samples=7 groups=98955 ranked=71955`。
4. Chromium desktop/mobile 全量稽核
   → 驗證：7 樣本 × 2 viewport = 14/14 樣本狀態通過；
   64 張 PNG 與 receipt 完整。

## 執行命令與 authority drift 處理

直接執行 builder 時，fail-closed 正確攔到 LongLineage 工作樹內尚未提交的
`seed_incumbent` 修改。candidate receipt 綁定的是 LongLineage Git `HEAD`
版本；兩個 committed blob 的 size/SHA 與 receipt 完全一致：

- `obligation_bnb.hpp`：1,054 bytes，
  `1d19dfb9a1074a8a4cbae3d9a55d6f68b7d89fe02aa9e47a6c0021195a2676a0`
- `obligation_bnb.cpp`：17,975 bytes，
  `ba4924bf3b7de443cc29de725569e3f3ef6a5fa62ec481fe769e1e86423d17e9`

為避免 stash、restore、覆寫或誤簽新 receipt，重建程序用 `bwrap
--ro-bind-data` 將上述 Git HEAD 內容只在該 process 的 mount namespace
投影到 receipt 記錄的原路徑；主機上的 LongLineage 三個未提交檔案完全未動。
本頁仍綁 pre-seed `all7_v2`，不宣稱已驗證新的 seed patch。

```bash
bwrap --die-with-parent --bind / / \
  --ro-bind-data 3 /big7_disk/liaoyoyo2001/LongLineage/include/longlineage/solver/obligation_bnb.hpp \
  --ro-bind-data 4 /big7_disk/liaoyoyo2001/LongLineage/src/solver/obligation_bnb.cpp \
  --chdir /big7_disk/liaoyoyo2001/InterSubMod -- \
  env PYTHONDONTWRITEBYTECODE=1 python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py \
  3< <(git -C /big7_disk/liaoyoyo2001/LongLineage show HEAD:include/longlineage/solver/obligation_bnb.hpp) \
  4< <(git -C /big7_disk/liaoyoyo2001/LongLineage show HEAD:src/solver/obligation_bnb.cpp)
```

實際輸出片段：

```text
OK exact-PS authority verified: samples=7 groups=98955 ranked=71955
OUTPUT=/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/layered_workstation/index.html
OUTPUT_SAMPLE_PAGES=7
```

Chromium 稽核命令：

```bash
python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py \
  --output-dir \
  research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10
```

實際輸出片段：

```text
{"all_pass": true, "all_pages_audited": true,
 "no_console_errors": true, "no_external_requests": true,
 "no_horizontal_overflow": true, "audited_inputs_stable": true,
 "screenshots_recorded": true}
```

## 數值回歸結果

- 所有 14 個樣本 viewport 的
  `chr22_width / chr1_width = 0.2041259574336267`，精確等於
  `50,818,468 / 248,956,422`。
- 每個 viewport 都有 22 條 chr1–22 骨架。
- 197,910 個初始 rendered points（98,955 × 2 viewport）全部位於各自
  染色體骨架內。
- 無 console error、外部網路請求或頁面水平溢位。
- 最終 receipt：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10/receipt.json`。

## Claim ceiling

共同尺度只修正 GRCh38 物理位置與染色體相對長度的視覺語意；點大小仍只為
可見性，不代表 region span、read support、VAF、clone fraction 或生物重要性。
