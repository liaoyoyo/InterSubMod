<!--
建立時間: 2026-07-25 07:00
目標: 記錄 ALT/REF 甲基證據整合 exact-PS layered workstation 的設計決定、折衷與驗證
處理範圍: 資料 authority、sample/index HTML、候選樹 overlay、Chromium regression
關聯檔案:
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/pre-decision-audit.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
-->

# ALT/REF 甲基證據 × exact-PS 工作站 implementation notes

狀態：validated
Spec：`methyl-alt-ref-topology-overlay`
上游 audit：GO 70/100

## 🔵 設計決定

### D-01 — 三種 grain 分開

- **Status:** Accepted
- **Tier:** L1
- **Decision:** We will preserve `formal pair`、`focal ALT/REF site control`、`exact-PS×HP lane` as separate payload levels.
- **Rationale:** focal control 的 joint clustering 與 formal G1 ALT-only methyl groups 不是同一 read population。
- **Revisit if:** 上游產生單一 receipt-bound normalized schema。

### D-02 — Overlay 不改 topology authority

- **Status:** Accepted
- **Tier:** L1
- **Decision:** We will render methyl evidence as annotation/filter/visual overlay only; it will not alter candidates, AF score, determinacy, or selected exemplar.
- **Rationale:** 目前證據是 read-level association 與 model-conditional relation，不是 topology objective。
- **Revisit if:** 有預先註冊、獨立驗證的 joint methyl-topology objective。

### D-03 — 缺少 overlay 不等於陰性

- **Status:** Accepted
- **Tier:** L1
- **Decision:** We will label non-target rows `NOT_IN_FORMAL_OVERLAY` and explicitly state that this is not a methyl-negative call.
- **Rationale:** 本次只有 7 個 formal-positive pairs，不是全 98,955 groups 的同口徑普查。
- **Revisit if:** 全 cohort 同 gate 的 positive/negative sidecar 完成。

### D-04 — Direct read support 拆分 signal 與 RR-only background

- **Status:** Accepted
- **Tier:** L1
- **Decision:** 首頁與樣本頁分別顯示 638 formal-signal reads 與 152 paired RR-only background reads；790 只保留為守恆總量。
- **Rationale:** 雖然原介面已標示 7 signal／3 background lanes，但把 790 當主數仍可能被誤讀成全部都是關聯訊號。
- **Revisit if:** lane role 或 formal-pair selection contract 改版。

## 🟠 偏離

### DV-01 — 原 validation report 綁 all7_v1

- **Status:** Closed
- **Tier:** L1
- **Decision:** Report 已改綁 current `all7_v2` 並重算；10/10 lanes 均以 exact `sample+region_id` 對到目前 topology，current source SHA 與 receipt-bound SHA byte-identical。
- **Result:** 7 pairs／10 lanes 的 active loci、resolution、ties、candidate relation 與重算結果一致；H2009 chr5 background lane 明確標為 `PAIR_NOT_ACTIVE_IN_GLOBAL_BEST`。
- **Revisit if:** `all7_v2/receipt.json`、factorization source 或任一 join TSV identity 改變。

### DV-02 — Portable HTML delivery receipt 未由 delivery 腳本同步

- **Status:** Closed
- **Tier:** L1
- **Decision:** `deliver_portable_with_topbar_fix.mjs` 現在只在 32 blocks／6 charts／5 metrics／6 tables、1440/390 viewport 與 source dialog 全 PASS 後，重寫 output HTML 與 `html_delivery_receipt.json`。
- **Result:** receipt schema 1.1.0 的 artifact SHA、HTML SHA 與 size 均逐位元組對到目前檔案；不再沿用舊封裝 receipt。
- **Revisit if:** portable builder 或驗證器版本改變。

## 🟡 折衷

### T-01 — Pair-level marker 展開到多個 HP lanes

- **Status:** Accepted
- **Tier:** L2
- **Decision:** We will make every matched lane discoverable, while the detail panel labels focal-control statistics as site-level and strict/candidate fields as lane-level.
- **Rationale:** 隱藏 RR-only background lane 會失去 `W_containers` 與 exact PS×HP 完整脈絡；不分層則會誤導。
- **Revisit if:** UI 能在 genome view 原生顯示 pair-level bracket，而非 group points。

## 🔴 未決問題

### Q-01 — 是否擴成全 cohort 甲基普查

- **Status:** Open
- **Tier:** L5
- **Question:** 是否要把 upstream 407,738 evaluated rows 轉成所有 final groups 都有明確 `evaluated/pass/fail/not-testable` 的 sidecar？
- **本輪處理:** 不擴 scope；只把 7 formal positives 清楚標示。

## 📚 Lore

- H2009 chr4:2,307,521 是唯一 focal ALT/REF `ALIGNED`，但其 mutation-bearing HP2 topology 因 resource guard ABSTAIN；不能畫成已解析 candidate branch。
- H2009 chr18 full tree 有 6 個 AF-best ties，但 focal→partner relation 在 6/6 trees 一致；局部 relation 可解不等於 full tree unique。
- HCC1395 focal ALT=98、REF=2，focal ALT/REF control 不可檢定；其可標的是 focal-ALT methyl group 與 linked-partner allele 的 read-level enrichment。

## ✅ 完成與驗證

### 產出

- `InterSubMod/docs/methodology/_assets/layered_workstation/index.html`
- `InterSubMod/docs/methodology/_assets/layered_workstation/{7 datasets}.html`
- `InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/20260725_ALTREF甲基差異與latest候選拓撲對應驗證_01.html`
- `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v4_all7_09/receipt.json`

### 執行命令

```bash
python3 research/20260725_methyl_alt_ref_topology_overlay_validation/scripts/build_alt_ref_topology_report.py
python3 research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
python3 -m unittest -v \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_layered_workstation.py
python3 research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py \
  --output-dir research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v4_all7_09
```

### 實際驗證片段

- sidecar：`all_pass=true`、32 checks、7 pairs、10 lanes、638 signal + 152 RR-only background = 790 direct reads、3 resolved pair relations、1 focal-axis aligned。
- renderer：`samples=7 groups=98955 ranked=71955`，輸出 7 sample pages + index。
- unit/regression：7/7 workstation tests + 28/28 gene/drug、candidate、cohort tests PASS。
- Chromium 148：16/16 page states、50 screenshots；無 console error、無 external request、無頁面水平 overflow，audited input SHA 前後一致。
- filter conservation：H2009 signal/aligned/resolved = 5/1/1；HCC1395_HKU = 1/0/1；HCC1954 = 1/0/1；其餘四個 dataset = 0/0/0。
- Claude Code `claude-opus-5` 唯讀對抗式複核：`VERDICT=PASS`、無 blocker；其 direct-support 拆分建議已由 D-04 落地，並另使 resolved filter 直接依 candidate relation 判斷，避免未來 aligned+resolved lane 被顯示優先序遮蔽。
