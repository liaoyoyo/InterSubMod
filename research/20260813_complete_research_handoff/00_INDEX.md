<!--
建立時間: 2026-08-13 10:40
目標: InterSubMod／LongLineage 完整研究交接實作的 pre-registration、scope 與證據鏈入口
處理範圍: B Comprehensive validation + D External handoff；不重跑 TiB science
cycle_id: cycle_20260813-1032-complete-research-handoff
status: in_progress
build_branch: agent/research-handoff-audit-evidence-20260813
build_commit: ddd8909a838318d8a77969313e9561c8ff9d01c2
worktree: /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813
data_sources: InterSubMod/research/20260813_complete_research_handoff/pre-decision-audit.md,InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json,InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
驗證方式: 每個 gate 以 machine-readable receipt、退出碼、hash、hostname 或 Git commit 驗證；沒有 receipt 即 blocked
證據等級: L2 ⭐⭐⭐⭐（handoff 工程進行中；release/public gate 未通過）
關聯檔案:
  - InterSubMod/research/20260813_complete_research_handoff/pre-decision-audit.md
  - InterSubMod/research/20260813_complete_research_handoff/implementation-notes.md
-->

# 完整研究交接與 GitHub 公開準備

> **PROBE / NOT RELEASE-READY**：本輪建立 research handoff snapshot；不把 frozen hash replay 誤稱 clean-source science rerun，也不在雙機、claims、license gates 前建立 release。

**Owner**: user + Codex
**Status**: in-progress
**Task types**: B + D
**Cycle**: `cycle_20260813-1032-complete-research-handoff`
**Serves**: G1、G3、G4、G5

## §1 Pre-registration

| H 預測（confirmatory） | 否證條件 | Decision threshold |
|---|---|---|
| H1 frozen authority 可在本輪保持逐位元一致 | 19 個登錄項任一 missing 或 mismatch | 19/19 MATCH 才可沿用 frozen authority；否則相關 claims 全部 blocked |
| H2 handoff registry 可消除「靠資料夾名稱猜 final」 | 任一 `FINAL_FOR_SCOPE` 缺 producer/input/command/schema/scope/hash/location | required fields 全通過 validator；否則該 artifact 降級 |
| H3 公開入口可收斂到同一 claim registry | 158 claims 中仍有 `NEEDS_CORRECTION` 或 `CONTRADICTED`，或 entrypoint SHA 分裂 | 158 rows 僅 `CONFIRMED`、`CONFIRMED_WITH_LIMITS` 或醒目 `UNVERIFIED` |
| H4 portable interface 能證明 bip7/bip8 可使用 | `--plan` 產生檔案，或任一真正 host 的 doctor/build/smoke/preflight 失敗 | 每台主機獨立 hostname receipt PASS；缺 bip8 receipt 就保持 blocked |
| H5 LongLineage 可達 public research preview 的法律與安全門檻 | 21 mappings 任一無合法 disposition，或 notices/SBOM/safety/CI 不通過 | 全 gate PASS 才能切 public；production `run` 仍維持 exit 6 與 NOT_READY |

**SMART**：

- [x] Specific：五個 gate 皆有可觀察輸入、輸出與禁止升格範圍。
- [x] Measurable：19 authority entries、158 claims、21 source mappings、兩個獨立 hostname receipts。
- [x] Achievable：不重跑 TiB science，只做 frozen replay、build/test、synthetic smoke 與 read-only preflight。
- [x] Relevant：直接服務 G1/G3/G4/G5。
- [x] Time-bound：本次 handoff implementation cycle 結束前產出 receipts；未完成項保留 BLOCKED，不延後改門檻。

## §2 Scope lock

<!-- BEGIN USER-SPECIFIED -->

- InterSubMod baseline 固定 `ddd8909a838318d8a77969313e9561c8ff9d01c2`。
- LongLineage preview candidate 固定 `b9aaa12a11fa00606bd174dabd0f172a5d112359`。
- 2026-08-01 authority manifest 與 denominator registry 維持科學數值 SoT。
- `confirmed cellular subclone=0`、`confirmed linear ancestry=0`；methylation association-only；CN/LOH not integrated。
- `88.2579%` 只表示 model-conditional graph shape，不是 accuracy 或 prevalence。
- `73afaeac`、08-13 drilldown/CNV、LongLineage 後續 commits 與 dirty work 不進本 snapshot。
- 不重寫 Git 歷史；只有真實 secret 才另案撤銷與清理。
- 未通過 gate 不建立 production release，也不把 research preview 升格。

<!-- END USER-SPECIFIED -->

## §3 Step → Verify

1. Freeze identity → 驗證：兩 worktree HEAD 精確匹配、dirty=0、include/exclude/dirty registries 存在。
2. Rebuild inventory → 驗證：18/19/35 與 missing runbook 均有明確 disposition。
3. Build registries/handoff → 驗證：schema validator、required-field、date-order、hash/finality checks 0 failure。
4. Hygiene → 驗證：tracked HEAD broken/absolute symlink、local settings、secret finding、raw payload皆 0。
5. Portable workflows → 驗證：`--plan` 前後 file census 相同；doctor 能定位 required inputs/index/tools。
6. Claims → 驗證：158 rows 無 contradicted/needs-correction；default README/Wiki/Pages source 同 commit。
7. LongLineage preview → 驗證：license/source mapping/notices/SBOM/safety/CI PASS；strict blockers與exit 6仍被正確辨識。
8. Release decision → 驗證：只有所有 blocking receipt PASS 才允許 tag/Release/public visibility。

## §4 Research-start five questions

- Thread D：相關；read-level epigenetic claim ceiling 必保留。
- Thread B retraction：相關；不得復活 whitelist/filter overclaim。
- KDE-corrected：本輪不重算數據；每個沿用 artifact 記錄其原資料版本，未知即 `UNVERIFIED`。
- VCF caller AF：需要描述 provenance；不得以 merged AF 代替 caller AF。
- 長計算/C++/搬移/NO-GO：不做 TiB 長算；會做 build/test、archive-first hygiene；公開發布 NO-GO 已鎖。

## §5 Reproducibility checklist

- [x] baseline commits fixed
- [x] authority and claim audit inputs identified
- [ ] include/exclude and dirty-work registries generated
- [ ] authority replay receipt committed
- [ ] environment/tool hashes recorded
- [ ] synthetic inputs and expected schema/hash committed
- [ ] bip7 independent receipt
- [ ] bip8 independent receipt
- [ ] InterSubMod CI receipt
- [ ] LongLineage preview CI/license receipt

## §6 Evidence state

| Gate | Current verdict | Evidence |
|---|---|---|
| clean worktree freeze | PASS | two exact HEADs; dirty=0 |
| frozen authority | PASS on bip7 data plane | fresh 19/19 hash replay；receipt待落檔 |
| handoff inventory | PROBE | 18/19/35 尚待逐項裁決 |
| public claims | BLOCKED | 58/158 need action；34 P0 |
| bip7 acceptance | PENDING | 尚未完整 build/smoke/preflight |
| bip8 acceptance | BLOCKED | 尚無 hostname=`bip8` receipt |
| LongLineage private preview work | GO | candidate clean；license audit待完成 |
| LongLineage public visibility | BLOCKED | source-origin/license/notices/SBOM 未過 |
| InterSubMod tag/Release | BLOCKED | release gates 未全過 |

## §7 Next spaced check

本 cycle 未結案前不設定 30-day recall；若以 blocked handoff 結束，建立 REFLECTION 與具體 reopen 條件。
