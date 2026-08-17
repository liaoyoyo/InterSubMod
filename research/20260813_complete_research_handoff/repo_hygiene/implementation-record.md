<!--
建立時間: 2026-08-13 11:01 +08:00
目標: archive-first 修復 release worktree 的不可攜／壞 symlink 與 tracked local settings
處理範圍: InterSubMod release hygiene；Task Type B Comprehensive validation + D External handoff
關聯檔案:
  - InterSubMod/scripts/handoff/repo_hygiene.py
  - InterSubMod/tests/test_repo_hygiene.py
  - InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_before_20260813.json
  - InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_archive_manifest_20260813.json
  - InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_after_20260813.json
驗證方式: exact inventory preflight、外部 archive checksum、worktree symlink scan、restore round-trip unit test、git diff --check
證據等級: L1（repository state 與 byte checksum）；不構成科學結果驗證
狀態: PASS
-->

# Repo hygiene 實作紀錄

**TL;DR：archive-first 清理已完成 — release worktree 由 512 個 tracked symlink 中移除 501 個壞連結、將 7 個 repo 內連結改為有效相對路徑、將 4 個外部衍生資料連結改為 provenance pointer；目前為 0 absolute、0 broken，local settings ignore gate 亦已閉合。（影響：高，信心：高）**

用 Claim–Evidence–Verdict：

- **Claim**：本次 hygiene 的每個破壞性 worktree 動作都有可復原的 archive metadata；tracked local settings 原檔亦已在 Git 外保存。
- **Evidence**：外部 archive manifest SHA-256=`93017a2d146641004903e41ab27ff1e94812d36ee475e21457b41853562d965d`，512 筆 link metadata 完整，settings payload SHA-256=`db69b7d1f9d29dbfe64c8bb44842ed50a2f3f14f8efae6d36a2bba49ce01bd3f`；6/6 單元測試通過。
- **Verdict**：hygiene gate PASS；`.claude/settings.local.json` 與 machine-local site profile 已加入精確 ignore 規則，after receipt 的 `pass=true`。

## 1. 輸入與固定邊界

| 項目 | 實際值 |
|---|---|
| Release worktree | `/big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813` |
| 執行時 branch | `agent/research-handoff-repo-hygiene-20260813` |
| 執行時 HEAD | `e3c0889100da88e6056e9268e07f6f94f8e9312b` |
| Release baseline ancestor | `ddd8909a838318d8a77969313e9561c8ff9d01c2` |
| 原始 tracked symlink | 512 |
| `/proc/self/fd/*` broken symlink | 500，全部位於 `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/` |
| 其他 absolute symlink | 11：7 個 repo 內程式、4 個外部衍生資料目錄 |
| 其他 relative broken symlink | 1 |
| Tracked local settings | `InterSubMod/.claude/settings.local.json`，32,786 bytes；內容未寫入 Git receipt |

本次未讀取或改寫科研 payload，未逐檔計算大型圖片內容 hash，未 stage、未 commit、未 push，也未觸碰原始 dirty workspace `/big7_disk/liaoyoyo2001/InterSubMod`。

## 2. Step → Verify

1. **精確掃描 512 個 tracked symlink**
   → 驗證：`tracked=512`、`remove_proc=500`、`remove_broken=1`、`convert_relative=7`、`pointer=4`；未知 path／target／數量一律 fail closed。
2. **先建立 Git 外 archive**
   → 驗證：archive manifest 512 筆、settings 原檔 SHA-256 相符、archive directory mode=`0700`、settings mode=`0600`。
3. **按 manifest 修改 worktree**
   → 驗證：7 個 repo 內 target 改為有效相對 symlink；4 個外部 target 改成 pointer；501 個無可用 payload 的壞連結移除；local settings 從 worktree 移除並留下安全 example。
4. **重掃 prospective release tree**
   → 驗證：有效相對 symlink=7、absolute=0、broken=0、pointer directory=4、settings local present=false。
5. **可復原與格式驗證**
   → 驗證：6/6 unit tests PASS，涵蓋 archive→apply→restore 與 partial-apply recovery；`git diff --check` exit code=0。

## 3. 實際動作與結果

| Action | 數量 | 結果 |
|---|---:|---|
| `remove_proc_self_fd` | 500 | 僅保存 path、Git blob OID、原 target、availability 與 claim ceiling 後移除 |
| `remove_broken` | 1 | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/figures` 已移除；archive 記為 `MISSING/HISTORICAL` |
| `convert_to_relative` | 7 | `InterSubMod/research/tpfp_loh_af_kde_discrimination/scripts/` 的 7 個程式連結改指向 `../../ng_kde_rescaling/scripts/...`，target 內容 SHA-256 已記錄 |
| `replace_with_pointer` | 4 | 原 symlink 改成含 `README.md` 與 `ARTIFACT_POINTER.json` 的目錄；不複製外部 payload |
| remove tracked local settings | 1 | 原檔只保存於外部 private archive；Git 端加入 `InterSubMod/.claude/settings.local.example.json` |

### 外部衍生目錄 inventory

下表 digest 策略固定為：依排序後的 `entry_type + relative_path + file_size/link_target` 做 SHA-256；**不是逐檔內容 hash，不可用來宣稱科研內容 byte-identical**。

| Pointer path | Files | Dirs | Bytes | Metadata digest | Claim ceiling |
|---|---:|---:|---:|---|---|
| `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/obs_ws/cpp_wg/figs/` | 34,736 | 0 | 2,839,780,790 | `375f5ae4c359240f9e5b43bfabf3c3630b07fcd6f6eb675b7e7bb3e309ebc656` | PARTIAL visualization；非 authority／生物驗證 |
| `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/obs_ws/cpp_wg/figs_tn/` | 34,736 | 0 | 1,903,028,876 | `dbd72273860b776434dc0c179302dec5b9e8d8f3b3f3a0228ae770eebf2efdbb` | PARTIAL visualization；非 authority／生物驗證 |
| `InterSubMod/docs/reports/validated/2026/04/20260406_LOH雙定義交叉分析報告/figures/visual_inspection/` | 212 | 46 | 15,449,707 | `b11898eeea596e7b12660d51edbcde4a4a83a2b0bc46704ee2983de426d2f673` | derived visual inspection；需回查報告與 upstream evidence |
| `InterSubMod/docs/reports/validated/2026/04/purity_figures/` | 6 | 0 | 788,001 | `861bb144bf77ffd702f5e49b29b52adb08ba1bd6e546aeaa2952f5782d1998d9` | derived visualization；非 authority／獨立驗證 |

## 4. 執行命令與實際輸出

### Dry scan

```bash
python3 scripts/handoff/repo_hygiene.py scan \
  --repo-root /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813 \
  --expected-baseline ddd8909a838318d8a77969313e9561c8ff9d01c2
```

實際輸出：

```text
[RESULT] tracked_symlinks=512 actions={"convert_to_relative": 7, "remove_broken": 1, "remove_proc_self_fd": 500, "replace_with_pointer": 4, "tracked_symlinks": 512}
```

### Archive-first apply

```bash
python3 scripts/handoff/repo_hygiene.py apply \
  --repo-root /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813 \
  --expected-baseline ddd8909a838318d8a77969313e9561c8ff9d01c2 \
  --archive-root /big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene \
  --receipt-dir research/20260813_complete_research_handoff/receipts
```

實際輸出：

```text
[ARCHIVE] verified_manifest=/big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene/archive_manifest.json
[RESULT] pass=false pass_excluding_parent_owned_ignore_rule=true absolute=0 broken=0
[GATE] IGNORE_RULE_PENDING: .claude/settings.local.json is not ignored; parent integration owns .gitignore
```

上述為 archive-first apply 當下的中間結果；主 agent 隨後補上精確 ignore rule，並以 `verify` 重生最終 after receipt。

### Archive、單測與 diff 驗證

```bash
python3 scripts/handoff/repo_hygiene.py verify-archive \
  --archive-root /big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene
python3 scripts/handoff/repo_hygiene.py verify \
  --repo-root /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813 \
  --expected-baseline ddd8909a838318d8a77969313e9561c8ff9d01c2 \
  --archive-root /big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene
python3 -m unittest tests.test_repo_hygiene -v
git diff --check
```

實際輸出摘要：

```text
archive pass=true; symlink_record_count=512
verify pass=true; pass_excluding_parent_owned_ignore_rule=true
Ran 6 tests in 0.535s — OK
git diff --check exit code=0
```

## 5. Archive 與 receipts

| 類型 | 路徑 | SHA-256／狀態 |
|---|---|---|
| 外部 archive root | `/big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene/` | directory mode `0700` |
| 外部 archive manifest | `.../repo_hygiene/archive_manifest.json` | `93017a2d146641004903e41ab27ff1e94812d36ee475e21457b41853562d965d` |
| 外部 archived settings | `.../repo_hygiene/private/.claude/settings.local.json` | `db69b7d1f9d29dbfe64c8bb44842ed50a2f3f14f8efae6d36a2bba49ce01bd3f`；mode `0600` |
| Before receipt | `InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_before_20260813.json` | `f779c1539b7c08cacfb86ecd88f198ed336406895566ed6d285f9310c0109c9e` |
| Archive manifest receipt | `InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_archive_manifest_20260813.json` | 與外部 manifest byte-identical |
| After receipt | `InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_after_20260813.json` | `38124ad1f5a4d4ee9e0abff14870f56944aca93a51bb6845b660b7131c01a9ee`；`pass=true` |

緊急復原介面如下；只有確定要回復 reviewed pre-cleanup state 時才可執行：

```bash
python3 scripts/handoff/repo_hygiene.py restore \
  --repo-root /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813 \
  --archive-root /big7_disk/liaoyoyo2001/release_archive/InterSubMod/research-handoff-2026.08.1/repo_hygiene \
  --confirm-restore
```

restore 僅會處理 manifest 的 512 個 exact paths 與 settings；pointer 內若出現未預期檔案或內容被修改，會拒絕刪除並 fail closed。

## 6. 整合狀態

1. `.gitignore` 已加入 `.claude/settings.local.json` 與 machine-local site profile 的精確規則。
2. 最終 after receipt 已重生，`pass=true`、`errors=[]`。
3. 既有 ignore 規則會遮蔽 8 個 pointer files；整合時只精確 force-add 各 pointer directory 的 `README.md` 與 `ARTIFACT_POINTER.json`，不 force-add payload。
4. hygiene commit 後必須從 HEAD 再掃一次，確認 0 absolute／0 broken，而非只依賴工作樹 receipt。

## 7. 科學與發布 ceiling

這次結果只證明 repository hygiene、archive 可復原性與路徑可攜性，不重算 science、不提高 artifact finality，也不改寫 frozen authority。尤其不得由本紀錄推導 cellular lineage、subclone、methylation association、CN/LOH 或 88.2579% 的任何新結論。
