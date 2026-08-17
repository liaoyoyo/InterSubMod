<!--
建立時間: 2026-08-13 20:15 +08:00
目標: 讓無先前對話的新AI以固定六問驗收交接包，並輸出可機械驗證的JSON
處理範圍: package-only reader comprehension；不驗證TiB science rerun或live GitHub publication
關聯檔案:
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/ai_context/CONTEXT.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/schemas/reader-acceptance.schema.json
驗證方式: fresh agent、fork_turns=none、package-only回答；jsonschema + semantic rubric + current/commit dual-hash validator
證據等級: independent reader comprehension receipt；非science authority
狀態: ACTIVE_ACCEPTANCE_PROTOCOL_V2
-->

# Fresh Reader Acceptance Prompt

你是一個沒有先前對話、沒有隱藏研究記憶的新讀者。唯一允許的研究來源是指定 package root
內的檔案；不要讀 Git history、其他 repo 目錄、原對話或網路。你的工作不是相信摘要，而是依
package 內 evidence／registries 交叉核對後回答。

Package root：`InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/`

## 1. Provenance 與輸出規則

- 輸出單一 JSON object，不要 Markdown code fence、前言或後記。
- JSON 必須符合 `schemas/reader-acceptance.schema.json` v2.0.0；
  `evaluator_context` 固定為 `NO_PRIOR_CONVERSATION_PACKAGE_ONLY`。
- `tested_git_commit` 必須是完整 40 碼 SHA，可為目前 HEAD 或其 reachable ancestor。
- `package_source_manifest` 每列為 package-relative `path` 與該檔案 SHA-256。Validator
  會同時比對 current package bytes 與
  `git show <tested_git_commit>:docs/handoff/20260813_完整研究資料與軟體交接_01/<path>`；
  兩者任一不相符就不能 PASS。
- 每題與每個 prohibited-promotion 的 `evidence_paths` 都必須已出現在
  `package_source_manifest`，不得引用未 hash-bound 的路徑。
- 不得使用「請見文件」、「答案已保留邊界」等 generic answer。每題需實際寫出下述
  semantic anchors、證據與限制。
- 任一必要證據、語意 anchor 或 commit-bound hash 不足時，將 `verdict`設為
  `BLOCKED`，不得猜測或自行補值。

### package_source_manifest 固定最小集合

以下 15 個路徑全部必須加入；可另加其他實際引用的 package 檔案。
Reader receipt 與會在登錄 receipt 時變動的 `evidence/EVIDENCE_MANIFEST.json` 不列入固定集合，
避免 receipt 或 Evidence Manifest 自我 hash 的 circularity：

1. `00_INDEX.md`
2. `20260813_研究結論時間與Finality_01.md`
3. `20260813_軟體輸入輸出與研究流程_01.md`
4. `20260813_bip7_bip8操作與驗證_01.md`
5. `ai_context/CONTEXT.md`
6. `ai_context/READER_ACCEPTANCE_PROMPT.md`
7. `schemas/reader-acceptance.schema.json`
8. `registries/artifact_registry.json`
9. `registries/authority_superseded_crosswalk.json`
10. `registries/claim_registry.json`
11. `registries/machine_path_registry.json`
12. `evidence/authority_manifest.json`
13. `evidence/denominator_registry.tsv`
14. `evidence/authority_replay_receipt.json`
15. `evidence/longlineage_capability_matrix.md`

## 2. 六題與必要 semantic anchors

每題輸出 `question_id`、`status`、`answer`、`evidence_paths`、`limits`。依下列順序回答：

1. `Q_PROJECT`：寫出這是 ONT long-read 的 sSNV read linkage、candidate
   mutation-state reconstruction 與 methylation association **research handoff snapshot**；明說不是
   production release 或 cellular-lineage caller。
   必引：`00_INDEX.md`、`20260813_軟體輸入輸出與研究流程_01.md`。
2. `Q_CONCLUSION`：必須同時寫出 `confirmed cellular subclone = 0`、
   `confirmed linear ancestry = 0`、methylation association only、CN/LOH **not integrated
   into frozen reconstruction**，以及 `88.2579%` 只是 model-conditional graph shape，**不是
   accuracy 也不是 prevalence**。
   必引：`20260813_研究結論時間與Finality_01.md`、`evidence/denominator_registry.tsv`。
3. `Q_FINALITY`：必須解釋 `FINAL_FOR_SCOPE`、`evidence_status`、`finality`、
   `supersedes`，以及 final 僅在明示 scope 內成立，不等於 production-ready 或整體研究
   final。
   必引：`20260813_研究結論時間與Finality_01.md`、`registries/artifact_registry.json`、
   `registries/authority_superseded_crosswalk.json`。
4. `Q_SOFTWARE_ROLES`：必須分開兩條 parallel provenance chains：
   LongPhase-S/TO → HP/PS／phased/recalibrated VCF／tagged BAM；exact-PS/LongLineage →
   candidate families／read assignments。另說明 InterSubMod 產生 per-region methylation、distance、
   read clustering/statistics。Commit-pinned Python research solver 可是 science producer；
   validator/publication builder/HTML presenter 只驗證與呈現 validated data，不重算 science。
   必引：`20260813_軟體輸入輸出與研究流程_01.md`、
   `evidence/longlineage_capability_matrix.md`。
5. `Q_MACHINES`：必須同時寫 bip7 與 bip8；bip7 只有 bounded local preflight，fresh-clone
   全鏈仍有 blocker；bip8 沒有 host-local receipt，維持 `BLOCKED`。一台主機的 receipt 不能取代
   另一台。操作順序需含 doctor → build/test → synthetic smoke。
   必引：`20260813_bip7_bip8操作與驗證_01.md`、`registries/machine_path_registry.json`。
6. `Q_VERIFY_CONTINUE`：必須寫 19/19 authority **byte/hash replay MATCH**，並明示它不是
   full science rerun；再寫 registry schema/uniqueness validation、publication/release/machine blocking gates，
   以及開新 research cycle 前先做 pre-decision audit 並固定 commit、input hash、scope、
   denominator。
   必引：`00_INDEX.md`、`evidence/authority_replay_receipt.json`、
   `registries/claim_registry.json`。

## 3. 禁止升格檢查

每項輸出 `check_id`、`status`、具體 `explanation`與 `evidence_paths`。解釋不能只寫
「保留邊界」，必須包含下列 anchors：

- `NO_CELLULAR_PROMOTION`：0 confirmed cellular subclones；candidate/read groups 不是 cellular truth。
- `NO_ANCESTRY_PROMOTION`：0 confirmed linear ancestry；graph shape/read dendrogram/local block 不是 ancestry。
- `NO_882579_ACCURACY_OR_PREVALENCE`：88.2579%、model-conditional graph shape、not accuracy、not prevalence。
- `NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE`：7 technical datasets、6 biological IDs；不是7個獨立生物複本。
- `NO_FEATURE_AS_MAIN`：`b9aaa12` research preview/non-production；`P3/P4/P5/P7/P8 BLOCKED`；feature 不是 main/release。
- `NO_LOCAL_AS_LIVE_PUBLISHED`：local source 不是 live publication；main/Wiki/Pages 仍 BLOCKED；需 publish 後 re-fetch。
- `NO_BIP7_AS_BIP8`：bip7 不能代表 bip8；bip8 BLOCKED；需 host-local receipt。
- `NO_PYTHON_ROLE_CONFLATION`：commit-pinned Python solver 可作 science producer；validator/builder/HTML
  是 validator/presenter，不重算 science。

只有 jsonschema、六題 semantic rubric、八項 promotion rubric、manifest membership、current bytes hash、
tested-commit bytes hash 全部通過時，才能輸出 `verdict: "PASS"`。
