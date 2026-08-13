<!--
建立時間: 2026-08-13
目標: 取代不可攜 absolute symlink，保留外部衍生資料的可追溯 pointer
處理範圍: docs/methodology/_assets/20260618_subcluster_pilot/obs_ws/cpp_wg/figs
關聯檔案: InterSubMod/research/20260813_complete_research_handoff/receipts/repo_hygiene_archive_manifest_20260813.json
驗證方式: 讀取同目錄 ARTIFACT_POINTER.json 的 path/count/size/metadata digest
證據等級: L3（provenance pointer；不含科研 payload）
-->

# 外部衍生資料指標

Whole-genome cpp_wg derived visualization directory.

- 原始本機路徑：`/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/figs_cpp_wg_full`
- 可用性：`LOCAL_DERIVED_UNTRACKED`
- 科學 claim ceiling：PARTIAL visualization evidence; not an authority artifact or biological validation.
- Git 政策：此目錄只保存 pointer 與 inventory；不複製大型衍生內容。

機器可讀的檔數、總大小及 metadata digest 請見 `ARTIFACT_POINTER.json`。該 digest
只涵蓋相對路徑、entry type 與檔案大小，不是逐檔內容 SHA-256。
