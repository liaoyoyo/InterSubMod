# SUPERSEDED／HISTORICAL — lineage Python 參考原型

> **不可作為執行入口。** 這些原型只保留 source-origin／歷史 parity provenance。
> LongLineage private public-preview candidate `b9aaa12` 仍為 **NOT_READY／non-production**，
> P3/P4/P5/P7/P8 `BLOCKED`；「當時 parity」不是 current、跨資料集或 release 驗證。
> 現行操作入口請看 [`QUICKSTART.md`](../../../QUICKSTART.md) 與
> [`完整交接索引`](../../../docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md)。

這三個檔是 LongLineage 三支 C++ 工具改寫前的原型。**它們不在任何生產路徑上**，
保留的唯一理由是：C++ 版本當初是以它們為基準逐位元比對驗證的，
留著才查得到「C++ 的行為為什麼長這樣」。

| 原型 | 對應的 C++（LongLineage） | 當時的 parity 驗證結果 |
|---|---|---|
| `build_lineage_paths.py` | `apps/lineage_paths_main.cpp` | **逐位元相同，含輸出列順序** |
| `build_read_lineage_assignments.py` | `apps/read_assign_main.cpp` | 8 項統計全部相同；內容排序後 100% 相同 |
| `ll_bam_tag.py` | `apps/tag_bam_main.cpp` | 統計逐項相同（25,700 / 9,694 / 9,585 / 355 / 29 / 326） |

**一項刻意保留的差異**：當一條 read 落在多個 block 時，兩邊的輸出順序不同 ——
Python 走 dict 插入順序，C++ 走 block index 升冪。C++ 的規則不依賴查找順序，
結果更確定，因此沒有為了對齊 Python 而改回去。

## 為什麼放在 InterSubMod 而不是 LongLineage

LongLineage 的 `scripts/ci/check_source_boundaries.sh` 規定 **`.py` 只能出現在
`presentation/`**，用意是讓「資料處理路徑純 C++」這件事有機械保證，而不是靠自律。
把參考原型留在該 repo 會讓這道 gate 失去意義，因此三個檔移回它們原本的出處。

## 對照時的注意事項

- **C++ 才是權威。** 兩邊若有出入，以 C++ 為準；本目錄不隨 C++ 更新而同步。
- 逐位元一致性是**當時**驗證的結果（含輸出列順序），不是持續保證。
- 這些檔案不再維護，也不應該被拿來跑生產資料。

本檔不再提供可複製的 LongLineage command。任何 private-preview 實跑須先由 capability matrix
確認 revision、具名 blocker set 與資料授權，並使用該 revision 自帶的 fail-closed驗收；不得用本原型當 production parity 證據。
