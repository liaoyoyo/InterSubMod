<!--
建立時間：2026-07-24
目標：記錄 exact-PS C++ topology/read-AF research CLI 的可重現 build 與語意邊界
處理範圍：research-only；不修改 production LongLineage 或 InterSubMod target
關聯檔案：exact_ps_topology_af.cpp、../tests/test_exact_ps_topology_af.py
-->

# exact-PS C++ topology/read-AF CLI

## Build

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
  -I/big7_disk/liaoyoyo2001/LongLineage/include \
  InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp \
  /big7_disk/liaoyoyo2001/LongLineage/src/solver/obligation_bnb.cpp \
  /big7_disk/liaoyoyo2001/LongLineage/src/solver/parent_mapping.cpp \
  -ljansson -lcrypto -o /tmp/exact_ps_topology_af
```

## Run

```bash
/tmp/exact_ps_topology_af \
  --input /absolute/path/exact_ps_mlhp_SAMPLE.json \
  --output /absolute/path/topology.jsonl \
  --receipt /absolute/path/topology.receipt.json \
  --max-family-size 0 \
  --max-search-nodes 0
```

`0` 表示不設限。既有 output 或 receipt 不會被覆寫。JSONL 與 receipt
都先寫入同層 temporary sibling，完成後才 atomic rename。

## Exact scope

- active bits：完整或 partial supported pattern 中至少一次為 `A` 的原始位點。
- full `R/A`：mandatory vertex；partial `R/A/X`：subcube terminal group。
- structure：LongLineage recurrence-allowed exact obligation B&B，`k_active <= 12`。
- 固定 vertex set：每個 child 的 Hamming-1 parent 可獨立選擇。
- AF edge score：
  `sum_{i in parent} [AF_i - AF_acquired]`，
  `AF=ALT/(REF+ALT)`；全程以 `cpp_int` exact rational 計算。
- 全域 best tie count：同分 minimum vertex sets 的 per-child argmax
  parent-choice products 相加，不展開每棵完整樹。
- 若所有 minimum vertex sets 都不通過 rooted three-gamete test，
  標為 `recurrence_required` 且不做 AF ranking；若至少一組相容，則和舊
  Python ambiguous scope 一樣，排名全部 recurrence-allowed minimum candidates。

## Fail-closed boundaries

- top-level schema/JSON 錯誤：exit 1，不發布 output/receipt。
- malformed unit：保留明確 abstention row，receipt `all_pass=false`。
- family cap：可保留 certified `h*`，但不發布 partial family 或 AF ranking。
- search-node cap：objective/family 都 abstain。
- coverage 缺漏或分母為零：structure 可完整，AF ranking abstain。
- `representative_best_morphology` 只描述一棵 deterministic best
  representative，不是全部 top-tied trees 的 exact shape census。

## Tests

```bash
python3 \
  InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/tests/test_exact_ps_topology_af.py
```

測試用 Python exhaustive oracle 僅展開 `q<=4` 合成立方體，獨立核對
minimum family、tree count、AF best score 與 exact tie count。
