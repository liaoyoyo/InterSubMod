---
allowed-tools: Bash(cd:*), Bash(cmake:*), Bash(make:*)
description: 編譯專案 (Release 模式)
---

編譯 InterSubMod 專案。

## 執行步驟

```bash
cd /big8_disk/liaoyoyo2001/InterSubMod/build && make -j$(nproc)
```

## 如果需要重新配置

```bash
cd /big8_disk/liaoyoyo2001/InterSubMod
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

編譯完成後報告：
1. 編譯是否成功
2. 是否有警告
3. 可執行檔位置
