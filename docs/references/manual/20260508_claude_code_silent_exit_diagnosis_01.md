---
title: "Claude Code 啟動 silent exit 完整診斷與 workaround"
date: 2026-05-08
status: ops_note
tags: [claude-code, networking, ipv6, troubleshooting]
---

# Claude Code 啟動 silent exit 完整診斷

## §0 症狀

新啟動 `claude` 命令立即 silent exit 0：
- 沒任何 stdout/stderr
- exit code = 0
- 但已啟動的 session 仍能用

## §1 根因分析（2026-05-08 完整診斷）

### Layer 1：bun runtime IPv6 happy-eyeballs bug
- **GitHub issue**: anthropics/claude-code#20240 (close as dup of #19301 / #16521)
- **症狀**: bun-compiled `claude` binary 在 dual-stack DNS 環境嘗試 IPv6 → ENETUNREACH → 不 graceful fall back IPv4 → child exit_group(1) → parent silent exit_group(0)
- **strace 證據**:
  ```
  connect(IPv6 → 2607:6bc0::10:443) = ENETUNREACH
  connect(IPv4 → 160.79.104.10:443) = EINPROGRESS（沒等完）
  child exit_group(1)
  parent silent exit_group(0)
  ```

### Layer 2：claude 2.1.133 startup-after-mcp silent exit
- **症狀**: LD_PRELOAD 修 IPv6 後仍 silent exit
- **strace 證據**: read=378 / write=223（資料有傳輸）但 stdout/stderr 0 bytes
- **debug log**: 停在 `[claudeai-mcp] Fetched 4 servers` 之後
- **推測**: cloud MCP server response 處理或某個 startup check 失敗，binary 內部 catch 後 silent exit

## §2 已驗證 NOT 是元兇（消去清單）

| 假設 | 證據 |
|---|---|
| Plugin/skill 過多 | 移走 anthropic-skills + superpowers 5.0.7 仍 silent |
| github plugin token 缺失 | 移除後 ERROR 消失但仍 silent |
| 版本 binary 損壞 | 重灌 2.1.133 仍 silent；--version 正常 |
| OAuth 過期 | 重新 login（Sat May 9 05:04:17 expires）仍 silent |
| HOME 設定錯 | HOME=/bip7_disk/liaoyoyo2001（同活 session）仍 silent |
| Cloud MCP 卡住 | --mcp-config empty.json 仍 silent |
| Plugin/hook 干擾 | --bare 跳過全部仍 silent |
| network IPv4 不通 | curl IPv4 0.13s 完成 (401) ✓ |

## §3 Layer 1 Workaround（部分修復）

LD_PRELOAD 攔截 getaddrinfo，強制 IPv4-only：

```bash
# 編譯 .so
cat > /tmp/no_ipv6.c << 'CEOF'
#define _GNU_SOURCE
#include <stdio.h>
#include <netdb.h>
#include <sys/socket.h>
#include <dlfcn.h>
#include <string.h>
typedef int (*getaddrinfo_t)(const char *, const char *, const struct addrinfo *, struct addrinfo **);
int getaddrinfo(const char *node, const char *service, const struct addrinfo *hints, struct addrinfo **res) {
    static getaddrinfo_t real = NULL;
    if (!real) real = (getaddrinfo_t)dlsym(RTLD_NEXT, "getaddrinfo");
    struct addrinfo h;
    if (hints) {
        memcpy(&h, hints, sizeof(h));
        if (h.ai_family == AF_UNSPEC || h.ai_family == AF_INET6) h.ai_family = AF_INET;
        return real(node, service, &h, res);
    } else {
        memset(&h, 0, sizeof(h));
        h.ai_family = AF_INET;
        return real(node, service, &h, res);
    }
}
CEOF
gcc -shared -fPIC -o ~/.local/lib/no_ipv6.so /tmp/no_ipv6.c -ldl

# Bash alias
echo 'alias claude="LD_PRELOAD=$HOME/.local/lib/no_ipv6.so /bip7_disk/liaoyoyo2001/.local/bin/claude"' >> ~/.bashrc
source ~/.bashrc
```

**結果**：strace 從 `connect IPv6 ENETUNREACH` → 全 IPv4，子進程從 `exit_group(1)` → `exit_group(0)`，但 **claude 仍 silent exit 0**（Layer 2 未解）。

## §4 Layer 2 推測根因（無法本地驗證）

3 個候選：

1. **anthropic API endpoint 改動**（5/8 起 server-side response handling）
2. **bun runtime worker thread internal exception**（startup setup() 內某個 lifecycle 拋 exception → catch 後 silent exit，無 stack trace）
3. **claude binary 2.1.133 在當前 anthropic API + 機器網路狀態下的 race condition**

無法在用戶端解決。

## §5 為何 4 個活著的 session 還能用

| Time | Event |
|---|---|
| 4/25-5/8 早上 | 4 session 啟動成功，TLS keep-alive 建立 |
| 5/8 某時段 | IPv6 路由失效 / server-side 改動 → **新啟動的全卡** |
| 之後 | 4 session 持有舊 connection 持續可用 |

PID 284504 / 2399574 / 3624725 / 3625299，HOME=/bip7_disk/liaoyoyo2001，全部 `claude -[cr] --allow-dangerously-skip-permissions`。

## §6 Fallback 方案

| 方案 | 可行性 | 限制 |
|---|---|---|
| 用 4 個活 session 收尾工作 | ✅ 立即 | 全部 close 後就沒了 |
| VSCode Claude Code Extension | ❌ 同樣卡（已驗證 5/8）| - |
| https://claude.ai/code Web 版 | ✅ 不同 IPC channel | 開發體驗差 |
| 等 anthropic 修 + 重啟機器 | 🟡 時間不定 | issue 4 個版本沒修 |
| 換網路（mobile hotspot）| 🟡 issue 20240 user OK | 不便 |
| 寫 Anthropic SDK + API key 自架 wrapper | ✅ 確定可繞 | 失去 claude code 全部 IDE 功能 |

## §7 關鍵 references

- [GitHub anthropic/claude-code#20240](https://github.com/anthropics/claude-code/issues/20240) — IPv6 ECONNRESET
- [GitHub anthropic/claude-code#19301](https://github.com/anthropics/claude-code/issues/19301) — Starlink IPv6 latency
- [GitHub anthropic/claude-code#16521](https://github.com/anthropics/claude-code/issues/16521) — Master IPv6 issue
- [GitHub anthropic/claude-code#24988](https://github.com/anthropics/claude-code/issues/24988) — Silent exit 0 RHEL 8

## §8 變更歷史

| 日期 | 動作 |
|---|---|
| 2026-05-08 21:55 | 初版完成診斷與紀錄；確認 Layer 1 + Layer 2 雙層問題 |

## §9 後續驗證 — 已排除嫌疑（2026-05-08 下午）

| 嫌疑 | 動作 | 結果 |
|---|---|---|
| `.agents/skills/` 與 `.claude/skills/` 重複載入衝突 | mv `.agents/` 到 `/tmp/dot_agents_backup_*` | 仍 silent exit 29 行 |
| `.claude/agents/` 內 SKILL.md 誤放 | 檢查 frontmatter 正常為 agent 格式 | 不是元兇 |
| methodology-reviewer.md 用過期 model `claude-opus-4-6` | 應該 WARN 而非 abort | 不是元兇 |

## §10 終極結論（2026-05-08 22:30）

每次測試 debug log 都精準停在 line 29：
```
[DEBUG] [claudeai-mcp] Fetched 4 servers
```

無論：plugin 多寡 / mcp 配置 / hooks 啟用 / skills 數量 / HOME 路徑 / OAuth 狀態 / IPv6 / IPv4 強制 / .agents 移除 / --bare / --debug-file。

→ silent exit 在 setup() 完成「fetch cloud MCP servers」之後立即發生。
→ 這 step 屬於 binary 內 application-layer 邏輯，不在 syscall/network 可見範圍。
→ 本地無法繼續排查。

### 暫時對策

1. 用既有 4 個 active session（PID 284504 / 2399574 / 3624725 / 3625299）繼續工作
2. 等 anthropic 修 bug 或 server-side 自然恢復
3. 若需臨時不用 claude code，改用 `claude.ai/code` web 版

### 變更歷史

| 日期 | 動作 |
|---|---|
| 2026-05-08 22:30 | 追加 §9 §10：.agents/ 不是元兇；最終結論為 server-side / binary internal bug |
