<!--
建立時間: 2026-04-14 10:30
目標: 收集 2025-2026 年 CLAUDE.md / AI coding agent instruction 最新進階技術
處理範圍: CLAUDE.md 結構、Skill 觸發、Memory 管理、Hook 工程、Context 預算、Multi-agent 通訊
關聯檔案:
  - .claude/CLAUDE.md
  - .claude/settings.local.json
-->

# Prompt Harness Advanced Patterns — 2026-04 Survey

## 1. CLAUDE.md Structure Optimization

### 1.1 Maturity Level Framework (L0-L6)

| Level | Name | Core Technique |
|-------|------|---------------|
| L0 | Absent | No file |
| L1 | Basic | File exists, boilerplate |
| L2 | Scoped | Prescriptive rules ("MUST use X") |
| L3 | Structured | Multi-file split by concern; main file = routing doc |
| L4 | Abstracted | Path-scoped rule loading via `.claude/rules/` + frontmatter |
| L5 | Maintained | Staleness tracking, periodic reviews, backbone mapping |
| L6 | Adaptive | Dynamic skills + MCP; conditional capability loading |

Source: [DEV Community — From Basic to Adaptive](https://dev.to/cleverhoods/claudemd-best-practices-from-basic-to-adaptive-9lm)

**Our project is L3-L4.** Gaps: no `.claude/rules/` path-scoped loading, no staleness tracking.

### 1.2 The 150-Line Rule

Multiple sources converge on a hard ceiling: **150 lines per CLAUDE.md file**. Beyond this, compliance drops sharply because rules get buried in noise.

> "I wrote 200 lines of rules for Claude Code. It ignored them all." — [DEV Community](https://dev.to/minatoplanb/i-wrote-200-lines-of-rules-for-claude-code-it-ignored-them-all-4639)

**Our CLAUDE.md is ~400+ lines.** This is a concrete risk. The mitigation strategy:

1. **Keep CLAUDE.md lean** (<150 lines): project overview, build commands, top 5 critical rules
2. **Migrate domain rules** to `.claude/rules/` with path-scoped frontmatter
3. **Migrate enforceable rules** to hooks (see Section 4)
4. **Migrate reference info** to skills or linked docs

### 1.3 Path-Scoped Rules (`.claude/rules/`)

Files in `.claude/rules/` with YAML frontmatter load only when Claude works with matching files:

```yaml
---
globs:
  - "src/**/*.cpp"
  - "include/**/*.hpp"
---
# C++ Development Rules
- Use clang-format before commit
- All new functions need unit tests in tests/
```

**Known issue**: `paths:` frontmatter is documented but `globs:` actually works reliably. Path-scoped rules trigger on Read but NOT on Write/Create operations.

Source: [GitHub Issue #17204](https://github.com/anthropics/claude-code/issues/17204), [GitHub Issue #23478](https://github.com/anthropics/claude-code/issues/23478)

### 1.4 Compaction Protection

CLAUDE.md survives compaction (it loads as system prompt). Conversation context does NOT. Two strategies:

1. **Custom compaction prompt** via `~/.claude/settings.json` `"compactPrompt"` field — tell it what to preserve
2. **CLAUDE.md directive**: `"When compacting, always preserve: modified file list, test commands, current task state, all CLAUDE.md rules verbatim"`

Source: [Compaction Deep Dive](https://okhlopkov.com/claude-code-compaction-explained/)

---

## 2. Skill Triggering Reliability

### 2.1 Activation Rate Benchmarks

| Approach | Success Rate |
|----------|-------------|
| No optimization | ~20% |
| Optimized description | 50% |
| CLAUDE.md cross-reference | 60-70% |
| LLM pre-eval hook | 80% |
| Forced eval hook | **84%** |
| Structured evals + A/B testing | **94%** |

Source: [Skill Structure Guide (Gist)](https://gist.github.com/mellanon/50816550ecb5f3b239aa77eef7b8ed8d), [Towards AI — Evaluate and Tune](https://medium.com/@richardhightower/claude-code-how-to-build-evaluate-and-tune-ai-agent-skills-34afa808d1c9)

### 2.2 Description Writing Rules

1. **Third person only**: "Processes Excel files" not "I can help you"
2. **Include "USE WHEN" trigger**: explicit activation conditions
3. **Name file types and formats**: "Use when working with .cpp, .hpp, or CMakeLists.txt"
4. **Keep under 50 lines**: split long skills into focused sub-skills
5. **Include correct/incorrect examples**: reduces ambiguity for edge cases

### 2.3 Hook-Based Auto-Activation

When description optimization plateaus (~50%), use hooks for deterministic triggering:

```json
{
  "hooks": {
    "UserPromptSubmit": [{
      "matcher": "",
      "hooks": [{
        "type": "command",
        "command": "python3 .claude/hooks/skill-router.py"
      }]
    }]
  }
}
```

The router script checks current working files against a `skill-rules.json` mapping and injects skill references into context. This bypasses LLM semantic matching entirely.

Source: [paddo.dev — Skills Auto-Activation](https://paddo.dev/blog/claude-skills-hooks-solution/)

**Note**: Claude Code 2.0.64+ ships native path matching via `.claude/rules/`, partially obsoleting this hook pattern. However, hooks still provide more flexibility (e.g., multi-condition triggers, dynamic skill composition).

---

## 3. Memory System Optimization

### 3.1 Four-Layer Memory Architecture (Native)

| Layer | Storage | Persistence | Load Timing | Limit |
|-------|---------|-------------|-------------|-------|
| CLAUDE.md | Project root | Permanent (git) | Session start | ~4K chars/file, ~12K total |
| Auto Memory | `MEMORY.md` | Permanent (git) | Session start | 200 lines (index) |
| Session Memory | Conversation | Volatile | During session | Until compaction |
| AutoDream | Memory topic files | Permanent | On trigger | Consolidation cycle |

### 3.2 AutoDream Consolidation (4 Phases)

1. **Orient**: Read current MEMORY.md + topic files, build mental map
2. **Gather Signal**: Daily logs > drifted memories > JSONL session grep (never exhaustive)
3. **Consolidate**: Merge into existing topics, convert relative dates to absolute, delete contradicted facts
4. **Prune and Index**: Update MEMORY.md to <200 lines

Trigger conditions (ALL must be true): 24+ hours since last dream AND 5+ sessions AND no concurrent dream running.

Source: [dream-skill (GitHub)](https://github.com/grandamenium/dream-skill), [Anthropic system prompt](https://github.com/Piebald-AI/claude-code-system-prompts/blob/main/system-prompts/agent-prompt-dream-memory-consolidation.md)

### 3.3 External Memory Systems

**claude-mem** (46K+ GitHub stars): dual-database architecture

| Component | Technology | Purpose |
|-----------|-----------|---------|
| Structured store | SQLite + FTS5 | ACID transactions, complex filtering |
| Semantic store | ChromaDB + all-MiniLM-L6-v2 | Vector embeddings, semantic search |
| Compression | Claude Agent SDK (Haiku) | AI-driven observation summarization |

Key insight: local embeddings via ONNX avoid external API calls. Observations are compressed before storage, keeping token costs manageable.

Source: [claude-mem (GitHub)](https://github.com/thedotmack/claude-mem)

**memsearch**: Milvus-based alternative with persistent vector memory and plugin architecture.

Source: [Milvus Blog](https://milvus.io/blog/adding-persistent-memory-to-claude-code-with-the-lightweight-memsearch-plugin.md)

### 3.4 Our Memory Assessment

Our MEMORY.md index is well-structured (Active / Pending / Concluded sections). Missing:
- No AutoDream or equivalent consolidation cycle
- No staleness tracking on concluded items
- No semantic search across session history

---

## 4. Hook Engineering Patterns

### 4.1 The Rule Migration Principle

> "If Claude already does it correctly, delete it from CLAUDE.md. If a single violation costs significant time, build a hook."

| Rule Type | Where | Why |
|-----------|-------|-----|
| Broad preferences, style | CLAUDE.md | Low cost of occasional non-compliance |
| Build-before-commit, no rm -rf | **Hooks (blocking)** | Single violation = hours of rework |
| Context injection per directory | **Hooks (context)** | Deterministic, not subject to LLM judgment |
| Session protocols, memory search | **Hooks (lifecycle)** | Must run every time, zero exceptions |

Source: [DEV Community — 500 Lines of Rules](https://dev.to/mikeadolan/i-wrote-500-lines-of-rules-for-claude-code-heres-how-i-made-it-actually-follow-them-3c8)

### 4.2 Context Injection Pattern

```python
# .claude/hooks/inject-rules.py
# Triggered by PreToolUse on Read|Edit|Write
import json, sys
tool_input = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
file_path = tool_input.get("file_path", "")

PACKAGE_RULES = {
    "src/core/": ".claude/rules/core-rules.md",
    "tests/": ".claude/rules/test-rules.md",
}
for pkg_path, rules_path in PACKAGE_RULES.items():
    if pkg_path in file_path:
        rules = open(rules_path).read()
        print(json.dumps({
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "additionalContext": f"Rules for {pkg_path}:\n{rules}"
            }
        }))
        break
```

Hook output is capped at **10,000 characters**. Exceeding this saves to file with a preview.

Source: [DEV Community — Guaranteed Context Injection](https://dev.to/sasha_podles/claude-code-using-hooks-for-guaranteed-context-injection-2jg)

### 4.3 Advanced Hook Patterns

| Pattern | Implementation | Use Case |
|---------|---------------|----------|
| StatusLine monitoring | Rust binary tracking context % | Alert at 30%/15%/5% remaining |
| Pre-compact backup | `PreCompact` hook saves session state | Prevents information loss |
| Frustration detection | NLP on user messages | Escalate or pause when frustration detected |
| Session restore | Analyze git diff + session logs | Multi-factor context recovery |
| Circuit breaker | Count consecutive tool failures | Prevent infinite retry loops |

**Parallel hook warning**: When multiple PreToolUse hooks return `updatedInput`, the last to finish wins (non-deterministic). Avoid multiple hooks modifying the same tool's input.

### 4.4 Our Hook Assessment

Our hooks cover: build-before-commit (blocking), dangerous commands (warning), C++ edit tracking (state). Missing:
- No context injection hooks (directory-aware rules)
- No pre-compact protection
- No context usage monitoring (StatusLine)
- No skill auto-activation routing

---

## 5. Context Window Budget Management

### 5.1 Token Budget Anatomy (200K Window)

| Component | Tokens | Notes |
|-----------|--------|-------|
| System prompt | ~2.7K | Fixed |
| System tools | ~16.8K | Pre-loaded |
| CLAUDE.md + rules | ~4-12K | Scales with file count |
| Autocompact buffer | **33K** | Hardcoded, non-configurable |
| Usable conversation | **~134-143K** | Grows until compaction |

Compaction triggers at ~83.5% usage (~167K tokens used). With 1M context (Opus 4.6 GA), the economics change dramatically: no pricing premium, buffer is proportionally smaller.

### 5.2 Optimization Techniques

1. **`CLAUDE_AUTOCOMPACT_PCT_OVERRIDE`**: Set 70-90% to control compaction timing (lower = more aggressive)
2. **Strategic manual `/compact`**: Compact after completing major features, before starting new components
3. **Progressive disclosure**: Load detailed docs via skills/links, not inline in CLAUDE.md
4. **Path-scoped rules**: Load only relevant rules per directory (Section 1.3)
5. **1M context window**: `sonnet[1m]` or `opus[1m]` eliminates 200K constraint entirely

Source: [ClaudeFast — Context Buffer Management](https://claudefa.st/blog/guide/mechanics/context-buffer-management)

### 5.3 Our Context Assessment

Our CLAUDE.md (~400 lines) + skills + hooks + MEMORY.md index likely consume 15-20K tokens at session start. This leaves ~115-125K for actual work — adequate but could be tighter for long research sessions. Splitting CLAUDE.md would recover ~5-8K tokens.

---

## 6. Agent-to-Agent Communication Patterns

### 6.1 Five Coordination Patterns

| Pattern | Communication | Best For | Complexity |
|---------|--------------|----------|-----------|
| Generator-Verifier | Sequential feedback loop | Quality-critical output | Low |
| Orchestrator-Subagent | Hub-spoke delegation | Decomposed bounded tasks | Medium |
| Agent Teams | Peer-to-peer mailbox + shared task list | Parallel independent work | High |
| Message Bus | Pub/sub event topics | Event-driven pipelines | High |
| Shared State | Read/write to persistent store | Collaborative research synthesis | Medium |

Source: [Anthropic Blog — Multi-Agent Coordination](https://claude.com/blog/multi-agent-coordination-patterns)

### 6.2 Agent Teams (Claude Code Native, Feb 2026)

- **Flat peer model**: no central orchestrator, all agents are equal
- **Shared task list**: `~/.claude/tasks/{team-name}/` with status, ownership, dependencies
- **Mailbox system**: `SendMessage` tool for direct peer-to-peer or broadcast
- **File locking**: prevents concurrent edit conflicts

Best for: codebase migrations, parallel long-running subtasks. NOT for sequential dependency chains.

### 6.3 Subagents vs Agent Teams

| Feature | Subagents | Agent Teams |
|---------|-----------|-------------|
| Context | Shared parent context | Independent per agent |
| Communication | Report to parent only | Peer-to-peer mailbox |
| Lifecycle | Within single session | Persistent across tasks |
| Coordination | Parent as bottleneck | Shared task list |

### 6.4 Our Multi-Agent Assessment

We use subagents (researcher, etc.) in orchestrator-subagent pattern. This is appropriate for our workflow. Agent Teams would add value for parallel research investigations (e.g., simultaneously investigating LOH patterns across multiple chromosomes).

---

## 7. Actionable Recommendations for InterSubMod

### Priority 1 — CLAUDE.md Restructuring (Immediate)

- [ ] Split current 400+ line CLAUDE.md into <150 line core + `.claude/rules/` files
- [ ] Create path-scoped rules: `cpp-dev.md` (globs: src/, include/), `test-rules.md` (globs: tests/), `docs-rules.md` (globs: docs/)
- [ ] Add compaction directive to CLAUDE.md: preserve modified files, test commands, current task state

### Priority 2 — Hook Hardening (Short-term)

- [ ] Add pre-compact hook to capture session state before auto-compaction
- [ ] Add context-injection hook for directory-aware rule loading
- [ ] Evaluate StatusLine context-usage monitor (Rust or Python)

### Priority 3 — Skill Triggering (Short-term)

- [ ] Audit all skills: add "USE WHEN" triggers and file type mentions to descriptions
- [ ] Keep each skill under 50 lines; split research-loop into sub-skills if needed
- [ ] Consider hook-based skill router for deterministic activation

### Priority 4 — Memory Consolidation (Medium-term)

- [ ] Implement AutoDream-equivalent consolidation cycle (or adopt dream-skill)
- [ ] Add staleness tracking to MEMORY.md concluded items
- [ ] Evaluate claude-mem for semantic search across session history

### Priority 5 — Context Budget (Medium-term)

- [ ] Measure actual token usage at session start (system + CLAUDE.md + MEMORY.md + hooks)
- [ ] Set `CLAUDE_AUTOCOMPACT_PCT_OVERRIDE` to optimal value based on typical session length
- [ ] Evaluate 1M context window for long research sessions

---

## Sources

- [Anthropic — Best Practices](https://code.claude.com/docs/en/best-practices)
- [Anthropic — Memory](https://code.claude.com/docs/en/memory)
- [Anthropic — Hooks Guide](https://code.claude.com/docs/en/hooks-guide)
- [Anthropic — Agent Teams](https://code.claude.com/docs/en/agent-teams)
- [Anthropic — Multi-Agent Coordination Patterns](https://claude.com/blog/multi-agent-coordination-patterns)
- [DEV Community — CLAUDE.md From Basic to Adaptive](https://dev.to/cleverhoods/claudemd-best-practices-from-basic-to-adaptive-9lm)
- [DEV Community — 500 Lines of Rules](https://dev.to/mikeadolan/i-wrote-500-lines-of-rules-for-claude-code-heres-how-i-made-it-actually-follow-them-3c8)
- [DEV Community — Guaranteed Context Injection](https://dev.to/sasha_podles/claude-code-using-hooks-for-guaranteed-context-injection-2jg)
- [paddo.dev — Skills Auto-Activation](https://paddo.dev/blog/claude-skills-hooks-solution/)
- [paddo.dev — Skills Controllability Problem](https://paddo.dev/blog/claude-skills-controllability-problem/)
- [ClaudeFast — Context Buffer Management](https://claudefa.st/blog/guide/mechanics/context-buffer-management)
- [ClaudeFast — Compaction Explained](https://okhlopkov.com/claude-code-compaction-explained/)
- [Gist — Skill Structure Guide](https://gist.github.com/mellanon/50816550ecb5f3b239aa77eef7b8ed8d)
- [dream-skill (GitHub)](https://github.com/grandamenium/dream-skill)
- [claude-mem (GitHub)](https://github.com/thedotmack/claude-mem)
- [awesome-claude-code (GitHub)](https://github.com/hesreallyhim/awesome-claude-code)
- [Milvus Blog — memsearch](https://milvus.io/blog/adding-persistent-memory-to-claude-code-with-the-lightweight-memsearch-plugin.md)
