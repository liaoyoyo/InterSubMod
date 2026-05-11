# Case Study: /tmp 800 GB Disaster (2026-05-08)

## What happened

2026-05-07/08 — InterSubMod ran longphase-to-mod haplotag on full-genome BAM. Default output paths: `/tmp/t12_genome_{baseline,v3f,v5}/tumor_tagged.bam`. Three files = 268 + 268 + 174 = **710 GB** written to /tmp.

This server's /tmp is **NOT tmpfs** — it's a regular directory on root volume `/dev/mapper/ubuntu--vg-ubuntu--lv` (1 TB). Root partition hit 100% full.

## Cascade failure

- Claude Code Bash tool fully broken (`/tmp/claude-${UID}/` mkdir ENOSPC, even `true`/`echo` returns exit 1)
- All 33 server accounts: cron jobs / docker / apt / journald simultaneously crashed
- Only diagnosis path: SSH in, run `df -h /` to see 100% root partition

## Root cause

Default tooling assumption: `/tmp` is RAM-based tmpfs (Linux tutorials usually assume this). On THIS server, /tmp is regular root-volume directory. Assumption mismatch → catastrophic outcome.

## 5 prevention rules (apply to EVERY pipeline)

### 1. Any pipeline writing >100 MB intermediate must explicit redirect

| Tool | Default | Override |
|------|---------|----------|
| longphase / haplotag | `/tmp/tumor_tagged.bam` | `-o /big7_disk/.../<output>` |
| samtools sort/view --threads | `$TMPDIR or /tmp` | `export TMPDIR=/big7_disk/liaoyoyo2001/tmp` |
| bcftools / minimap2 / star | `$TMPDIR or /tmp` | same |
| Python `tempfile.NamedTemporaryFile()` | `/tmp` | same |
| pandas read_csv compression decompress | `/tmp` | same |

### 2. Wrapper script standard preamble

```bash
#!/bin/bash
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
mkdir -p "$TMPDIR"
set -euo pipefail
```

### 3. Measure before full run

Run chr19/chr1 pilot first, measure single-chromosome intermediate size, multiply by 24 for upper bound. If upper bound > `df --output=avail /tmp` minus 50 GB margin, **abort immediately**.

### 4. Mount path safety table

| Path | Verdict | Why |
|------|---------|-----|
| `/tmp` `/var/tmp` `/root` `/var/log` | NO | root volume, 1 TB shared with all 33 accounts |
| `/big7_disk/liaoyoyo2001/` | YES | user disk, TB scale |
| `/big8_disk/liaoyoyo2001/` | YES | another user disk |
| `/bip7_disk/liaoyoyo2001/tmp_claude/` | YES | Claude Code symlink rescue dir |
| `/gpu_disk*` | CAUTION | NFS, slow; avoid heavy random I/O |

### 5. Bash silent failure → first check

```bash
df -h /     # root partition
df -h /tmp  # is /tmp independently mounted?
mount | grep tmp  # confirm /tmp filesystem type
```

Exit 1 + zero stdout from any command = 99% probability of disk full.

## Detection in retrospective

For first run on a new server, MUST run `mount | grep tmp` to confirm /tmp filesystem type before launching any GB-scale pipeline.
