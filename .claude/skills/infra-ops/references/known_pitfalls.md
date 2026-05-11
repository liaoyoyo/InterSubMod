# Known Infrastructure Pitfalls

## 1. /tmp cleanup (root >85% AND /tmp <1GB free)

**Trigger**: `df -h /` shows root partition >85%, AND `df -h /tmp` shows <1 GB free.
**Investigation**: `du -sh /tmp/* | sort -h | tail -20` to find biggest occupants.
**Fix**: First identify which process owns the file (`lsof /tmp/<biggest>`). If owned by killed/orphaned process, safe to remove. If active, wait for process to finish or kill explicitly.
**Pattern**: dry-run list → AskUserQuestion confirm → execute → audit log.
**Reference**: `tmp_disaster_2026_05_08.md`.

## 2. CMake cache reset

**Trigger**: `make` fails with mysterious "no rule to make target", `CMakeCache.txt` mtime older than `CMakeLists.txt`.
**Fix**: `rm -rf build && mkdir build && cd build && cmake ..` — full reset. Lower-impact alternative: `cmake --fresh ..` (CMake 3.24+).
**Caveat**: re-cmake erases `compile_commands.json` — IDE may need re-config.

## 3. Output archive (synthesis/ >500 GB)

**Trigger**: `output/synthesis/` exceeds 500 GB; new runs need space.
**Fix**: identify oldest concluded research rounds (per `docs/data_specs/output 信任度` spec), `mv` to `output/bip8_output_archive/<YYYYMM>/` first, then optionally `gzip` archived TSVs. NEVER `rm` without archive copy. Update `output/OUTPUT_INDEX.md`.

## 4. Conda env reset (broken activation)

**Trigger**: `conda activate <env>` fails or env reports missing core packages.
**Investigation**: `conda info --envs`, check for `.conda` lockfiles or partially-installed states.
**Fix**: `conda deactivate && conda env list` to confirm corruption, then `conda env remove -n <env>` and recreate from `environment.yml` if available.

## 5. Docker dangling images (>50 GB)

**Trigger**: `docker system df` shows >50 GB in dangling/unused images.
**Fix**: dry-run with `docker system df -v`, then `docker system prune -a --filter "until=720h"` (30+ day old images only). Be careful: `docker system prune` is destructive — confirm no critical images before.

## 6. Log rotation (/var/log >50 GB)

**Trigger**: `du -sh /var/log` exceeds 50 GB.
**Investigation**: `du -sh /var/log/* | sort -h | tail -10`.
**Fix**: usually requires sudo. Common offender: journald with `SystemMaxUse=` not set. Set in `/etc/systemd/journald.conf`, then `systemctl restart systemd-journald`. Old journal: `journalctl --vacuum-size=2G`.
