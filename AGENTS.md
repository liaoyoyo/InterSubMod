# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the C++ core (split into `core/`, `io/`, `utils/`); headers live in `include/`.
- `tests/` contains GoogleTest unit tests; `src/test/` contains phase-specific test drivers.
- `tools/` houses Python analysis/plotting utilities; `scripts/` contains shell workflows.
- `data/` stores example inputs; `output/` is generated analysis output; `docs/` and `images/` support documentation.
- `build/` is the out-of-tree build output created by CMake.

## Build, Test, and Development Commands
- Build: `mkdir -p build && cd build && cmake .. && make -j$(nproc)`; binary at `build/bin/inter_sub_mod`.
- Run core manually: `./build/bin/inter_sub_mod --tumor-bam data/tumor.bam --reference data/ref.fa --vcf data/somatic.vcf --output-dir results`.
- Full pipeline script: `./scripts/run_vcf_all_snv.sh --mode all-with-w1000 --plot-type distance` (see `--help` for options).
- Output checks: `./scripts/verify_output.sh` validates expected files and matrix dimensions.
- Python deps for plotting: `pip install -r requirements.txt`.
- Optional container: `docker build -f Dockerfile.dev -t intersubmod:dev .` and `docker run -it --rm -v $(pwd):/workspace intersubmod:dev`.

## Coding Style & Naming Conventions
- C++17 code with `.hpp` headers and `.cpp` sources; namespace is `InterSubMod`.
- Formatting follows `.clang-format` (Google base, 4-space indent, 120 column limit); run `clang-format` on touched C++ files.
- Naming patterns: `CamelCase` classes (e.g., `BamReader`), `snake_case` methods and files.

## Testing Guidelines
- Unit tests live in `tests/test_*.cpp`; run with `ctest --test-dir build` or `./build/bin/run_tests`.
- Phase tests compile to `build/bin/test_phase*` from `src/test/`; `scripts/run_random_snv_test.sh` provides a quick smoke test.
- No explicit coverage target is enforced; add GTest coverage for new core logic when feasible.

## Commit & Pull Request Guidelines
- Commit messages in recent history use short imperative summaries like `Add ...` or `Refactor: ...`; keep to one line and add a prefix (`Fix:`, `Docs:`) when helpful.
- PRs should include a concise summary, commands run, and sample outputs/logs or plots when analysis or visualization changes.

## Data, Outputs, and Configuration Tips
- Many scripts default to absolute `/big8_disk/...` paths; prefer overriding via flags like `--vcf`, `--out`, `--threads`, and `--plot-type`.
- Keep generated artifacts in `output/` and avoid committing large datasets unless explicitly requested.
