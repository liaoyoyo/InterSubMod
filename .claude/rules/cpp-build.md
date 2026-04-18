---
globs:
  - "src/**/*.cpp"
  - "src/**/*.hpp"
  - "include/**/*.hpp"
  - "include/**/*.h"
  - "CMakeLists.txt"
  - "tests/**/*.cpp"
---

# C++ Build & Change Checklist

## Build Commands

```bash
# Configure and build (from project root)
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Debug build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# Run all tests
cd build && ctest --output-on-failure

# Run specific test
./build/tests/test_<name>
```

## Post-Change Checklist

After modifying C++ code, complete these steps:

```bash
# 1. Compile
cd build && make -j$(nproc)

# 2. Run tests
./scripts/run_batch_vcf_analysis.sh

# 3. Verify results are reasonable

# 4. Update Docker config (if needed)
```

## Adding a New Module

1. Create header in `include/core/`
2. Create implementation in `src/core/`
3. Update `CMakeLists.txt` (if needed)
4. Write unit tests in `tests/`

## Code Standards

- **C++17** standard
- Follow `.clang-format` (Google style, 120 char line width, 4 space indent)
- Format before commit: `clang-format -i <file>`
- **Code comments**: English
