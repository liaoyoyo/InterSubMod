# Debug Notes - SNV Loading Module

## 2025-11-22: CMake DEBUG Macro Conflict

### Issue

Compilation failed with errors like `expected identifier before numeric constant` in `Logger.hpp` for the enum `LogLevel::DEBUG`.

### Cause

`CMakeLists.txt` defines `-DDEBUG` when `CMAKE_BUILD_TYPE` is `Debug`:

```cmake
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0)
    add_definitions(-DDEBUG)
```

This preprocessor definition replaced the token `DEBUG` in the C++ code with `1`, breaking the enum definition.

### Solution

Renamed `LogLevel` enum values to avoid collision with common macros:

- `DEBUG` -> `L_DEBUG`
- `INFO` -> `L_INFO`
- `WARNING` -> `L_WARNING`
- `ERROR` -> `L_ERROR`

This is a safer practice in C++ projects that might interact with legacy C headers or system definitions.
