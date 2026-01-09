#!/bin/bash
# =============================================================================
# Build InterSubMod inside Docker container
# =============================================================================
# Usage: ./docker/build.sh [debug|release]
# =============================================================================

set -e

BUILD_TYPE="${1:-Release}"

echo "=== InterSubMod Docker Build ==="
echo "Build Type: ${BUILD_TYPE}"
echo "================================"

# Create build directory
mkdir -p /workspace/build
cd /workspace/build

# Configure with CMake
echo "[1/3] Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"

# Build
echo "[2/3] Building..."
make -j$(nproc)

# Run tests
echo "[3/3] Running tests..."
ctest -R "^(run_tests|test_phase|Distance|Hierarchical|Significance|Global|Local|Structure|Bootstrap|BamReader|Config)" --output-on-failure

echo ""
echo "=== Build Complete ==="
echo "Executables: /workspace/build/bin/"
ls -la /workspace/build/bin/
