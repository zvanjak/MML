#!/usr/bin/env bash
# build.sh - Build MML_VisualizationApp
# Usage:
#   ./build.sh            # Release build (default)
#   ./build.sh Debug      # Debug build

set -euo pipefail

CONFIG="${1:-Release}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "Building MML_VisualizationApp ($CONFIG) ..."

cmake --build "$ROOT/build" \
      --config "$CONFIG" \
      --target MML_VisualizationApp \
      --parallel

echo "Build succeeded."
echo "  Executable: $ROOT/build/src/visualization_examples/MML_VisualizationApp"
