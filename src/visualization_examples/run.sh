#!/usr/bin/env bash
# run.sh - Run MML_VisualizationApp
# Usage:
#   ./run.sh                         # default demo
#   ./run.sh forms_3d                # typed forms single-vector scene
#   ./run.sh forms_field_3d          # vortex tube field scene
#   ./run.sh forms_em_3d             # EM dipole field scene
#   ./run.sh real_function           # y = f(x) example
#   ./run.sh field_3d                # 3D vector field
#   ./run.sh all                     # run all visualizations

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# On Linux/macOS single-config generators (Make/Ninja) the exe is directly
# under the build output directory, with no Debug/Release subfolder.
EXE="$ROOT/build/src/visualization_examples/MML_VisualizationApp"

if [ ! -f "$EXE" ]; then
    echo "Executable not found: $EXE"
    echo "Run ./build.sh first."
    exit 1
fi

if [ "${1:-}" = "" ]; then
    echo "Running MML_VisualizationApp (default demo) ..."
    "$EXE"
else
    echo "Running MML_VisualizationApp $1 ..."
    "$EXE" "$1"
fi
