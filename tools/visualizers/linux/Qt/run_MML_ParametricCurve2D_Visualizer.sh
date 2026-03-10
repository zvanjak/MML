#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LD_LIBRARY_PATH="$SCRIPT_DIR/lib:$LD_LIBRARY_PATH"
export QT_QPA_PLATFORM_PLUGIN_PATH="$SCRIPT_DIR/plugins/platforms"
export QT_PLUGIN_PATH="$SCRIPT_DIR/plugins"
exec "$SCRIPT_DIR/MML_ParametricCurve2D_Visualizer_Qt" "$@"
