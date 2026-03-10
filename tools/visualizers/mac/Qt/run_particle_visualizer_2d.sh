#!/bin/bash

# Script to run MML_ParticleVisualizer2D_Qt with all test data files

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VISUALIZER="$SCRIPT_DIR/MML_ParticleVisualizer2D_Qt"
DATA_DIR="$SCRIPT_DIR/../../../data/ParticleVisualizer2D"

if [ ! -f "$VISUALIZER" ]; then
    echo "Error: Visualizer not found at $VISUALIZER"
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory not found at $DATA_DIR"
    exit 1
fi

echo "Running MML_ParticleVisualizer2D_Qt with test data files..."
echo "Press Ctrl+C to stop"
echo ""

for datafile in "$DATA_DIR"/*.txt; do
    if [ -f "$datafile" ]; then
        filename=$(basename "$datafile")
        echo "======================================"
        echo "Loading: $filename"
        echo "======================================"
        "$VISUALIZER" "$datafile"
        echo ""
        echo "Press Enter to continue to next file, or Ctrl+C to exit..."
        read
    fi
done

echo "All test files completed!"
