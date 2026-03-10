#!/bin/bash

# Script to run MML_ParticleVisualizer2D_FLTK with all test data files

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VISUALIZER="$SCRIPT_DIR/MML_ParticleVisualizer2D_FLTK"
DATA_DIR="$SCRIPT_DIR/../../../data/ParticleVisualizer2D/Basic"

# Check if visualizer executable exists
if [ ! -f "$VISUALIZER" ]; then
    echo "Error: Visualizer not found at $VISUALIZER"
    echo "Please build the project first."
    exit 1
fi

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory not found at $DATA_DIR"
    exit 1
fi

echo "Running MML_ParticleVisualizer2D with test data files..."
echo "Press Ctrl+C to stop"
echo ""

# Loop through all .txt files in the data directory
for datafile in "$DATA_DIR"/*.txt; do
    if [ -f "$datafile" ]; then
        filename=$(basename "$datafile")
        echo "======================================"
        echo "Loading: $filename"
        echo "======================================"
        
        # Run the visualizer with the data file
        "$VISUALIZER" "$datafile"
        
        # Wait a moment between runs
        sleep 0.5
        echo ""
    fi
done

echo "All test files processed."
