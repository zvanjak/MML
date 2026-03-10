#!/bin/bash

# Script to visualize Particle Visualizer 2D data files
# Using Qt MML_ParticleVisualizer2D

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VISUALIZER="$SCRIPT_DIR/MML_ParticleVisualizer2D"
DATA_DIR="$SCRIPT_DIR/../../../data/ParticleVisualizer2D/Basic"

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================"
echo -e "Particle Visualizer 2D Runner (Qt)"
echo -e "========================================${NC}"
echo ""

# Check if visualizer exists
if [ ! -f "$VISUALIZER" ]; then
    echo -e "${RED}ERROR: Visualizer not found at: $VISUALIZER${NC}"
    echo -e "${YELLOW}Please ensure the Qt visualizer is built and deployed.${NC}"
    exit 1
fi

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo -e "${RED}ERROR: Data directory not found at: $DATA_DIR${NC}"
    exit 1
fi

# Get all .txt files
DATA_FILES=("$DATA_DIR"/*.txt)

if [ ${#DATA_FILES[@]} -eq 0 ] || [ ! -f "${DATA_FILES[0]}" ]; then
    echo -e "${YELLOW}No data files found in: $DATA_DIR${NC}"
    exit 0
fi

echo -e "${CYAN}Found ${#DATA_FILES[@]} data files in:${NC}"
echo -e "  $DATA_DIR"
echo ""
echo -e "${CYAN}Visualizer:${NC}"
echo -e "  $VISUALIZER"
echo ""
echo -e "${GREEN}========================================${NC}"
echo ""

# Counter
FILE_NUM=0

# Process each file
for FILE in "${DATA_FILES[@]}"; do
    ((FILE_NUM++))
    FILENAME=$(basename "$FILE")
    
    echo -e "${YELLOW}[$FILE_NUM/${#DATA_FILES[@]}] Processing: $FILENAME${NC}"
    
    "$VISUALIZER" "$FILE"
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -ne 0 ]; then
        echo -e "  ${RED}WARNING: Visualizer exited with code $EXIT_CODE${NC}"
    else
        echo -e "  ${GREEN}Completed successfully${NC}"
    fi
    
    echo ""
done

echo -e "${GREEN}========================================"
echo -e "All files processed: $FILE_NUM"
echo -e "========================================${NC}"
