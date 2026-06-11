#!/bin/bash

# Download a Figshare file
PROJECT_ID="52692896"
TARGET_FILE="GTDB_sliced_seqs_sliding_window.fna.gz"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if script exists
if [ ! -f "$SCRIPT_DIR/api.py" ]; then
    echo "ERROR: api.py not found in $SCRIPT_DIR"
    exit 1
fi

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python 3 is required but not installed"
    exit 1
fi

# Check if requests library is installed
python3 -c "import requests" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Python package: requests is required but not installed"
    exit 1
fi

echo "All required tools are available: downloading file..."
python3 "$SCRIPT_DIR/api.py" "$PROJECT_ID" "$TARGET_FILE"

if [ $? -eq 0 ]; then
    echo "Download successful!"
else
    echo "Download failed!"
    exit 1
fi
