#!/bin/bash
echo "Starting Hyperspectral Field Viewer..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hyperspectral
streamlit run app.py
