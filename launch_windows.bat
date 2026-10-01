@echo off
echo Starting Hyperspectral Field Viewer...
call conda activate hyperspectral
streamlit run app.py
pause
