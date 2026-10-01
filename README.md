# 🌱 Hyperspectral Field Viewer

A browser-based tool for exploring UAV hyperspectral imagery in agricultural settings.
Load ENVI BSQ/HDR files, visualize vegetation indices, inspect pixel and polygon spectra,
and pull field-wide statistics — no GIS software or coding required.

Built by the Micah Levinson as part of the [UC Davis Crop Sensing Group](https://cropscience.ucdavis.edu) (Department of Viticulture and Enology).

---

## Features

| Feature | Description |
|---|---|
| 🖼️ **RGB Composite** | True color, color infrared, red edge, or custom band assignments |
| 📊 **VI Browser** | 17 built-in vegetation indices across 5 categories |
| 🔀 **Swipe Compare** | Drag a curtain slider between any two computed views |
| ✏️ **Custom VI** | Write your own index formula using `B<wavelength>` syntax |
| 📈 **Spectra** | Click individual pixels or draw polygon ROIs to extract spectra |
| 📉 **Statistics** | Field-wide stats, histograms, violin plots, CSV export |
| 📂 **Large file support** | Load files directly from disk via path — no upload size limit |

---

## Quick Start

### Requirements
- Python 3.10 or newer — download at [python.org](https://www.python.org/downloads/)
- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (recommended)
- 8 GB RAM minimum · 16 GB recommended for large files

---

### Option A — Conda (recommended)

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/HyperspectralFieldViewer.git
cd HyperspectralFieldViewer

# 2. Create and activate the environment (one-time setup)
conda env create -f environment.yml
conda activate hyperspectral

# 3. Launch the app
streamlit run app.py
```

A browser tab opens automatically at **http://localhost:8501**

---

### Option B — pip

```bash
git clone https://github.com/<your-username>/HyperspectralFieldViewer.git
cd HyperspectralFieldViewer
pip install -r requirements.txt
streamlit run app.py
```

---

### Option C — Double-click launcher

After completing the conda setup above, you can launch the app anytime by
double-clicking the launcher script in the project folder:

- **Windows:** `launch_windows.bat`
- **Mac / Linux:** `launch_mac_linux.sh`

---

## Loading Your Data

The app reads **ENVI BSQ** format files — a binary image cube (`.bsq` or `.img`)
paired with a plain-text metadata header (`.hdr`). Both files must have the same
base name and live in the same folder:

```
field_july.hdr    ← metadata (wavelengths, dimensions, map info)
field_july.bsq    ← binary image data
```

**For large files (> 200 MB):**
1. Select **📁 File path** in the sidebar
2. Paste the full path to your `.hdr` file (e.g. `D:\data\field_july.hdr`)
3. Click **Load File** — the `.bsq` is found automatically

**For small files:**
1. Select **⬆️ Upload** in the sidebar
2. Drag in both files

> ⚠️ Data files (`.bsq`, `.hdr`, `.tif`, etc.) are excluded from this repository
> by `.gitignore`. Do not commit raw hyperspectral imagery to GitHub.

---

## Tabs

### 📋 Overview
File metadata summary, band list, quick-look RGB composite, and a catalog of
which vegetation indices are computable with your sensor's wavelength range.

### 📊 VI Browser
Select any vegetation index or RGB composite preset, adjust the color scale, and
click **Compute & View**. Results are cached — switching between already-computed
indices is instant. Quick stats (mean, median, std dev) and a histogram are shown
below each map.

### 🔀 Compare
Compute at least two views in the VI Browser first, then come here to place them
side by side. Drag the white divider left and right to reveal one or the other.
Each VI shows a colorbar legend with its value range pinned to its corner.

### ✏️ Custom VI
Write your own index using `B<wavelength>` references:

```python
# Examples
(B800 - B670) / (B800 + B670 + eps)       # NDVI variant
(B750 - B705) / (B750 + B705 + eps)       # red-edge ratio
np.log(B800 / (B670 + eps))               # log ratio
(B531 - B570) / (B531 + B570 + eps)       # PRI variant
```

`eps` (1×10⁻⁶) is pre-defined to prevent division by zero.
`np` (NumPy) is available for math functions.

### 📈 Spectra
**Click pixel mode:** Display a reference image, hover to read coordinates, then
enter a row/column and plot the full reflectance spectrum. Up to 5 pixels can be
overlaid for comparison.

**Polygon ROI mode:** Click directly on the image to place polygon vertices.
Save named ROIs and extract the mean ± 1 std dev spectrum for each zone —
useful for comparing a healthy tree crown vs. a stressed one.

### 📉 Statistics
Compute field-wide statistics for one or more VIs simultaneously:
min, max, mean, median, std dev, and percentiles. Download as CSV.
Also includes band-by-band reflectance statistics with a mean spectrum plot.

---

## Vegetation Index Reference

### Greenness / Biomass

| Index | Formula | Wavelengths | Best for |
|---|---|---|---|
| **NDVI** | (NIR − Red) / (NIR + Red) | 800, 670 nm | Overall crop health and density |
| **GNDVI** | (NIR − Green) / (NIR + Green) | 800, 550 nm | Chlorophyll; less saturation than NDVI |
| **EVI** | 2.5×(NIR−Red) / (NIR+6×Red−7.5×Blue+1) | 800, 670, 480 nm | Dense canopy; atmospheric correction |
| **SAVI** | 1.5×(NIR−Red) / (NIR+Red+0.5) | 800, 670 nm | Early season; sparse canopy / bare soil |
| **OSAVI** | (NIR−Red) / (NIR+Red+0.16) | 800, 670 nm | Broad-use SAVI with optimized soil factor |
| **VARI** | (Green−Red) / (Green+Red−Blue) | 550, 670, 480 nm | Visible-only greenness; atmospheric robust |

### Chlorophyll / Nitrogen

| Index | Formula | Wavelengths | Best for |
|---|---|---|---|
| **NDRE** | (NIR − RedEdge) / (NIR + RedEdge) | 790, 720 nm | Nitrogen status; more sensitive than NDVI in dense crops |
| **CIre** | (NIR / RedEdge) − 1 | 800, 720 nm | Chlorophyll content; N use efficiency |
| **MCARI** | [(RE−Red)−0.2×(RE−Green)]×(RE/Red) | 700, 670, 550 nm | Chlorophyll absorption; early stress |
| **REIP** | 700+40×[(R670+R780)/2−R700]/(R740−R700) | 670, 700, 740, 780 nm | Red-edge position; shifts with chlorophyll |
| **RENDVI** | (R750−R705) / (R750+R705) | 750, 705 nm | Canopy health; red-edge combined with structure |

### Water / Stress

| Index | Formula | Wavelengths | Best for |
|---|---|---|---|
| **NDWI** | (Green − NIR) / (Green + NIR) | 860, 1240 nm | Plant water content (requires SWIR) |
| **WBI** | R900 / R970 | 900, 970 nm | Canopy water stress; values < 1 indicate stress |
| **PSRI** | (Red − Green) / RedEdge | 678, 500, 750 nm | Senescence and carotenoid:chlorophyll ratio |

### Disease / Pigment

| Index | Formula | Wavelengths | Best for |
|---|---|---|---|
| **PRI** | (R531 − R570) / (R531 + R570) | 531, 570 nm | Light use efficiency; early drought detection |
| **SIPI** | (NIR − Blue) / (NIR − Red) | 800, 445, 680 nm | Carotenoid:chlorophyll; disease and senescence |

### LAI / Structure

| Index | Formula | Wavelengths | Best for |
|---|---|---|---|
| **MTVI2** | 1.5×[1.2×(NIR−Green)−2.5×(Red−Green)] / √[...] | 800, 670, 550 nm | Leaf area index; soil and chlorophyll resistant |

> **Note:** Indices requiring wavelengths outside your sensor's range are automatically
> disabled in the VI Browser.

---

## RGB Composite Presets

| Name | Red | Green | Blue | Best for |
|---|---|---|---|---|
| True color | ~660 nm | ~550 nm | ~480 nm | General visual overview |
| Color infrared (CIR) | ~800 nm | ~660 nm | ~550 nm | Vegetation health at a glance |
| Red edge | ~750 nm | ~700 nm | ~670 nm | Crop stress and nitrogen status |

All three channels are adjustable — override any preset with any band in your dataset.

---

## Tips for Large Files

- First load reads the full array into memory — expect 30–60 seconds for 800 MB files
- After the first load, switching between VIs is fast (array stays cached)
- Computed VI maps are also cached — re-rendering is instant
- Close other memory-heavy applications before loading very large files
- On a machine with 16 GB RAM, files up to ~800 MB load reliably

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'streamlit'`**
→ You're in the wrong environment. Run `conda activate hyperspectral` before `streamlit run app.py`.

**`FileNotFoundError` when loading via path**
→ Check that your `.bsq` / `.img` file has the exact same base name as the `.hdr` and is in the same folder.

**Blank metric boxes in Overview tab**
→ No file is loaded yet. Load your data first using the sidebar.

**VI is greyed out / unavailable**
→ That index requires wavelengths outside your sensor's range. Check the Overview tab for a full list of available and unavailable indices.

**App is slow on first VI compute**
→ Normal — the array is being read from disk. After that, VIs compute in seconds.

**Polygon ROI extracts no pixels**
→ The polygon may be too small at image scale. Try drawing a larger region, or zoom in on the reference image first to place vertices more precisely.

---

## Dependencies

| Package | Purpose |
|---|---|
| `streamlit` | Browser-based UI framework |
| `spectral` (SPy) | ENVI BSQ/HDR file reading |
| `numpy` | Array math and band calculations |
| `pandas` | Statistics tables and CSV export |
| `plotly` | Interactive maps, spectra, histograms |
| `matplotlib` | Colorbar rendering for swipe view |
| `Pillow` | Image resizing and PNG encoding |
| `scipy` | Scientific computing utilities |
| `rasterio` | Geospatial raster support |

---

## Roadmap

- [ ] GeoTIFF export for VI maps (preserving geospatial coordinates)
- [ ] Shapefile zone statistics (per-zone VI stats from a field boundary file)
- [ ] Multi-date change detection (load two flights, compute difference maps)
- [ ] Automated field report PDF export

---

## Contributing

Issues and pull requests are welcome. If you encounter a bug or want to request
a feature, please open an issue on GitHub with:
- Your operating system and Python version
- The error message or unexpected behavior
- Steps to reproduce (you don't need to share your data files)

---

## Citation / Acknowledgments

If you use this tool in published research, please acknowledge Micah Levinsona and the UC Davis Crop Sensing Group, Department of Viticulture and Enology.

---

## License

MIT — see [LICENSE](LICENSE) for details.
