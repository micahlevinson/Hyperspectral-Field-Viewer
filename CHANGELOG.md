# Changelog

## [Unreleased]
- GeoTIFF export for VI maps (preserving geospatial coordinates)
- Shapefile zone statistics (pull per-zone VI stats)
- Multi-date change detection (load two flights, visualize difference)
- Automated field report PDF export

## [1.0.0] — 2026
### Added
- BSQ + HDR file loading via file path (large files) or upload (small files)
- 17 built-in vegetation indices across 5 categories:
  - **Greenness/Biomass:** NDVI, GNDVI, EVI, SAVI, OSAVI, VARI
  - **Chlorophyll/Nitrogen:** NDRE, CIre, MCARI, REIP, RENDVI
  - **Water/Stress:** NDWI, WBI, PSRI
  - **Disease/Pigment:** PRI, SIPI
  - **LAI/Structure:** MTVI2
- RGB composite viewer with preset band combinations (true color, CIR, red edge)
- Custom VI formula builder using B\<wavelength\> syntax
- Swipe curtain comparison between any two computed views, with colorbar legends
- Polygon ROI spectrum inspector — click vertices on a Plotly image to define zones
- Single-pixel spectrum inspector with up to 5 pixels overlaid
- Field-wide VI statistics with violin/distribution plots and CSV export
- Band-by-band statistics with mean reflectance spectrum plot
- Sensor-aware VI availability — indices outside the sensor's wavelength range are disabled
- Conda environment file and one-command launch scripts for Windows and Mac/Linux
