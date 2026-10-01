"""
Hyperspectral Field Viewer
==========================
Load BSQ + HDR files, explore vegetation indices, pull statistics,
inspect pixel spectra, and add custom VIs — all in a browser.

Run with:
    streamlit run app.py
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import plotly.graph_objects as go
import plotly.express as px
import spectral
import spectral.io.envi as envi
import base64, io, os, tempfile, warnings
from pathlib import Path
from PIL import Image
warnings.filterwarnings("ignore")

# ─── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Hyperspectral Field Viewer",
    page_icon="🌱",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .main .block-container { padding-top: 1.5rem; }
    .vi-card { background: #f0f7ed; border-left: 4px solid #2d7a3a; padding: 0.6rem 0.9rem;
               border-radius: 0 8px 8px 0; margin-bottom: 0.4rem; }
    h1 { color: #1a4f24; }
    .stTabs [data-baseweb="tab"] { font-size: 1rem; }
</style>
""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
#  VI LIBRARY
# ══════════════════════════════════════════════════════════════════════════════
VI_LIBRARY = {
    "NDVI":   {"full_name":"Normalized Difference Vegetation Index","formula_display":"(NIR − Red) / (NIR + Red)","description":"Classic measure of green biomass and photosynthetic activity. Saturates in dense canopies.","category":"Greenness / Biomass","bands":{"NIR":800,"Red":670},"range":(-1,1),"colormap":"RdYlGn"},
    "GNDVI":  {"full_name":"Green NDVI","formula_display":"(NIR − Green) / (NIR + Green)","description":"More sensitive to chlorophyll than NDVI; less saturation in dense crops.","category":"Greenness / Biomass","bands":{"NIR":800,"Green":550},"range":(-1,1),"colormap":"YlGn"},
    "EVI":    {"full_name":"Enhanced Vegetation Index","formula_display":"2.5×(NIR−Red)/(NIR+6×Red−7.5×Blue+1)","description":"Reduces atmospheric and soil background effects. Better in dense canopies.","category":"Greenness / Biomass","bands":{"NIR":800,"Red":670,"Blue":480},"range":(-1,1),"colormap":"YlGn"},
    "SAVI":   {"full_name":"Soil-Adjusted Vegetation Index","formula_display":"1.5×(NIR−Red)/(NIR+Red+0.5)","description":"Minimizes soil brightness. Best early season when cover is sparse.","category":"Greenness / Biomass","bands":{"NIR":800,"Red":670},"range":(-1,1),"colormap":"YlGn"},
    "OSAVI":  {"full_name":"Optimized SAVI","formula_display":"(NIR−Red)/(NIR+Red+0.16)","description":"Improved SAVI with optimized soil factor. Robust across canopy densities.","category":"Greenness / Biomass","bands":{"NIR":800,"Red":670},"range":(-1,1),"colormap":"YlGn"},
    "VARI":   {"full_name":"Visible Atmospherically Resistant Index","formula_display":"(Green−Red)/(Green+Red−Blue)","description":"Green vegetation fraction using visible bands only. Robust to atmosphere.","category":"Greenness / Biomass","bands":{"Green":550,"Red":670,"Blue":480},"range":(-1,1),"colormap":"YlGn"},
    "NDRE":   {"full_name":"Normalized Difference Red Edge","formula_display":"(NIR − RedEdge) / (NIR + RedEdge)","description":"Sensitive to chlorophyll and nitrogen. More accurate than NDVI in dense canopies.","category":"Chlorophyll / Nitrogen","bands":{"NIR":790,"RedEdge":720},"range":(-1,1),"colormap":"RdYlGn"},
    "CIre":   {"full_name":"Chlorophyll Index Red Edge","formula_display":"(NIR / RedEdge) − 1","description":"Strong predictor of chlorophyll content and nitrogen use efficiency.","category":"Chlorophyll / Nitrogen","bands":{"NIR":800,"RedEdge":720},"range":(0,10),"colormap":"YlGn"},
    "MCARI":  {"full_name":"Modified Chlorophyll Absorption Ratio Index","formula_display":"[(RE−Red)−0.2×(RE−Green)]×(RE/Red)","description":"Tracks chlorophyll absorption relative to leaf reflectance. Sensitive to early stress.","category":"Chlorophyll / Nitrogen","bands":{"RedEdge":700,"Red":670,"Green":550},"range":(0,1),"colormap":"RdYlGn"},
    "REIP":   {"full_name":"Red Edge Inflection Point","formula_display":"700+40×[(R670+R780)/2−R700]/(R740−R700)","description":"Position of the red-edge inflection. Higher = more chlorophyll = healthier.","category":"Chlorophyll / Nitrogen","bands":{"R670":670,"R700":700,"R740":740,"R780":780},"range":(700,740),"colormap":"RdYlGn"},
    "RENDVI": {"full_name":"Red Edge NDVI","formula_display":"(R750−R705)/(R750+R705)","description":"Combines chlorophyll sensitivity of red edge with normalized structure.","category":"Chlorophyll / Nitrogen","bands":{"R750":750,"R705":705},"range":(-1,1),"colormap":"RdYlGn"},
    "NDWI":   {"full_name":"Normalized Difference Water Index","formula_display":"(Green−NIR)/(Green+NIR)","description":"Sensitive to plant water content. Higher = more moisture.","category":"Water / Stress","bands":{"Green":860,"NIR":1240},"range":(-1,1),"colormap":"RdBu"},
    "WBI":    {"full_name":"Water Band Index","formula_display":"R900/R970","description":"Direct measure of canopy water. <1 indicates stress.","category":"Water / Stress","bands":{"R900":900,"R970":970},"range":(0.8,1.2),"colormap":"RdBu"},
    "PSRI":   {"full_name":"Plant Senescence Reflectance Index","formula_display":"(Red−Green)/RedEdge","description":"Carotenoid:chlorophyll ratio. Rises with senescence or stress.","category":"Water / Stress","bands":{"Red":678,"Green":500,"RedEdge":750},"range":(-0.2,0.2),"colormap":"RdYlGn_r"},
    "PRI":    {"full_name":"Photochemical Reflectance Index","formula_display":"(R531−R570)/(R531+R570)","description":"Xanthophyll cycle activity — proxy for light use efficiency and drought stress.","category":"Disease / Pigment","bands":{"R531":531,"R570":570},"range":(-0.2,0.2),"colormap":"RdYlGn"},
    "SIPI":   {"full_name":"Structure Insensitive Pigment Index","formula_display":"(NIR−Blue)/(NIR−Red)","description":"Carotenoids to chlorophyll. Increases with stress, senescence, or disease.","category":"Disease / Pigment","bands":{"NIR":800,"Blue":445,"Red":680},"range":(0.5,1.8),"colormap":"RdYlGn_r"},
    "MTVI2":  {"full_name":"Modified Triangular Vegetation Index 2","formula_display":"1.5×[1.2×(NIR−Green)−2.5×(Red−Green)]/√[...]","description":"Sensitive to LAI, resistant to soil and chlorophyll variation.","category":"LAI / Structure","bands":{"NIR":800,"Red":670,"Green":550},"range":(0,1),"colormap":"YlGn"},
}

# ══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ══════════════════════════════════════════════════════════════════════════════
@st.cache_data
def find_band(wavelengths, target_nm):
    return int(np.argmin(np.abs(np.array(wavelengths, dtype=float) - target_nm)))

def get_wavelengths(img):
    meta = img.metadata
    if "wavelength" in meta:
        return [float(w) for w in meta["wavelength"]]
    return list(range(1, img.nbands + 1))

def get_meta_summary(img):
    m = img.metadata
    wl = get_wavelengths(img)
    return {
        "Rows": img.nrows, "Columns": img.ncols, "Bands": img.nbands,
        "Wavelength min (nm)": f"{min(wl):.1f}" if wl else "—",
        "Wavelength max (nm)": f"{max(wl):.1f}" if wl else "—",
        "Data type": m.get("data type","—"),
        "Interleave": str(m.get("interleave","—")).upper(),
        "Map info": m.get("map info","—"),
        "Sensor type": m.get("sensor type","—"),
        "Description": m.get("description","—"),
    }

@st.cache_data(show_spinner="Reading data into memory…")
def read_as_array(_img):
    return _img.load().astype(np.float32)

def compute_vi(arr, vi_key, wavelengths):
    vi = VI_LIBRARY[vi_key]; b = vi["bands"]; eps = 1e-6
    def band(nm):
        return arr[:,:,find_band(wavelengths,nm)].astype(np.float64)
    if vi_key=="NDVI":   nir,red=band(b["NIR"]),band(b["Red"]); r=(nir-red)/(nir+red+eps)
    elif vi_key=="GNDVI":nir,g=band(b["NIR"]),band(b["Green"]); r=(nir-g)/(nir+g+eps)
    elif vi_key=="EVI":  nir,red,bl=band(b["NIR"]),band(b["Red"]),band(b["Blue"]); r=2.5*(nir-red)/(nir+6*red-7.5*bl+1+eps)
    elif vi_key=="SAVI": nir,red=band(b["NIR"]),band(b["Red"]); r=1.5*(nir-red)/(nir+red+0.5+eps)
    elif vi_key=="OSAVI":nir,red=band(b["NIR"]),band(b["Red"]); r=(nir-red)/(nir+red+0.16+eps)
    elif vi_key=="VARI": g,red,bl=band(b["Green"]),band(b["Red"]),band(b["Blue"]); r=(g-red)/(g+red-bl+eps)
    elif vi_key=="NDRE": nir,re=band(b["NIR"]),band(b["RedEdge"]); r=(nir-re)/(nir+re+eps)
    elif vi_key=="CIre": nir,re=band(b["NIR"]),band(b["RedEdge"]); r=(nir/(re+eps))-1
    elif vi_key=="MCARI":re,red,g=band(b["RedEdge"]),band(b["Red"]),band(b["Green"]); r=((re-red)-0.2*(re-g))*(re/(red+eps))
    elif vi_key=="REIP": r670,r700,r740,r780=band(b["R670"]),band(b["R700"]),band(b["R740"]),band(b["R780"]); r=700+40*(((r670+r780)/2)-r700)/(r740-r700+eps)
    elif vi_key=="RENDVI":r750,r705=band(b["R750"]),band(b["R705"]); r=(r750-r705)/(r750+r705+eps)
    elif vi_key=="NDWI": g,nir=band(b["Green"]),band(b["NIR"]); r=(g-nir)/(g+nir+eps)
    elif vi_key=="WBI":  r900,r970=band(b["R900"]),band(b["R970"]); r=r900/(r970+eps)
    elif vi_key=="PSRI": red,g,re=band(b["Red"]),band(b["Green"]),band(b["RedEdge"]); r=(red-g)/(re+eps)
    elif vi_key=="PRI":  r531,r570=band(b["R531"]),band(b["R570"]); r=(r531-r570)/(r531+r570+eps)
    elif vi_key=="SIPI": nir,bl,red=band(b["NIR"]),band(b["Blue"]),band(b["Red"]); r=(nir-bl)/(nir-red+eps)
    elif vi_key=="MTVI2":
        nir,red,g=band(b["NIR"]),band(b["Red"]),band(b["Green"])
        num=1.5*(1.2*(nir-g)-2.5*(red-g))
        den=np.sqrt(np.maximum((2*nir+1)**2-(6*nir-5*np.sqrt(np.maximum(red,0)))-0.5,eps))
        r=num/den
    else: r=np.zeros((arr.shape[0],arr.shape[1]))
    return r.astype(np.float32)

def compute_custom_vi(arr, wavelengths, formula):
    import re
    bands_used = re.findall(r'B(\d+(?:\.\d+)?)', formula)
    local_vars = {f"B{nm}": arr[:,:,find_band(wavelengths,float(nm))].astype(np.float64) for nm in bands_used}
    local_vars["np"] = np; local_vars["eps"] = 1e-6
    try:
        return np.array(eval(formula, {"__builtins__":{}}, local_vars), dtype=np.float32)
    except Exception as e:
        st.error(f"Formula error: {e}"); return None

def make_rgb(arr, wavelengths, r_nm=660, g_nm=550, b_nm=480, percentile=2):
    def stretch(ch):
        lo,hi = np.percentile(ch[np.isfinite(ch)],[percentile,100-percentile])
        return np.clip((ch-lo)/(hi-lo+1e-6),0,1)
    r = stretch(arr[:,:,find_band(wavelengths,r_nm)])
    g = stretch(arr[:,:,find_band(wavelengths,g_nm)])
    b = stretch(arr[:,:,find_band(wavelengths,b_nm)])
    return np.stack([r,g,b],axis=-1)

def vi_statistics(vi_array):
    data = vi_array.flatten(); data = data[np.isfinite(data)]
    if len(data)==0: return {}
    return {"Min":float(np.min(data)),"Max":float(np.max(data)),"Mean":float(np.mean(data)),
            "Median":float(np.median(data)),"Std Dev":float(np.std(data)),
            "10th %ile":float(np.percentile(data,10)),"25th %ile":float(np.percentile(data,25)),
            "75th %ile":float(np.percentile(data,75)),"90th %ile":float(np.percentile(data,90)),
            "Valid pixels":int(len(data))}

def array_to_png_b64(arr_2d, colormap="RdYlGn", vmin=None, vmax=None):
    """Render a 2D float array to a base64 PNG using matplotlib colormap."""
    cmap = plt.get_cmap(colormap)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    rgba = cmap(norm(arr_2d))
    rgba_uint8 = (rgba * 255).astype(np.uint8)
    img_pil = Image.fromarray(rgba_uint8, "RGBA")
    buf = io.BytesIO(); img_pil.save(buf, format="PNG"); buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def rgb_array_to_png_b64(rgb_float):
    """Convert float [0,1] RGB array to base64 PNG."""
    rgb_uint8 = (np.clip(rgb_float,0,1)*255).astype(np.uint8)
    img_pil = Image.fromarray(rgb_uint8, "RGB")
    buf = io.BytesIO(); img_pil.save(buf, format="PNG"); buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def vi_to_plotly_fig(vi_array, vi_name, colormap="RdYlGn", vmin=None, vmax=None, click_enabled=False):
    vi_info = VI_LIBRARY.get(vi_name,{})
    r = vi_info.get("range",(-1,1))
    vmin = vmin if vmin is not None else r[0]
    vmax = vmax if vmax is not None else r[1]
    fig = px.imshow(vi_array, color_continuous_scale=colormap, zmin=vmin, zmax=vmax,
                    aspect="equal", labels={"color":vi_name},
                    title=f"{vi_name} — {vi_info.get('full_name','')}")
    fig.update_layout(margin=dict(l=0,r=0,t=40,b=0),
                      coloraxis_colorbar=dict(title=vi_name,thickness=15), height=500)
    if click_enabled:
        fig.update_traces(hovertemplate="Row: %{y}<br>Col: %{x}<br>Value: %{z:.4f}<extra></extra>")
    return fig

def rgb_to_plotly_fig(rgb_array, click_enabled=False):
    fig = px.imshow(rgb_array, aspect="equal", title="RGB Composite")
    fig.update_layout(margin=dict(l=0,r=0,t=40,b=0), height=500)
    if click_enabled:
        fig.update_traces(hovertemplate="Row: %{y}<br>Col: %{x}<extra></extra>")
    return fig

# ══════════════════════════════════════════════════════════════════════════════
#  SESSION STATE
# ══════════════════════════════════════════════════════════════════════════════
for key, default in [("custom_vis",{}), ("vi_cache",{}),
                     ("loaded_hdr",None), ("loaded_bsq",None),
                     ("click_row",None), ("click_col",None)]:
    if key not in st.session_state:
        st.session_state[key] = default

# ══════════════════════════════════════════════════════════════════════════════
#  SIDEBAR — FILE LOADING
# ══════════════════════════════════════════════════════════════════════════════
st.title("🌱 Hyperspectral Field Viewer")
st.caption("Load your BSQ + HDR file pair to explore vegetation indices, spectra, and field statistics.")

with st.sidebar:
    st.header("📂 Load Data")
    load_mode = st.radio("Load method", ["📁 File path (large files)", "⬆️ Upload (small files)"])

    if load_mode == "📁 File path (large files)":
        st.caption("Paste the full path to your HDR file. BSQ is found automatically.")
        hdr_path_input = st.text_input("Path to .HDR file", placeholder=r"D:\data\field_july.hdr")
        if st.button("🔍 Load File", type="primary"):
            hdr_path_input = hdr_path_input.strip().strip('"').strip("'")
            if not os.path.isfile(hdr_path_input):
                st.error(f"HDR not found:\n`{hdr_path_input}`")
            else:
                base = os.path.splitext(hdr_path_input)[0]
                found = next((base+ext for ext in [".bsq",".BSQ",".img",".IMG",".dat",".raw"] if os.path.isfile(base+ext)), None)
                if not found:
                    st.error("Could not find matching BSQ/img file alongside HDR.")
                else:
                    st.session_state["loaded_hdr"] = hdr_path_input
                    st.session_state["loaded_bsq"] = found
                    st.session_state["vi_cache"] = {}
        hdr_file = bsq_file = None
    else:
        hdr_file = st.file_uploader("HDR file", type=["hdr","HDR"])
        bsq_file = st.file_uploader("BSQ file", type=["bsq","BSQ","img","dat","raw"])

    st.divider()
    st.header("⚙️ Display Settings")
    stretch_pct = st.slider("RGB stretch percentile", 0, 10, 2)
    st.divider()
    st.markdown("**Built with** `spectral` · `streamlit` · `plotly`")

# ── Resolve ready state ────────────────────────────────────────────────────────
ready = False
if load_mode == "📁 File path (large files)" and st.session_state["loaded_hdr"]:
    ready = True
elif load_mode == "⬆️ Upload (small files)" and hdr_file and bsq_file:
    ready = True

if not ready:
    st.info("👈 Use the sidebar to load your HDR + BSQ files to get started.")
    with st.expander("ℹ️ What this tool does"):
        st.markdown("""
- 🖼️ **Overview** — file metadata + RGB quick-look
- 📊 **VI Browser** — 17 vegetation indices + RGB composite, with histogram and stats
- 🔀 **Compare** — side-by-side swipe slider between any two views
- ✏️ **Custom VI** — write your own index formula
- 📈 **Spectra** — click a pixel on the image or draw a polygon to extract spectra
- 📉 **Statistics** — field-wide stats and distribution plots

**VIs included:** NDVI · NDRE · GNDVI · EVI · SAVI · OSAVI · CIre · MCARI · REIP · NDWI · WBI · PSRI · PRI · SIPI · VARI · MTVI2 · RENDVI
        """)
    st.stop()

# ── Open image ─────────────────────────────────────────────────────────────────
@st.cache_resource(show_spinner="Opening file…")
def open_from_path(hdr, bsq):
    return envi.open(hdr, bsq)

@st.cache_resource(show_spinner="Preparing uploaded files…")
def open_from_upload(hdr_bytes, bsq_bytes):
    tmpdir = tempfile.mkdtemp()
    hdr_path = os.path.join(tmpdir,"image.hdr")
    bsq_path = os.path.join(tmpdir,"image.bsq")
    open(hdr_path,"wb").write(hdr_bytes)
    open(bsq_path,"wb").write(bsq_bytes)
    return envi.open(hdr_path, bsq_path), hdr_path, bsq_path

if load_mode == "📁 File path (large files)":
    img = open_from_path(st.session_state["loaded_hdr"], st.session_state["loaded_bsq"])
    st.sidebar.success(f"✅ `{Path(st.session_state['loaded_hdr']).name}`")
else:
    img, _, _ = open_from_upload(hdr_file.read(), bsq_file.read())

wavelengths = tuple(get_wavelengths(img))
wl_min, wl_max = min(wavelengths), max(wavelengths)

def vi_available(vi_key):
    return all(wl_min-30 <= nm <= wl_max+30 for nm in VI_LIBRARY[vi_key]["bands"].values())

available_vis = [k for k in VI_LIBRARY if vi_available(k)]
unavailable_vis = [k for k in VI_LIBRARY if not vi_available(k)]

# ── All renderable options (VIs + RGB presets) ────────────────────────────────
RGB_PRESETS = {
    "RGB — True color (~660/550/480 nm)":   (660, 550, 480),
    "RGB — Color infrared (~800/660/550 nm)": (800, 660, 550),
    "RGB — Red edge (~750/700/670 nm)":      (750, 700, 670),
}
all_view_options = list(RGB_PRESETS.keys()) + available_vis + list(st.session_state["custom_vis"].keys())

# ══════════════════════════════════════════════════════════════════════════════
#  TABS
# ══════════════════════════════════════════════════════════════════════════════
tab_overview, tab_vi, tab_compare, tab_custom, tab_spectra, tab_stats = st.tabs([
    "📋 Overview", "📊 VI Browser", "🔀 Compare", "✏️ Custom VI", "📈 Spectra", "📉 Statistics"
])

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 1 — OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════
with tab_overview:
    meta = get_meta_summary(img)

    # Metrics — only shown when data is real (not zero)
    c1,c2,c3,c4 = st.columns(4)
    c1.metric("Rows", f"{meta['Rows']:,}")
    c2.metric("Columns", f"{meta['Columns']:,}")
    c3.metric("Bands", meta["Bands"])
    c4.metric("Wavelength range", f"{meta['Wavelength min (nm)']}–{meta['Wavelength max (nm)']} nm")

    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("**File metadata**")
        st.dataframe(pd.DataFrame(list(meta.items()),columns=["Field","Value"]),
                     use_container_width=True, hide_index=True)
    with col_b:
        st.markdown("**Band list (first 100)**")
        st.dataframe(pd.DataFrame({"Band #":range(1,min(101,img.nbands+1)),
                                   "Wavelength (nm)":[f"{w:.2f}" for w in wavelengths[:100]]}),
                     use_container_width=True, hide_index=True, height=320)

    st.divider()
    st.markdown("**Quick-look RGB** — rendered from the loaded data")
    rgb_preset = st.selectbox("Composite preset", list(RGB_PRESETS.keys()), key="ov_rgb")
    if st.button("🖼️ Render RGB", key="ov_rgb_btn"):
        with st.spinner("Rendering…"):
            arr = read_as_array(img)
            r_nm,g_nm,b_nm = RGB_PRESETS[rgb_preset]
            rgb = make_rgb(arr, wavelengths, r_nm, g_nm, b_nm, stretch_pct)
            st.plotly_chart(rgb_to_plotly_fig(rgb), use_container_width=True)

    st.divider()
    st.markdown("**Available vegetation indices**")
    cats = sorted(set(VI_LIBRARY[k]["category"] for k in available_vis))
    for cat in cats:
        st.markdown(f"*{cat}*")
        vis_in_cat = [k for k in available_vis if VI_LIBRARY[k]["category"]==cat]
        cols = st.columns(min(len(vis_in_cat),4))
        for i,vi_key in enumerate(vis_in_cat):
            with cols[i%4]:
                vi = VI_LIBRARY[vi_key]
                st.markdown(f"<div class='vi-card'><strong>{vi_key}</strong><br>"
                            f"<span style='font-size:0.78rem;color:#444'>{vi['full_name']}</span></div>",
                            unsafe_allow_html=True)
    if unavailable_vis:
        with st.expander(f"⚠️ {len(unavailable_vis)} VIs unavailable (outside sensor range)"):
            for k in unavailable_vis:
                req = ", ".join(f"{nm} nm" for nm in VI_LIBRARY[k]["bands"].values())
                st.markdown(f"**{k}** — requires {req}")

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 2 — VI BROWSER  (includes RGB options)
# ══════════════════════════════════════════════════════════════════════════════
with tab_vi:
    st.subheader("VI Browser")

    all_options = list(RGB_PRESETS.keys()) + available_vis + list(st.session_state["custom_vis"].keys())

    cats_for_filter = ["All", "RGB Composites"] + sorted(set(
        VI_LIBRARY[k]["category"] if k in VI_LIBRARY else "Custom"
        for k in (available_vis + list(st.session_state["custom_vis"].keys()))
    ))
    cat_filter = st.selectbox("Filter by category", cats_for_filter, key="vi_cat")

    if cat_filter == "All":
        shown = all_options
    elif cat_filter == "RGB Composites":
        shown = list(RGB_PRESETS.keys())
    else:
        shown = [k for k in (available_vis + list(st.session_state["custom_vis"].keys()))
                 if (VI_LIBRARY[k]["category"] if k in VI_LIBRARY else "Custom") == cat_filter]

    col_left, col_right = st.columns([1,3])
    with col_left:
        selected = st.radio("Select view", shown, key="vi_select")

    with col_right:
        is_rgb = selected in RGB_PRESETS
        is_custom = selected in st.session_state["custom_vis"]

        if is_rgb:
            r_nm,g_nm,b_nm = RGB_PRESETS[selected]
            st.markdown(f"**{selected}**")
            st.caption(f"Bands: R={r_nm} nm · G={g_nm} nm · B={b_nm} nm")
            custom_r = st.selectbox("🔴 Override Red band", range(len(wavelengths)),
                                    index=find_band(wavelengths,r_nm),
                                    format_func=lambda i: f"{wavelengths[i]:.1f} nm", key="ov_r")
            custom_g = st.selectbox("🟢 Override Green band", range(len(wavelengths)),
                                    index=find_band(wavelengths,g_nm),
                                    format_func=lambda i: f"{wavelengths[i]:.1f} nm", key="ov_g")
            custom_b = st.selectbox("🔵 Override Blue band", range(len(wavelengths)),
                                    index=find_band(wavelengths,b_nm),
                                    format_func=lambda i: f"{wavelengths[i]:.1f} nm", key="ov_b")
            render_key = f"rgb_{selected}_{custom_r}_{custom_g}_{custom_b}"
        elif is_custom:
            formula = st.session_state["custom_vis"][selected]
            st.markdown(f"**Custom: {selected}**")
            st.code(formula)
            cmap_choice = st.selectbox("Colormap",["RdYlGn","YlGn","RdBu","RdYlGn_r","plasma","viridis"],key="vi_cmap_c")
            vmin_val = st.number_input("Min",value=-1.0,step=0.05,key="vi_vmin_c")
            vmax_val = st.number_input("Max",value=1.0,step=0.05,key="vi_vmax_c")
            render_key = f"custom_{selected}"
        else:
            vi_info = VI_LIBRARY[selected]
            st.markdown(f"**{vi_info['full_name']}**")
            st.markdown(f"`{vi_info['formula_display']}`")
            st.caption(vi_info["description"])
            st.caption("Bands used: " + ", ".join(f"{nm} nm" for nm in vi_info["bands"].values()))
            cmap_choice = st.selectbox("Colormap",["RdYlGn","YlGn","RdBu","RdYlGn_r","plasma","viridis"],key="vi_cmap")
            default_range = vi_info.get("range",(-1,1))
            col_vmin, col_vmax = st.columns(2)
            with col_vmin: vmin_val = st.number_input("Min",value=float(default_range[0]),step=0.05,key="vi_vmin")
            with col_vmax: vmax_val = st.number_input("Max",value=float(default_range[1]),step=0.05,key="vi_vmax")
            render_key = selected

        if st.button(f"📊 Compute & View", type="primary", key="vi_render"):
            with st.spinner("Computing…"):
                arr = read_as_array(img)
                if is_rgb:
                    r_out = wavelengths[custom_r]; g_out = wavelengths[custom_g]; b_out = wavelengths[custom_b]
                    result_rgb = make_rgb(arr, wavelengths, r_out, g_out, b_out, stretch_pct)
                    st.session_state["vi_cache"][render_key] = ("rgb", result_rgb)
                    fig = rgb_to_plotly_fig(result_rgb, click_enabled=True)
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    if is_custom:
                        vi_map = compute_custom_vi(arr, wavelengths, st.session_state["custom_vis"][selected])
                    else:
                        vi_map = compute_vi(arr, selected, wavelengths)
                    if vi_map is not None:
                        st.session_state["vi_cache"][render_key] = ("vi", vi_map, cmap_choice, vmin_val, vmax_val)
                        fig = vi_to_plotly_fig(vi_map, selected, cmap_choice, vmin_val, vmax_val, click_enabled=True)
                        st.plotly_chart(fig, use_container_width=True)

                        stats = vi_statistics(vi_map)
                        cc1,cc2,cc3,cc4 = st.columns(4)
                        cc1.metric("Mean", f"{stats.get('Mean',0):.4f}")
                        cc2.metric("Median", f"{stats.get('Median',0):.4f}")
                        cc3.metric("Std Dev", f"{stats.get('Std Dev',0):.4f}")
                        cc4.metric("Valid pixels", f"{stats.get('Valid pixels',0):,}")

                        data_flat = vi_map.flatten()
                        data_flat = data_flat[np.isfinite(data_flat) & (data_flat>=vmin_val) & (data_flat<=vmax_val)]
                        fig_hist = px.histogram(x=data_flat, nbins=100,
                                                title=f"{selected} distribution",
                                                labels={"x":selected,"y":"Count"},
                                                color_discrete_sequence=["#2d7a3a"])
                        fig_hist.update_layout(height=260, margin=dict(t=40,b=0))
                        st.plotly_chart(fig_hist, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 3 — SWIPE COMPARE
# ══════════════════════════════════════════════════════════════════════════════
with tab_compare:
    st.subheader("🔀 Side-by-side Swipe Comparison")
    st.caption("Select two views, render them, then drag the slider to reveal one or the other.")

    cached_keys = list(st.session_state["vi_cache"].keys())

    if len(cached_keys) < 2:
        st.info("You need at least **2 computed views** in the VI Browser first — then come back here to compare them.")
        st.markdown("Computed so far: " + (", ".join(f"`{k}`" for k in cached_keys) if cached_keys else "*none yet*"))
    else:
        col_l, col_r = st.columns(2)
        with col_l:
            left_key = st.selectbox("⬅️ Left image", cached_keys, index=0, key="cmp_left")
        with col_r:
            right_key = st.selectbox("➡️ Right image", cached_keys,
                                     index=min(1,len(cached_keys)-1), key="cmp_right")

        display_w = st.slider("Display width (px)", 400, 1400, 800, step=50)

        if st.button("🔀 Generate Swipe View", type="primary"):
            def render_to_b64(cache_entry, target_w):
                kind = cache_entry[0]
                if kind == "rgb":
                    rgb_float = cache_entry[1]
                    h0,w0 = rgb_float.shape[:2]
                    new_h = int(h0 * target_w / w0)
                    pil = Image.fromarray((np.clip(rgb_float,0,1)*255).astype(np.uint8),"RGB")
                    pil = pil.resize((target_w, new_h), Image.LANCZOS)
                    buf = io.BytesIO(); pil.save(buf,"PNG"); buf.seek(0)
                    return base64.b64encode(buf.read()).decode(), new_h
                else:
                    vi_map, cmap, vmin, vmax = cache_entry[1], cache_entry[2], cache_entry[3], cache_entry[4]
                    b64 = array_to_png_b64(vi_map, cmap, vmin, vmax)
                    h0,w0 = vi_map.shape
                    new_h = int(h0 * target_w / w0)
                    raw = base64.b64decode(b64)
                    pil = Image.open(io.BytesIO(raw)).resize((target_w, new_h), Image.LANCZOS)
                    buf = io.BytesIO(); pil.save(buf,"PNG"); buf.seek(0)
                    return base64.b64encode(buf.read()).decode(), new_h

            def make_colorbar_b64(cmap_name, vmin, vmax, label, bar_w=28, bar_h=220):
                """Render a vertical colorbar as a base64 PNG."""
                fig_cb, ax = plt.subplots(figsize=(1.1, bar_h/72), dpi=96)
                norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
                cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap_name)),
                                  cax=ax, orientation="vertical")
                cb.set_label(label, fontsize=8, color="white")
                cb.ax.yaxis.set_tick_params(color="white", labelsize=7, labelcolor="white")
                ax.set_facecolor("#111"); fig_cb.patch.set_facecolor("#111")
                buf = io.BytesIO(); fig_cb.savefig(buf, format="png", bbox_inches="tight",
                                                   facecolor="#111", dpi=96)
                plt.close(fig_cb); buf.seek(0)
                return base64.b64encode(buf.read()).decode()

            def entry_legend_html(entry, key, side):
                """Return an HTML snippet with label + colorbar (or just label for RGB)."""
                kind = entry[0]
                label_html = (f'<div style="background:rgba(0,0,0,0.6);color:white;padding:3px 8px;'
                               f'border-radius:4px;font-size:12px;font-family:sans-serif;'
                               f'white-space:nowrap;">{key}</div>')
                if kind == "rgb":
                    return label_html, ""
                # entry is ("vi", vi_map, cmap, vmin, vmax)
                cmap, vmin, vmax = entry[2], entry[3], entry[4]
                cb_b64 = make_colorbar_b64(cmap, vmin, vmax, key)
                if side == "left":
                    cb_style = "position:absolute;bottom:8px;left:8px;z-index:10;"
                else:
                    cb_style = "position:absolute;bottom:8px;right:8px;z-index:10;"
                cb_html = (f'<div style="{cb_style}">'
                           f'<img src="data:image/png;base64,{cb_b64}" '
                           f'style="height:180px;border-radius:4px;opacity:0.92;" />'
                           f'</div>')
                return label_html, cb_html

            with st.spinner("Rendering swipe view…"):
                left_entry = st.session_state["vi_cache"][left_key]
                right_entry = st.session_state["vi_cache"][right_key]
                ref_arr = left_entry[1]
                ref_h, ref_w = ref_arr.shape[0], ref_arr.shape[1]
                disp_h = int(ref_h * display_w / ref_w)

                lb64, _ = render_to_b64(left_entry, display_w)
                rb64, _ = render_to_b64(right_entry, display_w)

                l_label_html, l_cb_html = entry_legend_html(left_entry, left_key, "left")
                r_label_html, r_cb_html = entry_legend_html(right_entry, right_key, "right")

            swipe_html = f"""
<div style="position:relative;width:{display_w}px;height:{disp_h}px;
            overflow:hidden;cursor:col-resize;user-select:none;border:1px solid #ccc;border-radius:6px;"
     id="swipe-container">
  <img src="data:image/png;base64,{rb64}"
       style="position:absolute;top:0;left:0;width:{display_w}px;height:{disp_h}px;display:block;" />
  <div style="position:absolute;top:8px;right:8px;z-index:10;">{r_label_html}</div>
  {r_cb_html}
  <div id="left-clip"
       style="position:absolute;top:0;left:0;width:{display_w//2}px;height:{disp_h}px;overflow:hidden;">
    <img src="data:image/png;base64,{lb64}"
         style="width:{display_w}px;height:{disp_h}px;display:block;" />
    <div style="position:absolute;top:8px;left:8px;z-index:10;">{l_label_html}</div>
    {l_cb_html}
  </div>
  <div id="divider"
       style="position:absolute;top:0;left:{display_w//2}px;width:3px;height:{disp_h}px;
              background:white;box-shadow:0 0 6px rgba(0,0,0,0.6);cursor:col-resize;z-index:20;">
    <div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);
                background:white;border-radius:50%;width:32px;height:32px;
                display:flex;align-items:center;justify-content:center;
                box-shadow:0 0 4px rgba(0,0,0,0.5);font-size:16px;font-weight:bold;">⟺</div>
  </div>
</div>
<script>
(function(){{
  const container = document.getElementById('swipe-container');
  const clip      = document.getElementById('left-clip');
  const divider   = document.getElementById('divider');
  let dragging = false;
  function setPos(x) {{
    const rect = container.getBoundingClientRect();
    let pos = Math.max(0, Math.min(x - rect.left, {display_w}));
    clip.style.width    = pos + 'px';
    divider.style.left  = pos + 'px';
  }}
  divider.addEventListener('mousedown',  e => {{ dragging=true; e.preventDefault(); }});
  document.addEventListener('mouseup',   () => dragging=false);
  document.addEventListener('mousemove', e => {{ if(dragging) setPos(e.clientX); }});
  divider.addEventListener('touchstart', e => {{ dragging=true; e.preventDefault(); }});
  document.addEventListener('touchend',  () => dragging=false);
  document.addEventListener('touchmove', e => {{ if(dragging) setPos(e.touches[0].clientX); }});
  container.addEventListener('click',    e => setPos(e.clientX));
}})();
</script>
"""
            st.components.v1.html(swipe_html, height=disp_h+20, scrolling=False)
            st.caption(f"⬅️ {left_key}  |  {right_key} ➡️  — drag the white line or click to reveal")

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 4 — CUSTOM VI
# ══════════════════════════════════════════════════════════════════════════════
with tab_custom:
    st.subheader("✏️ Add a Custom Vegetation Index")
    st.markdown("""
Write your formula using **`B<wavelength>`** (e.g. `B800` = band nearest 800 nm).
`np` and `eps` (1e-6) are available.

**Examples:**
- `(B800 - B670) / (B800 + B670 + eps)` — custom NDVI
- `np.log(B800 / (B670 + eps))` — log ratio
""")
    st.caption(f"Sensor range: **{wl_min:.1f} – {wl_max:.1f} nm**")

    with st.form("custom_vi_form"):
        vi_name_in = st.text_input("Index name (no spaces)", placeholder="MyVI")
        vi_formula_in = st.text_area("Formula", placeholder="(B800 - B670) / (B800 + B670 + eps)", height=70)
        if st.form_submit_button("💾 Save"):
            name = vi_name_in.strip()
            if not name or " " in name:
                st.error("Name required, no spaces.")
            elif not vi_formula_in.strip():
                st.error("Formula required.")
            else:
                st.session_state["custom_vis"][name] = vi_formula_in.strip()
                st.success(f"Saved **{name}** — find it in VI Browser under Custom.")

    if st.session_state["custom_vis"]:
        st.divider()
        for name, formula in list(st.session_state["custom_vis"].items()):
            ca, cb = st.columns([4,1])
            with ca:
                st.markdown(f"<div class='vi-card'><strong>{name}</strong><br>"
                            f"<code style='font-size:0.8rem'>{formula}</code></div>",
                            unsafe_allow_html=True)
            with cb:
                if st.button("🗑️", key=f"del_{name}"):
                    del st.session_state["custom_vis"][name]
                    st.rerun()

    with st.expander("📋 Wavelength reference"):
        st.dataframe(pd.DataFrame({
            "Band": range(1,img.nbands+1),
            "Wavelength (nm)": [f"{w:.2f}" for w in wavelengths],
            "Use as": [f"B{round(w)}" for w in wavelengths],
        }), use_container_width=True, hide_index=True, height=280)

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 5 — SPECTRA  (click-to-inspect + polygon ROI)
# ══════════════════════════════════════════════════════════════════════════════
with tab_spectra:
    st.subheader("📈 Spectrum Inspector")

    spectra_mode = st.radio("Mode", ["🖱️ Click pixel", "🔲 Draw polygon / ROI"],
                            horizontal=True)

    if spectra_mode == "🖱️ Click pixel":
        st.caption("Use the VI Browser to render a map, then enter the coordinates shown in the hover tooltip below — or type them directly.")

        # Show a quick reference image
        ref_mode = st.selectbox("Reference image", list(RGB_PRESETS.keys()) + available_vis[:3], key="sp_ref")
        if st.button("🖼️ Load reference image", key="sp_load_ref"):
            with st.spinner("Rendering…"):
                arr = read_as_array(img)
                if ref_mode in RGB_PRESETS:
                    r_nm,g_nm,b_nm = RGB_PRESETS[ref_mode]
                    rgb = make_rgb(arr, wavelengths, r_nm, g_nm, b_nm, stretch_pct)
                    fig_ref = rgb_to_plotly_fig(rgb, click_enabled=True)
                else:
                    vi_map = compute_vi(arr, ref_mode, wavelengths)
                    vi_info = VI_LIBRARY[ref_mode]
                    fig_ref = vi_to_plotly_fig(vi_map, ref_mode, vi_info["colormap"],
                                               vi_info["range"][0], vi_info["range"][1], click_enabled=True)
                st.plotly_chart(fig_ref, use_container_width=True,
                                key="spectra_ref_plot")
                st.caption("📍 Hover over the image to read Row/Col coordinates, then enter them below.")

        col1, col2 = st.columns(2)
        with col1:
            px_row = st.number_input("Row", 0, img.nrows-1, img.nrows//2, key="sp_row")
        with col2:
            px_col = st.number_input("Column", 0, img.ncols-1, img.ncols//2, key="sp_col")

        # Multi-pixel comparison
        n_compare = st.number_input("Extra pixels to compare (0–4)", 0, 4, 0)
        compare_pixels = []
        for i in range(int(n_compare)):
            cc1,cc2 = st.columns(2)
            with cc1: r2 = st.number_input(f"Row {i+2}", 0, img.nrows-1, img.nrows//3, key=f"spr{i}")
            with cc2: c2 = st.number_input(f"Col {i+2}", 0, img.ncols-1, img.ncols//3, key=f"spc{i}")
            compare_pixels.append((int(r2),int(c2)))

        if st.button("📈 Plot Spectrum", type="primary", key="sp_plot"):
            with st.spinner("Extracting…"):
                arr = read_as_array(img)
                fig = go.Figure()
                colors_list = ["#2d7a3a","#e05c2f","#2f6ee0","#c44db5","#e0b52f"]
                all_pts = [(int(px_row),int(px_col))] + compare_pixels
                for idx,(rr,cc) in enumerate(all_pts):
                    sp = arr[rr,cc,:]
                    fig.add_trace(go.Scatter(x=list(wavelengths),y=sp.tolist(),mode="lines",
                                             name=f"({rr},{cc})",
                                             line=dict(color=colors_list[idx%5],width=2)))
                key_wls = {"Blue":480,"Green":550,"Red":670,"Red Edge":720,"NIR":800}
                for wname,nm in key_wls.items():
                    if wl_min<=nm<=wl_max:
                        fig.add_vline(x=nm,line_dash="dot",line_color="gray",
                                      annotation_text=wname,line_width=0.8)
                fig.update_layout(title="Reflectance spectrum",xaxis_title="Wavelength (nm)",
                                  yaxis_title="Reflectance",height=430,
                                  legend=dict(orientation="h",yanchor="bottom",y=1.02),
                                  margin=dict(t=60))
                st.plotly_chart(fig, use_container_width=True)

                # Key band values table
                key_data = {}
                for bn,nm in key_wls.items():
                    if wl_min<=nm<=wl_max:
                        bidx = find_band(wavelengths,nm)
                        key_data[f"{bn} ({nm} nm)"] = f"{arr[int(px_row),int(px_col),bidx]:.4f}"
                st.dataframe(pd.DataFrame(list(key_data.items()),columns=["Band","Reflectance"]),
                             use_container_width=True, hide_index=True)

    else:  # ── Polygon / ROI mode ──────────────────────────────────────────
        st.caption(
            "Click points on the image to build a polygon ROI, then extract the mean spectrum "
            "of all pixels inside. Click **Add ROI** to save it and start a new one."
        )

        # ── session state for polygon drawing ─────────────────────────────────
        if "roi_current_pts" not in st.session_state:
            st.session_state["roi_current_pts"] = []   # list of (col, row) in image coords
        if "roi_saved" not in st.session_state:
            st.session_state["roi_saved"] = []         # list of {"label":str, "pts":[(c,r),...]}

        # ── pick reference background ──────────────────────────────────────────
        roi_ref = st.selectbox("Background image", list(RGB_PRESETS.keys()) + available_vis[:3], key="roi_ref")
        render_col, clear_col = st.columns([3,1])
        with render_col:
            render_bg = st.button("🖼️ Load / refresh image", key="roi_render")
        with clear_col:
            if st.button("🗑️ Clear all ROIs", key="roi_clear_all"):
                st.session_state["roi_current_pts"] = []
                st.session_state["roi_saved"] = []
                st.rerun()

        if render_bg or "roi_bg_key" not in st.session_state or st.session_state.get("roi_bg_key") != roi_ref:
            with st.spinner("Rendering background…"):
                arr_bg = read_as_array(img)
                if roi_ref in RGB_PRESETS:
                    r_nm,g_nm,b_nm = RGB_PRESETS[roi_ref]
                    bg_rgb = make_rgb(arr_bg, wavelengths, r_nm, g_nm, b_nm, stretch_pct)
                    st.session_state["roi_bg_img"] = bg_rgb
                    st.session_state["roi_bg_kind"] = "rgb"
                else:
                    vi_map_bg = compute_vi(arr_bg, roi_ref, wavelengths)
                    st.session_state["roi_bg_img"] = vi_map_bg
                    st.session_state["roi_bg_kind"] = "vi"
                    st.session_state["roi_bg_vi"] = roi_ref
                st.session_state["roi_bg_key"] = roi_ref

        # ── build Plotly figure with background + drawn points/polygon ─────────
        if "roi_bg_img" in st.session_state:
            bg = st.session_state["roi_bg_img"]
            bg_kind = st.session_state["roi_bg_kind"]

            if bg_kind == "rgb":
                fig_roi = px.imshow(bg, aspect="equal")
            else:
                vi_info_bg = VI_LIBRARY.get(st.session_state.get("roi_bg_vi",""), {})
                r_bg = vi_info_bg.get("range",(-1,1))
                fig_roi = px.imshow(bg, color_continuous_scale=vi_info_bg.get("colormap","RdYlGn"),
                                    zmin=r_bg[0], zmax=r_bg[1], aspect="equal")

            # Draw saved ROIs
            colors_roi = ["#FF6600","#00AAFF","#CC00FF","#00CC66","#FFCC00"]
            for si, saved in enumerate(st.session_state["roi_saved"]):
                pts = saved["pts"]
                if len(pts) >= 2:
                    xs = [p[0] for p in pts] + [pts[0][0]]
                    ys = [p[1] for p in pts] + [pts[0][1]]
                    c = colors_roi[si % len(colors_roi)]
                    fig_roi.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers",
                                                  line=dict(color=c, width=2),
                                                  marker=dict(color=c, size=6),
                                                  name=saved["label"], showlegend=True))

            # Draw current in-progress points
            cur = st.session_state["roi_current_pts"]
            if cur:
                xs_c = [p[0] for p in cur]
                ys_c = [p[1] for p in cur]
                fig_roi.add_trace(go.Scatter(x=xs_c, y=ys_c, mode="lines+markers",
                                              line=dict(color="white", width=1.5, dash="dot"),
                                              marker=dict(color="white", size=7, symbol="circle-open"),
                                              name="Current ROI (in progress)", showlegend=True))

            fig_roi.update_layout(
                height=520, margin=dict(l=0,r=0,t=30,b=0),
                title="Click to add polygon vertices — use buttons below to save or undo",
                dragmode="zoom",
            )
            fig_roi.update_traces(hovertemplate="Col: %{x}<br>Row: %{y}<extra></extra>")

            clicked = st.plotly_chart(fig_roi, use_container_width=True,
                                       on_select="rerun", key="roi_plot")

            # Handle click → add point
            if clicked and clicked.get("selection") and clicked["selection"].get("points"):
                for pt in clicked["selection"]["points"]:
                    col_coord = int(round(pt["x"]))
                    row_coord = int(round(pt["y"]))
                    col_coord = np.clip(col_coord, 0, img.ncols-1)
                    row_coord = np.clip(row_coord, 0, img.nrows-1)
                    new_pt = (col_coord, row_coord)
                    if new_pt not in st.session_state["roi_current_pts"]:
                        st.session_state["roi_current_pts"].append(new_pt)
                        st.rerun()

            # Show current point list
            cur = st.session_state["roi_current_pts"]
            if cur:
                pts_str = " → ".join(f"({p[0]},{p[1]})" for p in cur)
                st.caption(f"**Current polygon vertices ({len(cur)}):** {pts_str}")

            # Controls row
            btn1, btn2, btn3, btn4 = st.columns(4)
            with btn1:
                roi_label = st.text_input("ROI label", value=f"ROI {len(st.session_state['roi_saved'])+1}",
                                           key="roi_label_in", label_visibility="collapsed",
                                           placeholder="ROI label…")
            with btn2:
                if st.button("✅ Add ROI", type="primary", key="roi_add",
                             disabled=len(cur) < 3):
                    st.session_state["roi_saved"].append({"label": roi_label, "pts": list(cur)})
                    st.session_state["roi_current_pts"] = []
                    st.rerun()
            with btn3:
                if st.button("↩️ Undo last point", key="roi_undo", disabled=len(cur)==0):
                    st.session_state["roi_current_pts"].pop()
                    st.rerun()
            with btn4:
                if st.button("🗑️ Clear current", key="roi_clear_cur", disabled=len(cur)==0):
                    st.session_state["roi_current_pts"] = []
                    st.rerun()

            # Extract spectra from saved ROIs
            if st.session_state["roi_saved"]:
                st.divider()
                if st.button("📈 Extract & Plot ROI Spectra", type="primary", key="roi_extract"):
                    with st.spinner("Extracting spectra…"):
                        from PIL import ImageDraw as PIDraw
                        arr_sp = read_as_array(img)
                        fig_sp = go.Figure()
                        wl_list = list(wavelengths)

                        for si, saved in enumerate(st.session_state["roi_saved"]):
                            pts = saved["pts"]
                            label = saved["label"]
                            color = colors_roi[si % len(colors_roi)]

                            if len(pts) < 3:
                                st.warning(f"{label}: need ≥ 3 points, skipping.")
                                continue

                            # Build mask
                            mask_img = Image.new("L", (img.ncols, img.nrows), 0)
                            PIDraw.Draw(mask_img).polygon(pts, fill=255)
                            mask = np.array(mask_img) > 0
                            n_px = int(mask.sum())

                            if n_px == 0:
                                st.warning(f"{label}: polygon encloses no pixels (too small?).")
                                continue

                            pixels = arr_sp[mask, :]
                            mean_sp = np.nanmean(pixels, axis=0)
                            std_sp  = np.nanstd(pixels, axis=0)

                            # Shaded ±1 std band
                            fig_sp.add_trace(go.Scatter(
                                x=wl_list + wl_list[::-1],
                                y=(mean_sp+std_sp).tolist() + (mean_sp-std_sp).tolist()[::-1],
                                fill="toself",
                                fillcolor=color + "33",
                                line=dict(color="rgba(0,0,0,0)"),
                                showlegend=False, hoverinfo="skip"
                            ))
                            fig_sp.add_trace(go.Scatter(
                                x=wl_list, y=mean_sp.tolist(), mode="lines",
                                name=f"{label} ({n_px:,} px)",
                                line=dict(color=color, width=2)
                            ))

                        key_wls = {"Blue":480,"Green":550,"Red":670,"Red Edge":720,"NIR":800}
                        for wname,nm in key_wls.items():
                            if wl_min<=nm<=wl_max:
                                fig_sp.add_vline(x=nm, line_dash="dot", line_color="gray",
                                                  annotation_text=wname, line_width=0.8)
                        fig_sp.update_layout(
                            title="Mean ROI spectra  (shaded = ±1 std dev)",
                            xaxis_title="Wavelength (nm)", yaxis_title="Reflectance",
                            height=450,
                            legend=dict(orientation="h", yanchor="bottom", y=1.02),
                            margin=dict(t=60)
                        )
                        st.plotly_chart(fig_sp, use_container_width=True)
        else:
            st.info("Click **Load / refresh image** to display the background, then click on it to add polygon vertices.")

# ══════════════════════════════════════════════════════════════════════════════
#  TAB 6 — STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
with tab_stats:
    st.subheader("📉 Field Statistics")

    st.markdown("### Vegetation Index Statistics")
    vis_to_stat = st.multiselect(
        "Select VIs",
        available_vis + list(st.session_state["custom_vis"].keys()),
        default=(["NDVI","NDRE"] if "NDVI" in available_vis and "NDRE" in available_vis
                 else available_vis[:min(2,len(available_vis))])
    )

    if st.button("📊 Compute Statistics", type="primary"):
        if not vis_to_stat:
            st.warning("Select at least one VI.")
        else:
            arr = read_as_array(img)
            all_stats = {}
            prog = st.progress(0,"Computing…")
            for i,vi_key in enumerate(vis_to_stat):
                prog.progress((i+1)/len(vis_to_stat), text=f"Computing {vi_key}…")
                vi_map = (compute_custom_vi(arr, wavelengths, st.session_state["custom_vis"][vi_key])
                          if vi_key in st.session_state["custom_vis"]
                          else compute_vi(arr, vi_key, wavelengths))
                if vi_map is not None:
                    all_stats[vi_key] = vi_statistics(vi_map)
                    st.session_state["vi_cache"][vi_key] = ("vi", vi_map,
                        VI_LIBRARY.get(vi_key,{}).get("colormap","RdYlGn"),
                        *VI_LIBRARY.get(vi_key,{}).get("range",(-1,1)))
            prog.empty()

            if all_stats:
                df_stats = pd.DataFrame(all_stats).T.round(5)
                st.dataframe(df_stats, use_container_width=True)
                st.download_button("💾 Download CSV", df_stats.to_csv().encode(), "vi_statistics.csv","text/csv")

                st.markdown("### Distributions")
                fig_box = go.Figure()
                for vi_key in vis_to_stat:
                    entry = st.session_state["vi_cache"].get(vi_key)
                    if entry and entry[0]=="vi":
                        vi_map = entry[1]
                        data_flat = vi_map.flatten(); data_flat = data_flat[np.isfinite(data_flat)]
                        r = VI_LIBRARY.get(vi_key,{}).get("range",(-1,1))
                        data_flat = data_flat[(data_flat>=r[0])&(data_flat<=r[1])]
                        if len(data_flat)>100_000: data_flat=np.random.choice(data_flat,100_000,replace=False)
                        fig_box.add_trace(go.Violin(y=data_flat,name=vi_key,box_visible=True,meanline_visible=True))
                fig_box.update_layout(height=400,title="VI distributions",yaxis_title="Index value",margin=dict(t=50))
                st.plotly_chart(fig_box, use_container_width=True)

    st.divider()
    st.markdown("### Band Statistics")
    n_sample = st.slider("Bands to sample", 5, img.nbands, min(20,img.nbands))
    if st.button("📊 Compute Band Statistics"):
        with st.spinner("Computing…"):
            arr = read_as_array(img)
            bidxs = np.linspace(0,img.nbands-1,n_sample,dtype=int)
            band_stats = []
            for idx in bidxs:
                bd = arr[:,:,idx].flatten(); bd = bd[np.isfinite(bd)]
                band_stats.append({"Band":idx+1,"Wavelength (nm)":f"{wavelengths[idx]:.1f}",
                                   "Mean":float(np.mean(bd)),"Std":float(np.std(bd)),
                                   "Min":float(np.min(bd)),"Max":float(np.max(bd))})
            df_b = pd.DataFrame(band_stats).round(4)
            st.dataframe(df_b, use_container_width=True, hide_index=True)
            fig_m = px.line(x=[wavelengths[i] for i in bidxs],y=[r["Mean"] for r in band_stats],
                            title="Mean reflectance spectrum",
                            labels={"x":"Wavelength (nm)","y":"Mean reflectance"})
            fig_m.update_traces(line_color="#2d7a3a"); fig_m.update_layout(height=320)
            st.plotly_chart(fig_m, use_container_width=True)
            st.download_button("💾 Download CSV", df_b.to_csv(index=False).encode(), "band_statistics.csv","text/csv")
