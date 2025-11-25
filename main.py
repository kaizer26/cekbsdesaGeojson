# app.py — Full script (fixed, no bounds_intersects)
import streamlit as st
import geopandas as gpd
import pandas as pd
import math
import os
import tempfile
import fiona
from io import BytesIO

# POP UP
import streamlit as st

with st.popover("📘 Ringkasan Proses Analisis"):
    st.markdown("""
### 📝 **Ringkasan Proses Script Analisis BS – Desa – SLS**

Script ini merupakan aplikasi **berbasis Streamlit** yang digunakan untuk melakukan analisis spasial antara:

- **BS (Blok Sensus)**  
- **Desa/Kelurahan**
- **SLS/SUBSLS (Satuan Lingkungan Setempat / RT)**

Tujuan utama proses adalah untuk menghitung persentase overlay antara BS dan Desa, mengidentifikasi ketidaksesuaian encoding, serta mengambil daftar SLS yang berada di dalam masing-masing BS.

---

### 🔧 **1. Input & Preprocessing Data**

Pengguna mengunggah:
- BS (GeoJSON)
- Desa (GeoJSON)
- (Opsional) Final SLS (GeoJSON) — berisi SLS dan SubSLS

Setelah file diunggah:
##### **1.1 Fix Geometri**
- Memperbaiki polygon rusak dengan `buffer(0)`
- Memastikan ID tetap string (menghindari hilangnya nol di depan)

##### **1.2 Penyeragaman CRS**
Script otomatis menghitung EPSG UTM berbasis centroid BS. Tujuannya:
- Proses overlay dan perhitungan area harus menggunakan CRS meter (UTM)
- Output bisa dikembalikan ke EPSG:4326

---

### 📐 **2. Analisis Overlay BS – Desa**

#### **2.1 Hitung luas BS**
Untuk tiap BS:
- `area_bs` dihitung dalam satuan meter persegi
- `encoded_iddesa` diambil dari 10 digit pertama IDBS

#### **2.2 Intersection**
Dengan `gpd.overlay()`, script membuat potongan area BS yang berada dalam masing-masing desa:
- Hasilnya adalah layer *AB*: potongan BS × desa
- Setiap potongan punya luas `area_part`

#### **2.3 Hitung Persentase**
`pct_of_bs = area_part / area_bs`
Ini menunjukkan berapa persen wilayah BS berada di desa tertentu.

---

### 🎯 **3. Identifikasi Desa Signifikan (≥ 20%)**
Potongan desa disaring:
`pct_of_bs >= ambang batas (misal 20%)`

#### **3.1 Hitung jumlah desa signifikan per BS**
`n_iddesa_ge20`

#### **3.2 Detail desa signifikan**
Format:
`iddesa1:60%; iddesa2:40%`

### **3.3 Desa Dominan**
Desa dengan persentase terbesar (`pct_of_bs` tertinggi)..

### **3.4 Validasi Encoding**
Membandingkan:
- `encoded_iddesa` (yang seharusnya) 
- `dominant_iddesa` (hasil overlay)

Jika encoded_iddesa ≠ dominant_iddesa → `mismatch_dominant = True`

---

### 🏘️ **4. Ekstraksi Nama SLS berdasarkan Lokasi (Opsional)**

Jika final SLS diunggah:

#### **4.1 Fix geometri & CRS**
SLS juga diperbaiki dan diubah ke CRS yang sama dengan BS.

#### **4.2 Spatial Join**
Script melakukan join spasial per BS, mengambil semua SLS yang berada di dalam setiap BS.

#### **4.3 Penamaan lengkap SLS-SUBSLS**
Jika `kdsubsls != "00"`:  
`nmsls_with_subsls = nmsls + "[" + kdsubsls + "]"`

#### **4.4 Gabungkan menjadi satu kolom**
`all_nmsls_with_subsls = "SLS01; SLS02 [01]; SLS03; ..."`

---

### 📦 **5. Output bs_flag**

Berisi antara lain:

- IDBS + metadata  
- area_bs  
- encoded_iddesa  
- n_iddesa_ge20  
- detail_iddesa_ge20  
- dominant_iddesa  
- pct_dominant  
- pct_encoded  
- mismatch_dominant  
- all_nmsls_with_subsls  

Tersedia ekspor:
- CSV  
- Excel  
- GPKG dengan 3 layer: `parts_ge20`, `bs_flag`, `desa_final`

---

### 🎯 **6. Tujuan Utama Script**

- Validasi kecocokan BS terhadap desa
    - Memeriksa apakah BS berada pada desa yang benar berdasarkan overlay geometri.
- Deteksi BS bermasalah
    1. BS berada di lebih dari satu desa signifikan
    2. BS dominan tidak sama dengan encoded desa
- Memetakan BS → SLS/SubSLS 
    Mengambil semua SLS yang ada dalam batas BS.
- Menyediakan analisis untuk updating wilayah & validasi lapangan
    """)



# ---------- Helpers ----------
def fix(gdf, idcol=None):
    """
    Copy GeoDataFrame, ensure idcol is string (if provided and exists),
    and try to fix common geometry problems using buffer(0).
    """
    gdf = gdf.copy()
    # safe check: idcol must be a string and in columns
    if isinstance(idcol, str) and (idcol in gdf.columns):
        gdf[idcol] = gdf[idcol].astype(str)
    # attempt to fix invalid geometry (buffer(0))
    try:
        gdf["geometry"] = gdf.geometry.buffer(0)
    except Exception:
        # if buffer fails, just keep original geometry
        pass
    return gdf

def utm_epsg_from_lonlat(lon, lat):
    zone = int(math.floor((lon + 180) / 6) + 1)
    south = lat < 0
    return f"EPSG:{32700 + zone if south else 32600 + zone}"

def to_metric_crs(gdf):
    """
    Returns (gdf_in_metric, epsg_code) using UTM zone based on centroid.
    """
    # guard: if geometry empty, raise
    if gdf.empty:
        raise ValueError("GeoDataFrame kosong saat mencoba menentukan CRS metric.")
    c = gdf.unary_union.centroid
    epsg = utm_epsg_from_lonlat(c.x, c.y)
    return gdf.to_crs(epsg), epsg

def safe_drop_geometry_col(df):
    """
    Return df without geometry column if exists (for displays / Excel).
    """
    if "geometry" in df.columns:
        return df.drop(columns="geometry")
    return df

def bbox_intersects(g1, g2):
    """
    Manual bounding-box intersection test (works for all shapely versions).
    Returns True if bbox(g1) intersects bbox(g2).
    """
    minx1, miny1, maxx1, maxy1 = g1.bounds
    minx2, miny2, maxx2, maxy2 = g2.bounds
    return not (maxx1 < minx2 or maxx2 < minx1 or maxy1 < miny2 or maxy2 < miny1)

# ---------- Core processing: BS vs Desa ----------
def run_process(bs, desa, sls, IDBS="idbs", IDDESA="iddesa", THRESH=0.20):
    """
    Main pipeline:
      - bs, desa, sls: GeoDataFrames (sls can be None)
      - IDBS, IDDESA: column names (strings)
      - THRESH: fraction (e.g. 0.20)
    Returns:
      AB_ge: GeoDataFrame of significant parts (in metric CRS)
      bs_flag: GeoDataFrame BS with computed attributes (in metric CRS)
      counts: DataFrame summary per BS
      desa_m: desa reprojected to metric CRS (for export)
    """
    # 1) fix inputs
    bs = fix(bs, IDBS)
    desa = fix(desa, IDDESA)

    # 2) project to metric CRS (UTM) for area computations
    bs_m, epsg = to_metric_crs(bs)
    desa_m = desa.to_crs(epsg)

    # 3) compute area of BS
    bs_m["area_bs"] = bs_m.geometry.area
    # encoded iddesa (first 10 digits)
    bs_m["encoded_iddesa"] = bs_m[IDBS].astype(str).str.slice(0, 10)

    # 4) overlay intersection between BS and Desa
    AB = gpd.overlay(
        bs_m[[IDBS, "area_bs", "geometry"]],
        desa_m[[IDDESA, "geometry"]],
        how="intersection",
        keep_geom_type=True
    )

    # 5) clean AB: remove empty/None geometries and non-polygon types
    AB = AB[AB.geometry.notna() & (~AB.geometry.is_empty)]
    AB = AB[AB.geometry.type.isin(["Polygon", "MultiPolygon"])].copy()

    # 6) compute area of parts and percentage relative to BS
    AB["area_part"] = AB.geometry.area
    AB["pct_of_bs"] = AB["area_part"] / AB["area_bs"].replace({0: pd.NA})

    # 7) filter significant parts using THRESH
    AB_ge = AB[(AB["pct_of_bs"] >= THRESH) & (AB["area_part"] > 0)].copy()

    # 8) counts: number of unique iddesa per idbs among significant parts
    counts = (
        AB_ge.groupby(IDBS)[IDDESA]
        .nunique()
        .reset_index(name="n_iddesa_ge")
    )

    # 9) detail table: list "iddesa:XX.XX%" ordered desc
    def details(df):
        tmp = (
            df.assign(pct=(df["pct_of_bs"] * 100).round(2))
              .sort_values("pct_of_bs", ascending=False)
              .apply(lambda r: f"{r[IDDESA]}:{r['pct']}%", axis=1)
        )
        return "; ".join(tmp.tolist())

    if not AB_ge.empty:
        detail_tbl = AB_ge.groupby(IDBS).apply(details).reset_index(name="detail_iddesa_ge")
    else:
        detail_tbl = pd.DataFrame(columns=[IDBS, "detail_iddesa_ge"])

    # 10) dominant iddesa per idbs (largest pct_of_bs)
    if not AB.empty:
        idx = AB.groupby(IDBS)["pct_of_bs"].idxmax()
        dominant = (
            AB.loc[idx, [IDBS, IDDESA, "pct_of_bs"]]
            .rename(columns={IDDESA: "dominant_iddesa", "pct_of_bs": "pct_dominant"})
        )
    else:
        dominant = pd.DataFrame(columns=[IDBS, "dominant_iddesa", "pct_dominant"])

    # 11) pct_encoded: if the encoded_iddesa matches a polygon iddesa, what's the pct
    pct_encoded = (
        AB[[IDBS, IDDESA, "pct_of_bs"]]
        .merge(bs_m[[IDBS, "encoded_iddesa"]], on=IDBS)
        .query(f"{IDDESA} == encoded_iddesa")
        .rename(columns={"pct_of_bs": "pct_encoded"})
        [[IDBS, "pct_encoded"]]
    )

    # 12) Merge into bs_m -> bs_flag
    bs_flag = (
        bs_m.merge(counts, on=IDBS, how="left")
            .merge(detail_tbl, on=IDBS, how="left")
            .merge(dominant, on=IDBS, how="left")
            .merge(pct_encoded, on=IDBS, how="left")
    )

    # ensure defaults and nice numeric columns
    bs_flag["n_iddesa_ge"] = bs_flag["n_iddesa_ge"].fillna(0).astype(int)
    bs_flag["pct_encoded"] = (bs_flag["pct_encoded"] * 100).round(2).fillna(0)
    bs_flag["pct_dominant"] = (bs_flag["pct_dominant"] * 100).round(2).fillna(0)
    bs_flag["mismatch_dominant"] = (bs_flag["encoded_iddesa"] != bs_flag["dominant_iddesa"])

    # 13) If sls provided: spatial intersection of BS with each SLS
    if sls is None:
        bs_flag["all_nmsls_with_subsls"] = None
        bs_flag["nmsls_ge20"] = None
        bs_flag["idsls_ge20"] = None

    else:
        # prepare sls: fix geometry, ensure columns exist
        gdf_sls = fix(sls)
        try:
            gdf_sls = gdf_sls.to_crs(epsg)
        except:
            pass

        if "nmsls" not in gdf_sls.columns:
            gdf_sls["nmsls"] = None
        if "idsubsls" not in gdf_sls.columns:
            gdf_sls["idsubsls"] = None

        if "kdsubsls" not in gdf_sls.columns:
            gdf_sls["kdsubsls"] = gdf_sls["idsubsls"].astype(str).str[-2:].fillna("00")

        # label with subsls
        def label_row(r):
            n = r.get("nmsls", None)
            k = r.get("kdsubsls", None)
            if n is None or pd.isna(n):
                return None
            if k in [None, "00"] or pd.isna(k):
                return str(n)
            return f"{n} [{k}]"

        gdf_sls["nmsls_with_subsls"] = gdf_sls.apply(label_row, axis=1)

        # spatial index
        try:
            sindex = gdf_sls.sindex
        except:
            sindex = None

        # dictionary hasil
        idbs_all_labels = {}
        idbs_nmsls_ge20 = {}
        idbs_idsls_ge20 = {}

        # --- LOOP PER BS ---
        for idx_bs, bs_row in bs_m[[IDBS, "geometry", "area_bs"]].iterrows():
            idbs_val = bs_row[IDBS]
            bs_geom = bs_row["geometry"]
            bs_area = bs_row["area_bs"]

            if bs_geom is None or bs_geom.is_empty or bs_area <= 0:
                idbs_all_labels[idbs_val] = None
                idbs_nmsls_ge20[idbs_val] = None
                idbs_idsls_ge20[idbs_val] = None
                continue

            # kandidat SLS
            if sindex is not None:
                try:
                    candidate_idxs = list(sindex.intersection(bs_geom.bounds))
                except:
                    candidate_idxs = list(range(len(gdf_sls)))
            else:
                candidate_idxs = list(range(len(gdf_sls)))

            # collect
            label_list = []
            nmsls_set = set()
            idsls_set = set()

            for ci in candidate_idxs:
                sls_row = gdf_sls.iloc[ci]
                sls_geom = sls_row.geometry
                if sls_geom is None or sls_geom.is_empty:
                    continue

                # bbox check
                if not bbox_intersects(bs_geom, sls_geom):
                    continue

                try:
                    inter = bs_geom.intersection(sls_geom)
                except:
                    continue
                if inter.is_empty:
                    continue

                inter_area = inter.area
                pct = inter_area / bs_area

                if pct >= THRESH:
                    # label lengkap
                    lab = sls_row["nmsls_with_subsls"]
                    if lab:
                        label_list.append(str(lab))

                    # untuk nmsls_ge20
                    nmsls_set.add(sls_row["nmsls"])

                    # untuk idsls_ge20
                    idsls_set.add(str(sls_row["idsls"]))

            # simpan hasil
            idbs_all_labels[idbs_val] = "; ".join(sorted(set(label_list))) if label_list else None
            idbs_nmsls_ge20[idbs_val] = "; ".join(sorted(nmsls_set)) if nmsls_set else None
            idbs_idsls_ge20[idbs_val] = "; ".join(sorted(idsls_set)) if idsls_set else None

        # mapping ke bs_flag
        bs_flag["all_nmsls_with_subsls"] = bs_flag[IDBS].map(idbs_all_labels)
        bs_flag["nmsls_ge20"] = bs_flag[IDBS].map(idbs_nmsls_ge20)
        bs_flag["idsls_ge20"] = bs_flag[IDBS].map(idbs_idsls_ge20)

    # done
    return AB_ge, bs_flag, counts, desa_m

# ---------- Export multi-layer GPKG safely ----------
def export_gpkg(AB_ge, bs_flag, desa_m, out_crs="EPSG:4326"):
    """
    Write AB_ge (parts), bs_flag, and desa_m into a single GPKG file.
    Reproject outputs to out_crs for consistency.
    Uses tempfile.mkstemp + fiona.Env to avoid WinError 32 issues.
    """
    fd, path = tempfile.mkstemp(suffix=".gpkg")
    os.close(fd)

    # reproject to out_crs (if geometry exists)
    try:
        AB_out = AB_ge.to_crs(out_crs) if not AB_ge.empty else AB_ge
    except Exception:
        AB_out = AB_ge

    try:
        bs_out = bs_flag.to_crs(out_crs) if "geometry" in bs_flag.columns else bs_flag
    except Exception:
        bs_out = bs_flag

    try:
        desa_out = desa_m.to_crs(out_crs)
    except Exception:
        desa_out = desa_m

    with fiona.Env():
        if not AB_out.empty:
            AB_out.to_file(path, layer="parts_geTHRESH", driver="GPKG")
        # bs_out may include geometry
        bs_out.to_file(path, layer="bs_flag", driver="GPKG")
        desa_out.to_file(path, layer="final_desa", driver="GPKG")

    return path

# ---------- Streamlit UI ----------
st.set_page_config(page_title="BS vs Desa vs SLS", layout="centered")
st.markdown("### 🔎 Analisis Blok Sensus vs Desa (dengan optional Final SLS join by location)")

st.markdown("""
- Upload BS (GeoJSON/GPKG), Desa (GeoJSON/GPKG).  
- Optional: upload Final SLS (GeoJSON/GPKG).  
- THRESH = fraction (misal 0.20 = 20%) dipakai untuk:
  - memilih potongan DESA yang signifikan terhadap luas BS (parts >= THRESH)
  - memilih SLS yang berpotongan signifikan terhadap BS (intersection / BS area >= THRESH)
""")

col1, col2, col3 = st.columns([1,1,1])
with col1:
    f_bs = st.file_uploader("Upload BS (GeoJSON/GPKG/SHAPE)", type=["geojson","gpkg","shp"])
with col2:
    f_desa = st.file_uploader("Upload Desa (GeoJSON/GPKG/SHAPE)", type=["geojson","gpkg","shp"])
with col3:
    f_sls = st.file_uploader("Upload Final SLS (GeoJSON/GPKG/SHAPE) — optional", type=["geojson","gpkg","shp"])

IDBS = st.text_input("Kolom ID BS (nama kolom)", value="idbs")
IDDESA = st.text_input("Kolom ID Desa (nama kolom)", value="iddesa")
THRESH_PCT = st.slider("Threshold (%) - dipakai untuk parts & SLS join", 1, 50, 20)
THRESH = THRESH_PCT / 100.0

# Process button
if st.button("🚀 Proses Semua"):
    if (f_bs is None) or (f_desa is None):
        st.error("Harap upload file BS dan Desa terlebih dahulu.")
    else:
        with st.spinner("Memproses — ini bisa memakan waktu tergantung ukuran data..."):
            # read files with geopandas
            try:
                bs = gpd.read_file(f_bs)
                desa = gpd.read_file(f_desa)
            except Exception as e:
                st.error(f"Gagal membaca BS/Desa: {e}")
                st.stop()

            # read sls only if provided
            sls = None
            if f_sls is not None:
                try:
                    sls = gpd.read_file(f_sls)
                except Exception as e:
                    st.warning(f"Warning: Gagal membaca final_sls: {e}. Lanjut tanpa SLS.")
                    sls = None

            # run main pipeline
            try:
                AB_ge, bs_flag, counts, desa_m = run_process(bs, desa, sls, IDBS=IDBS, IDDESA=IDDESA, THRESH=THRESH)
            except Exception as e:
                st.error(f"Error saat pemrosesan: {e}")
                raise

            # store results in session state so downloads don't re-run heavy computations
            st.session_state["AB_ge"] = AB_ge
            st.session_state["bs_flag"] = bs_flag
            st.session_state["desa_m"] = desa_m
            st.session_state["counts"] = counts
            st.session_state["THRESH"] = THRESH

        st.success("Selesai memproses data. Hasil disimpan di session.")

# ---------- show results and downloads (if present) ----------
if "bs_flag" in st.session_state:
    bs_flag = st.session_state["bs_flag"]
    AB_ge = st.session_state["AB_ge"]
    desa_m = st.session_state["desa_m"]
    counts = st.session_state.get("counts", pd.DataFrame())

    # st.markdown("### 📍 Ringkasan (counts)")
    # st.dataframe(counts)

    # st.markdown("### 📋 BS flag (tabel)")
    display_tbl = safe_drop_geometry_col(bs_flag)
    # st.dataframe(display_tbl)

    # show problematic rows
    problematic = bs_flag[(bs_flag["n_iddesa_ge"] > 1) | (bs_flag["mismatch_dominant"])]
    if not problematic.empty:
        st.markdown("### ⚠️ BS bermasalah (n_iddesa_ge > 1 atau mismatch_dominant)")
        st.dataframe(safe_drop_geometry_col(problematic))

    # ===============================
    # Kolom yang akan diekspor
    # ===============================
    export_cols = [
        "idbs", "nmprov", "nmkab", "nmkec", "nmdesa", "area_bs", "encoded_iddesa",
        "n_iddesa_ge", "detail_iddesa_ge20", "dominant_iddesa",
        "pct_dominant", "pct_encoded", "mismatch_dominant",
        "all_nmsls_with_subsls", "nmsls_ge20", "idsls_ge20"
    ]
    

    # Pastikan kolom yang tidak ada tidak menyebabkan error
    export_cols = [c for c in export_cols if c in display_tbl.columns]

    export_df = display_tbl[export_cols].copy()

    # Tentukan prefix nama file
    try:
        # Jika idbs tersedia → pakai 4 digit pertama
        prefix = str(export_df["idbs"].iloc[0])[:4]
    except:
        # fallback ke kdprov + kdkab
        prefix = kdprov + kdkab

    # ===============================
    # DOWNLOAD CSV
    # ===============================
    csv_data = export_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        f"⬇ Download bs_flag_{prefix}.csv",
        csv_data,
        f"bs_flag_{prefix}.csv",
        "text/csv"
    )

    # ===============================
    # DOWNLOAD EXCEL
    # ===============================
    excel_buf = BytesIO()
    with pd.ExcelWriter(excel_buf, engine="xlsxwriter") as writer:
        export_df.to_excel(writer, index=False, sheet_name="bs_flag")

    st.download_button(
        f"⬇ Download bs_flag_{prefix}.xlsx",
        excel_buf.getvalue(),
        f"bs_flag_{prefix}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

    # Export GPKG (reproject to EPSG:4326 inside)
    gpkg_path = export_gpkg(AB_ge, bs_flag, desa_m, out_crs="EPSG:4326")
    with open(gpkg_path, "rb") as f:
        st.download_button(f"⬇ Download Output GPKG {prefix} (parts, bs_flag, desa)",
                           f,
                           f"bs_vs_desa_sls_output_{prefix}.gpkg",
                           mime="application/octet-stream")

    st.info("Semua layer di GPKG direproject ke EPSG:4326 (WGS84).")
