import streamlit as st
import geopandas as gpd
import pandas as pd
import math
import tempfile
import os
import io
import requests
from datetime import datetime

# ==============================
# 1. Fetch API + Cache
# ==============================
@st.cache_data(ttl=3600)  # cache 1 jam
def fetch_bps_domains():
    url = "https://webapi.bps.go.id/v1/api/domain/type/all/key/687e204db62094de46edbcd7ed7cb204/"
    r = requests.get(url, timeout=10)

    if r.status_code != 200:
        st.error("Gagal memuat data API BPS")
        return []

    j = r.json()

    # data ada pada index [1]
    if "data" not in j:
        st.error("Format API berubah")
        return []

    return j["data"][1]


# ==============================
# 2. Helper: Ambil provinsi & kabupaten
# ==============================
def get_provinces(domains):
    return [
        d for d in domains
        if d["domain_id"].endswith("00") and d["domain_id"] != "0000"
    ]


def get_kabupaten(domains, prov_id):
    prefix = prov_id[:2]
    return [
        d for d in domains
        if d["domain_id"].startswith(prefix)
        and not d["domain_id"].endswith("00")
    ]


# ==============================
# 3. Init session_state
# ==============================
if "provinsi" not in st.session_state:
    st.session_state["provinsi"] = None

if "kabupaten" not in st.session_state:
    st.session_state["kabupaten"] = None


# ==============================
# 4. UI Dropdown
# ==============================

domains = fetch_bps_domains()
provinces = get_provinces(domains)

# Dropdown Provinsi
selected_prov = st.selectbox(
    "Pilih Provinsi",
    provinces,
    format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
    index=(
        next((i for i, p in enumerate(provinces)
              if st.session_state["provinsi"] and p["domain_id"] == st.session_state["provinsi"]), 0)
    ),
    key="provinsi_select",
)

# Simpan ke session_state
st.session_state["provinsi"] = selected_prov["domain_id"]


# Dropdown Kabupaten dinamis
kabupaten = get_kabupaten(domains, st.session_state["provinsi"])

selected_kab = st.selectbox(
    "Pilih Kabupaten/Kota",
    kabupaten,
    format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
    index=(
        next((i for i, k in enumerate(kabupaten)
              if st.session_state["kabupaten"] and k["domain_id"] == st.session_state["kabupaten"]), 0)
    ),
    key="kabupaten_select",
)

# Simpan ke session_state
st.session_state["kabupaten"] = selected_kab["domain_id"]

# ====== KONFIGURASI ======
folder_path = 'G:\\IPDS\\STREAMLIT\\'
BS_PATH = f"{folder_path}Final_BS_202416310sp2020.geojson"
SLS_PATH = f"{folder_path}peta_sls_202516310.geojson"
IDBS = "idbs"
IDDESA = "iddesa"
IDSLS = "idsls"
IDSUBSLS = "idsubsls"
THRESH = 0.20
OUT_GPKG = f"{folder_path}bs_vs_sls.gpkg"
domains = load_domains()

# ====== FUNGSI ======
def fix(gdf, idcol=None):
    gdf = gdf.copy()
    if idcol and idcol in gdf.columns:
        gdf[idcol] = gdf[idcol].astype(str)
    gdf["geometry"] = gdf.buffer(0)
    return gdf

def utm_epsg_from_lonlat(lon, lat):
    zone = int(math.floor((lon + 180) / 6) + 1)
    south = lat < 0
    return f"EPSG:{32700 + zone if south else 32600 + zone}"

def to_metric_crs(gdf):
    c = gdf.unary_union.centroid
    epsg = utm_epsg_from_lonlat(c.x, c.y)
    return gdf.to_crs(epsg), epsg

# ====== STREAMLIT UI ======
st.title("Pengecekan Desa Blok Sensus dan SLS")
st.sidebar.header("Unggah File GeoJSON")

bs_file = st.sidebar.file_uploader("Unggah File BS GeoJSON", type="geojson")
sls_file = st.sidebar.file_uploader("Unggah File SLS GeoJSON", type="geojson")
threshold = st.sidebar.slider("Threshold Overlap Signifikan (%)", 0, 100, 20)

# Ambil daftar provinsi
provinces = get_provinces(domains)

# Dropdown Provinsi
selected_prov = st.selectbox(
    "Pilih Provinsi",
    provinces,
    format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
    key="provinsi"
)

# Ambil kabupaten sesuai provinsi terpilih
kabupaten = get_kabupaten(domains, selected_prov["domain_id"])

# Dropdown Kabupaten otomatis
selected_kab = st.selectbox(
    "Pilih Kabupaten/Kota",
    kabupaten,
    format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
    key="kabupaten"
)

# Tombol run → hanya set flag
if st.sidebar.button("Jalankan Analisis"):
    st.session_state["run_analysis"] = True
    st.session_state["analysis_done"] = False   # reset agar bisa run ulang
    # hapus hasil lama
    for k in ["export_df", "export_df1", "gpkg_bytes"]:
        st.session_state.pop(k, None)


# ======================================================================
# ====== PROSES ANALISIS (hanya dijalankan sekali) ======
# ======================================================================
if (
    st.session_state.get("run_analysis") 
    and not st.session_state.get("analysis_done") 
    and bs_file 
    and sls_file
):

    # --- Load Data ---
    bs = fix(gpd.read_file(bs_file), IDBS)
    sls = fix(gpd.read_file(sls_file), IDSUBSLS)

    bs_m, epsg = to_metric_crs(bs)
    sls_m = sls.to_crs(epsg)

    # --- Luas BS ---
    bs_m["area_bs"] = bs_m.geometry.area
    bs_m["encoded_iddesa"] = bs_m[IDBS].str.slice(0, 10)
    sls_m["iddesa"] = sls_m[IDSLS].str.slice(0, 10)

    # --- Overlay ---
    ABC = gpd.overlay(
        bs_m[[IDBS, "area_bs", "geometry"]],
        sls_m[[IDSUBSLS, IDSLS, IDDESA, "nmsls", "geometry"]],
        how="intersection",
        keep_geom_type=True
    )
    ABC["area_part"] = ABC.geometry.area

    # --- Group by ---
    ABC_grouped = ABC.groupby([IDBS, "area_bs", "iddesa"], as_index=False).agg({
        "area_part": "sum",
    })
    ABC_grouped1 = ABC.groupby([IDBS, "area_bs", "idsubsls", "nmsls"], as_index=False).agg({
        "area_part": "sum",
    })
    ABC_grouped2 = ABC.groupby([IDBS, "area_bs", "idsls", "nmsls"], as_index=False).agg({
        "area_part": "sum",
    })

    ABC_grouped["pct_of_bs"] = ABC_grouped["area_part"] / ABC_grouped["area_bs"]
    ABC_grouped1["pct_of_bs"] = ABC_grouped1["area_part"] / ABC_grouped1["area_bs"]
    ABC_grouped2["pct_of_bs"] = ABC_grouped2["area_part"] / ABC_grouped2["area_bs"]

    # --- >20% ---
    ABC_ge20 = ABC_grouped[ABC_grouped["pct_of_bs"] >= threshold / 100].copy()
    ABC_ge20_g1 = ABC_grouped1[ABC_grouped1["pct_of_bs"] >= threshold / 100].copy()
    ABC_ge20_g2 = ABC_grouped2[ABC_grouped2["pct_of_bs"] >= threshold / 100].copy()

    # --- Count iddesa signifikan ---
    countsC = ABC_ge20.groupby(IDBS)[IDDESA].nunique().reset_index(name="n_iddesa_ge20")

    # --- Detail persentase iddesa signifikan ---
    detail_tblC = ABC_ge20.assign(
        pct=lambda d: (d["pct_of_bs"] * 100).round(2)
    ).sort_values(
        [IDBS, "pct"], ascending=[True, False]
    ).groupby(
        IDBS
    ).apply(
        lambda d: "; ".join(f"{r[IDDESA]}:{r['pct']}%" for _, r in d.iterrows())
    ).reset_index(
        name="detail_iddesa_ge20"
    )

    # --- nmsls & kdsls ---
    nmsls1 = ABC_ge20_g1.assign(
        pct=lambda d: (d["pct_of_bs"] * 100).round(2)
    ).sort_values(
        ["idbs", "pct"], ascending=[True, False]
    ).groupby(IDBS).apply(
        lambda d: " , ".join(
            f"{r['nmsls']} - [{r[IDSUBSLS][-2:]}]" if r[IDSUBSLS][-2:] != '00' else f"{r['nmsls']}"
            for _, r in d.iterrows()
        )
    ).reset_index(name="nmsls_g1")

    kdsls1 = ABC_ge20_g1.assign(
        pct=lambda d: (d["pct_of_bs"] * 100).round(2)
    ).sort_values(
        ["idbs", IDSUBSLS]
    ).groupby(IDBS).apply(
        lambda d: " , ".join(f"{r[IDSUBSLS][-6:]}" for _, r in d.iterrows())
    ).reset_index(name="kdsls_g1")

    kdsls2 = ABC_ge20_g2.assign(
        pct=lambda d: (d["pct_of_bs"] * 100).round(2)
    ).sort_values(
        ["idbs", IDSLS]
    ).groupby(IDBS).apply(
        lambda d: " , ".join(f"{r[IDSLS][-4:]}" for _, r in d.iterrows())
    ).reset_index(name="kdsls_g2")

    nmsls2 = ABC_ge20_g2.assign(
        pct=lambda d: (d["pct_of_bs"] * 100).round(2)
    ).sort_values(
        ["idbs", "pct"], ascending=[True, False]
    ).groupby(IDBS).apply(
        lambda d: " , ".join(f"{r['nmsls']}" for _, r in d.iterrows())
    ).reset_index(name="nmsls_g2")

    result = nmsls1.merge(kdsls1, on="idbs", how="outer") \
                   .merge(kdsls2, on="idbs", how="outer") \
                   .merge(nmsls2, on="idbs", how="outer")

    # --- Dominant iddesa ---
    idx = ABC_grouped.groupby(IDBS)["pct_of_bs"].idxmax()
    dominant = ABC_grouped.loc[idx, [IDBS, IDDESA, "pct_of_bs"]].rename(
        columns={IDDESA: "dominant_iddesa", "pct_of_bs": "pct_dominant"}
    )

    # --- encoded iddesa %
    pct_encoded = (
        ABC_grouped[[IDBS, IDDESA, "pct_of_bs"]]
        .merge(bs_m[[IDBS, "encoded_iddesa"]], on=IDBS)
        .query(f"{IDDESA} == encoded_iddesa")
        .rename(columns={"pct_of_bs": "pct_encoded"})
        [[IDBS, "pct_encoded"]]
    )

    # --- Gabungkan ke BS Flag ---
    bs_flag = (
        bs_m.merge(countsC, on=IDBS, how="left")
            .merge(detail_tblC, on=IDBS, how="left")
            .merge(dominant, on=IDBS, how="left")
            .merge(pct_encoded, on=IDBS, how="left")
            .merge(result, on=IDBS, how="left")
    )

    bs_flag["n_iddesa_ge20"] = bs_flag["n_iddesa_ge20"].fillna(0).astype(int)
    bs_flag["pct_encoded"] = (bs_flag["pct_encoded"] * 100).round(2).fillna(0)
    bs_flag["pct_dominant"] = (bs_flag["pct_dominant"] * 100).round(2).fillna(0)
    bs_flag["mismatch_dominant"] = bs_flag["encoded_iddesa"] != bs_flag["dominant_iddesa"]

    # Filter hasil analisis final
    hasil_analisis = bs_flag[
        bs_flag["mismatch_dominant"] | (bs_flag["n_iddesa_ge20"] > 1)
    ]

    export_cols = [
        "idbs", "nmprov", "nmkab", "nmkec", "nmdesa", "area_bs", "encoded_iddesa",
        "n_iddesa_ge20", "detail_iddesa_ge20", "dominant_iddesa", "pct_dominant",
        "pct_encoded", "nmsls_g1", "kdsls_g1", "kdsls_g2", "nmsls_g2", "mismatch_dominant"
    ]

    export_cols = [c for c in export_cols if c in hasil_analisis.columns]

    export_df = hasil_analisis[export_cols].copy()
    export_df1 = bs_flag[export_cols].copy()

    # --- GPKG Output ---
    bs_4326 = bs.to_crs(4326)
    sls_4326 = sls.to_crs(4326)
    bs_flag_4326 = bs_flag.to_crs(4326)
    # Nama layer mengikuti nama file upload 
    bs_layer_name = bs_file.name.replace(".geojson", "") 
    sls_layer_name = sls_file.name.replace(".geojson", "")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        gpkg_path = os.path.join(tmpdir, "output.gpkg")

        bs_flag_4326.to_file(gpkg_path, layer="bs_flag", driver="GPKG")
        bs_4326.to_file(gpkg_path, layer=bs_layer_name, driver="GPKG", mode="a")
        sls_4326.to_file(gpkg_path, layer=sls_layer_name, driver="GPKG", mode="a")

        with open(gpkg_path, "rb") as f:
            gpkg_bytes = f.read()

    # Simpan hasil ke session_state
    st.session_state["export_df"] = export_df
    st.session_state["export_df1"] = export_df1
    st.session_state["gpkg_bytes"] = gpkg_bytes
    st.session_state["analysis_done"] = True

else: st.warning("Unggah file GeoJSON untuk BS dan SLS terlebih dahulu.")

# ======================================================================
# ====== TAMPILKAN OUTPUT BILA ANALISIS SUDAH SELESAI ======
# ======================================================================
if st.session_state.get("analysis_done"):

    st.markdown(f"#### Jumlah IDBS Bermasalah: {len(st.session_state['export_df'])}")
    st.dataframe(st.session_state["export_df"])
    tanggal_ymdhms = datetime.now().strftime("%Y%m%d%H%M%S")
    kab = st.session_state["kabupaten"]

    # Download CSV
    st.download_button(
        "Unduh Hasil CSV",
        data=st.session_state["export_df1"].to_csv(index=False),
        file_name=f"{tanggal_ymdhms}_Hasil pengecekan desa blok sensus dan sls_{kab}.csv",
        mime="text/csv"
    )

    # Download Excel
    excel_data = io.BytesIO()
    with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
        st.session_state["export_df1"].to_excel(writer, index=False)
    excel_data.seek(0)

    st.download_button(
        "Unduh Hasil Excel",
        data=excel_data,
        file_name=f"{tanggal_ymdhms}_Hasil pengecekan desa blok sensus dan sls_{kab}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

    # Download GPKG
    st.download_button(
        "Unduh GeoPackage (EPSG:4326)",
        data=st.session_state["gpkg_bytes"],
        file_name=f"{tanggal_ymdhms}_bs_vs_sls_epsg4326_{kab}.gpkg",
        mime="application/geopackage+sqlite3"
    )

    st.success("Analisis selesai!")
    st.info("Gunakan bs_flag hanya sebagai alat bantu pengecekan, bukan untuk perbaikan BS.")
