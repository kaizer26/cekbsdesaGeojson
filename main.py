import streamlit as st
import geopandas as gpd
import pandas as pd
import math
import tempfile
import os
import io
import requests
from datetime import datetime
from zoneinfo import ZoneInfo

# ==============================
# 0. Config dasar
# ==============================
IDBS = "idbs"
IDDESA = "iddesa"
IDSLS = "idsls"
IDSUBSLS = "idsubsls"

# ==============================
# 1. Fetch API + Cache
# ==============================
@st.cache_data(ttl=3600)
def fetch_bps_domains():
    """
    Ambil seluruh domain dari API BPS.
    Mengembalikan list of dict (domain entries) atau [] jika error.
    """
    url = "https://webapi.bps.go.id/v1/api/domain/type/all/key/687e204db62094de46edbcd7ed7cb204/"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        j = r.json()
        if "data" not in j or len(j["data"]) < 2:
            st.error("Respons API BPS tidak mengandung daftar domain (format tidak sesuai).")
            return []
        return j["data"][1]
    except Exception as e:
        st.error(f"Gagal memuat data API BPS: {e}")
        return []

# ==============================
# 2. Helpers provinsi/kabupaten
# ==============================
def get_provinces(domains):
    # Provinsi = domain_id berakhiran "00", kecuali "0000"
    return [d for d in domains if d.get("domain_id", "").endswith("00") and d.get("domain_id") != "0000"]

def get_kabupaten(domains, prov_id):
    # Kab/Kota = domain_id yang punya prefix prov (2 digit) dan tidak berakhiran 00
    if not prov_id:
        return []
    prefix = prov_id[:2]
    return [d for d in domains if d.get("domain_id", "").startswith(prefix) and not d.get("domain_id", "").endswith("00")]

# ==============================
# 3. Utilitas GIS kecil
# ==============================
def safe_groupby_apply(df, group_col, func, colname):
    """
    Mengembalikan DataFrame aman, meskipun df kosong.
    group_col = string nama kolom untuk groupby
    func      = lambda function
    colname   = nama kolom output setelah apply
    """
    if df.empty:
        return pd.DataFrame(columns=[group_col, colname])

    out = df.groupby(group_col).apply(func)
    return out.reset_index(name=colname)


def fix(gdf, idcol=None):
    gdf = gdf.copy()
    if idcol and idcol in gdf.columns:
        gdf[idcol] = gdf[idcol].astype(str)
    # perbaiki geometri sederhana
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

def bersihkan_kolom_gpkg(gdf):
    forbidden = ["fid", "ogr_fid", "ogc_fid", "rowid", "oid"]
    for f in forbidden:
        if f in gdf.columns:
            gdf = gdf.drop(columns=[f])
    return gdf


# ==============================
# 4. Inisialisasi session_state
# ==============================
if "provinsi" not in st.session_state:
    st.session_state["provinsi"] = "6300"
if "kabupaten" not in st.session_state:
    st.session_state["kabupaten"] = "6310"
if "run_analysis" not in st.session_state:
    st.session_state["run_analysis"] = False
if "analysis_done" not in st.session_state:
    st.session_state["analysis_done"] = False
# hasil disimpan di session_state: export_df, export_df1, gpkg_bytes, timestamp

# ==============================
# 5. Sidebar UI (semua kontrol di sini)
# ==============================
# st.sidebar.title("Kontrol Aplikasi")
st.sidebar.header("Unggah File GeoJSON")
bs_file = st.sidebar.file_uploader("Unggah File BS GeoJSON", type="geojson")
sls_file = st.sidebar.file_uploader("Unggah File SLS GeoJSON", type="geojson")
threshold = st.sidebar.slider("Threshold Overlap Signifikan (%)", 0, 100, 20)

# Ambil domains dari API (cached)
domains = fetch_bps_domains()

if not domains:
    st.sidebar.error("Daftar domain BPS gagal dimuat. Cek koneksi atau API.")
    # tetap lanjutkan UI tapi dropdown kosong
provinces = get_provinces(domains)

# Dropdown Provinsi (sidebar)
if provinces:
    prov_choices = provinces
    # tentukan index sesuai session_state jika ada
    prov_index = 0
    if st.session_state["provinsi"]:
        for i, p in enumerate(prov_choices):
            if p.get("domain_id") == st.session_state["provinsi"]:
                prov_index = i
                break

    selected_prov = st.sidebar.selectbox(
        "Pilih Provinsi",
        prov_choices,
        format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
        index=prov_index,
        key="provinsi_select"
    )
    # simpan kode provinsi (domain_id) ke session_state
    st.session_state["provinsi"] = selected_prov["domain_id"]
else:
    selected_prov = None

# Dropdown Kabupaten (sidebar) - dinamis sesuai provinsi
kabupaten_list = get_kabupaten(domains, st.session_state["provinsi"]) if selected_prov else []
if kabupaten_list:
    kab_index = 0
    if st.session_state["kabupaten"]:
        for i, k in enumerate(kabupaten_list):
            if k.get("domain_id") == st.session_state["kabupaten"]:
                kab_index = i
                break

    selected_kab = st.sidebar.selectbox(
        "Pilih Kabupaten/Kota",
        kabupaten_list,
        format_func=lambda x: f"{x['domain_id']} - {x['domain_name']}",
        index=kab_index,
        key="kabupaten_select"
    )
    st.session_state["kabupaten"] = selected_kab["domain_id"]
else:
    selected_kab = None

# Tombol Jalankan di sidebar
if st.sidebar.button("Jalankan"):
    # reset hasil lama dan set trigger run
    st.session_state["run_analysis"] = True
    st.session_state["analysis_done"] = False
    for k in ["export_df", "export_df1", "gpkg_bytes", "timestamp"]:
        st.session_state.pop(k, None)

# ==============================
# 6. Judul aplikasi (main area)
# ==============================
st.title("Pengecekan Desa Blok Sensus dan SLS")

# Informasi singkat di main
st.markdown("Gunakan sidebar untuk memilih provinsi/kabupaten, mengunggah file, dan menjalankan analisis.")

# ==============================
# 7. Proses Analisis (hanya dijalankan sekali per klik)
# ==============================
if st.session_state.get("run_analysis") and not st.session_state.get("analysis_done"):
    # validasi input file
    if not (bs_file and sls_file):
        st.warning("Unggah file GeoJSON untuk BS dan SLS terlebih dahulu (di sidebar).")
    else:
        # set timestamp konsisten untuk semua output
        st.session_state["timestamp"] = datetime.now().strftime("%Y%m%d%H%M%S")

        # --- Load Data ---
        try:
            bs = fix(gpd.read_file(bs_file), IDBS)
            sls = fix(gpd.read_file(sls_file), IDSUBSLS)
        except Exception as e:
            st.error(f"Gagal membaca file GeoJSON: {e}")
            st.session_state["run_analysis"] = False

        else:
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
            ABC_grouped = ABC.groupby([IDBS, "area_bs", "iddesa"], as_index=False).agg({"area_part": "sum"})
            ABC_grouped1 = ABC.groupby([IDBS, "area_bs", "idsubsls", "nmsls"], as_index=False).agg({"area_part": "sum"})
            ABC_grouped2 = ABC.groupby([IDBS, "area_bs", "idsls", "nmsls"], as_index=False).agg({"area_part": "sum"})

            ABC_grouped["pct_of_bs"] = ABC_grouped["area_part"] / ABC_grouped["area_bs"]
            ABC_grouped1["pct_of_bs"] = ABC_grouped1["area_part"] / ABC_grouped1["area_bs"]
            ABC_grouped2["pct_of_bs"] = ABC_grouped2["area_part"] / ABC_grouped2["area_bs"]

            # --- threshold filter ---
            ABC_geX    = ABC_grouped[ABC_grouped["pct_of_bs"] >= (threshold / 100)].copy()
            ABC_geX_g1 = ABC_grouped1[ABC_grouped1["pct_of_bs"] >= (threshold / 100)].copy()
            ABC_geX_g2 = ABC_grouped2[ABC_grouped2["pct_of_bs"] >= (threshold / 100)].copy()
            
            # --- Count iddesa signifikan ---
            if ABC_geX.empty:
                countsC = pd.DataFrame(columns=[IDBS, "n_iddesa_ge20"])
            else:
                countsC = ABC_geX.groupby(IDBS)[IDDESA] \
                    .nunique().reset_index(name="n_iddesa_ge20")

            # --- Detail persen per iddesa ---
            detail_tblC = safe_groupby_apply(
                ABC_geX.assign(pct=lambda d: (d["pct_of_bs"] * 100).round(2))
                    .sort_values([IDBS, "pct"], ascending=[True, False]),
                IDBS,
                lambda d: "; ".join(f"{r[IDDESA]}:{r['pct']}%" for _, r in d.iterrows()),
                "detail_iddesa_ge20"
            )

            # --- nmsls_g1 ---
            nmsls1 = safe_groupby_apply(
                ABC_geX_g1.assign(pct=lambda d: (d["pct_of_bs"] * 100).round(2))
                        .sort_values(["idbs", "pct"], ascending=[True, False]),
                IDBS,
                lambda d: " , ".join(
                    f"{r['nmsls']} - [{r[IDSUBSLS][-2:]}]" if r[IDSUBSLS][-2:] != '00'
                    else f"{r['nmsls']}"
                    for _, r in d.iterrows()
                ),
                "nmsls_g1"
            )

            # --- kdsls_g1 ---
            kdsls1 = safe_groupby_apply(
                ABC_geX_g1.assign(pct=lambda d: (d["pct_of_bs"] * 100).round(2))
                        .sort_values(["idbs", IDSUBSLS]),
                IDBS,
                lambda d: " , ".join(f"{r[IDSUBSLS][-6:]}" for _, r in d.iterrows()),
                "kdsls_g1"
            )

            # --- kdsls_g2 ---
            kdsls2 = safe_groupby_apply(
                ABC_geX_g2.assign(pct=lambda d: (d["pct_of_bs"] * 100).round(2))
                        .sort_values(["idbs", IDSLS]),
                IDBS,
                lambda d: " , ".join(f"{r[IDSLS][-4:]}" for _, r in d.iterrows()),
                "kdsls_g2"
            )

            # --- nmsls_g2 ---
            nmsls2 = safe_groupby_apply(
                ABC_geX_g2.assign(pct=lambda d: (d["pct_of_bs"] * 100).round(2))
                        .sort_values(["idbs", "pct"], ascending=[True, False]),
                IDBS,
                lambda d: " , ".join(f"{r['nmsls']}" for _, r in d.iterrows()),
                "nmsls_g2"
            )

            # --- merge final ---
            result = (
                nmsls1.merge(kdsls1, on="idbs", how="outer")
                    .merge(kdsls2, on="idbs", how="outer")
                    .merge(nmsls2, on="idbs", how="outer")
            )

            # --- Dominant iddesa ---
            idx = ABC_grouped.groupby(IDBS)["pct_of_bs"].idxmax()
            dominant = ABC_grouped.loc[idx, [IDBS, IDDESA, "pct_of_bs"]].rename(
                columns={IDDESA: "dominant_iddesa", "pct_of_bs": "pct_dominant"}
            )

            # --- encoded iddesa %
            pct_encoded = ABC_grouped[[IDBS, IDDESA, "pct_of_bs"]].merge(
                bs_m[[IDBS, "encoded_iddesa"]], on=IDBS
            ).query(f"{IDDESA} == encoded_iddesa").rename(columns={"pct_of_bs": "pct_encoded"})[[IDBS, "pct_encoded"]]

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

            sls_4326 = bersihkan_kolom_gpkg(sls_4326)
            bs_4326 = bersihkan_kolom_gpkg(bs_4326)
            bs_flag_4326 = bersihkan_kolom_gpkg(bs_flag_4326)
            
            # --- Layer hasil analisis (GeoDataFrame) ---
            hasil_analisis_4326 = hasil_analisis.to_crs(4326)
            hasil_analisis_4326 = bersihkan_kolom_gpkg(hasil_analisis_4326)
            
            # Nama layer mengikuti nama file upload 
            bs_layer_name = bs_file.name.replace(".geojson", "") 
            sls_layer_name = sls_file.name.replace(".geojson", "")
            
            with tempfile.TemporaryDirectory() as tmpdir:
                gpkg_path = os.path.join(tmpdir, "output.gpkg")

                hasil_analisis_4326.to_file(gpkg_path, layer="bs_flag", driver="GPKG")
                bs_4326.to_file(gpkg_path, layer=bs_layer_name, driver="GPKG", mode="a")
                sls_4326.to_file(gpkg_path, layer=sls_layer_name, driver="GPKG", mode="a")
                # hasil_pemeriksaan_4326.to_file(gpkg_path, layer="hasil_pemeriksaan_bs_desa", driver="GPKG", mode="a")

                with open(gpkg_path, "rb") as f:
                    gpkg_bytes = f.read()

            # Simpan hasil ke session_state
            st.session_state["export_df"] = export_df
            st.session_state["export_df1"] = export_df1
            st.session_state["gpkg_bytes"] = gpkg_bytes
            st.session_state["analysis_done"] = True
            st.session_state["run_analysis"] = False

# ==============================
# 8. Tampilkan hasil & tombol download
# ==============================
if st.session_state.get("analysis_done"):
    st.markdown(f"#### Jumlah IDBS yang Bermasalah: {len(st.session_state['export_df'])}")
    st.dataframe(st.session_state["export_df"])

    # gunakan timestamp yang tersimpan
    # datetime.now(ZoneInfo("Asia/Jakarta")).strftime("%Y%m%d%H%M%S")
    tanggal_ymdhms = datetime.now(ZoneInfo("Asia/Jakarta")).strftime("%Y%m%d%H%M%S")
    kab = st.session_state.get("kabupaten", "")

    # Download CSV
    st.download_button(
        f"Unduh Hasil CSV ({tanggal_ymdhms})",
        data=st.session_state["export_df1"].to_csv(index=False),
        file_name=f"{tanggal_ymdhms}_Hasil_pengecekan_{kab}.csv",
        mime="text/csv"
    )

    # Download Excel
    excel_data = io.BytesIO()
    with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
        st.session_state["export_df1"].to_excel(writer, index=False)
    excel_data.seek(0)

    st.download_button(
        f"Unduh Hasil Excel ({tanggal_ymdhms})",
        data=excel_data,
        file_name=f"{tanggal_ymdhms}_Hasil_pengecekan_{kab}.xlsx",
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
