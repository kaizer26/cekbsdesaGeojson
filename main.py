import streamlit as st
import geopandas as gpd
import pandas as pd
import math
import tempfile
import os


# ====== KONFIGURASI ======
folder_path = 'G:\\IPDS\\STREAMLIT\\'
BS_PATH = f"{folder_path}Final_BS_202416310sp2020.geojson"
SLS_PATH = f"{folder_path}peta_sls_202516310.geojson"
IDBS = "idbs"       # 14 digit
IDDESA = "iddesa"   # 12 digit
IDSLS = "idsls"     # 14 digit
IDSUBSLS = "idsubsls" # 16 digit
THRESH = 0.20       # 20%
OUT_GPKG = f"{folder_path}bs_vs_sls.gpkg"

# ====== FUNGSI PEMBANTU ======
def fix(gdf, idcol=None):
    gdf = gdf.copy()
    if idcol and idcol in gdf.columns:
        gdf[idcol] = gdf[idcol].astype(str)
    gdf["geometry"] = gdf.buffer(0)  # Perbaiki self-intersections, dll.
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
st.title("Analisis BS vs SLS")
st.sidebar.header("Unggah File GeoJSON")

# Pengunggah file untuk GeoJSON
bs_file = st.sidebar.file_uploader("Unggah File BS GeoJSON", type="geojson")
sls_file = st.sidebar.file_uploader("Unggah File SLS GeoJSON", type="geojson")

# Slider untuk Threshold
threshold = st.sidebar.slider("Threshold Overlap Signifikan (%)", 0, 100, 20)
run_analysis = st.sidebar.button("Jalankan")


# Memuat dan mempersiapkan data
if run_analysis:
    st.session_state["run"] = True
    
if st.session_state.get("run") and bs_file and sls_file:
    if "export_df" not in st.session_state:
        # Membaca file GeoJSON yang diunggah
        bs = fix(gpd.read_file(bs_file), IDBS)
        sls = fix(gpd.read_file(sls_file), IDSUBSLS)
        bs_m, epsg = to_metric_crs(bs)
        sls_m = sls.to_crs(epsg)
    
        # Menghitung Luas BS
        bs_m["area_bs"] = bs_m.geometry.area
        bs_m["encoded_iddesa"] = bs_m[IDBS].str.slice(0, 10)  # 12 digit pertama
        sls_m["iddesa"] = sls_m[IDSLS].str.slice(0, 10)  # 12 digit pertama
    
        # Intersection (tumpang tindih) antara BS dan SLS
        ABC = gpd.overlay(bs_m[[IDBS, "area_bs", "geometry"]], 
                          sls_m[[IDSUBSLS, IDSLS, IDDESA, "nmsls", "geometry"]],
                          how="intersection", keep_geom_type=True)
        ABC["area_part"] = ABC.geometry.area
    
        # Groupby dan agregasi
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
    
        # Potongan signifikan (>=20%)
        ABC_ge20 = ABC_grouped[ABC_grouped["pct_of_bs"] >= threshold / 100].copy()
        ABC_ge20_g1 = ABC_grouped1[ABC_grouped1["pct_of_bs"] >= threshold / 100].copy()
        ABC_ge20_g2 = ABC_grouped2[ABC_grouped2["pct_of_bs"] >= threshold / 100].copy()
    
        # Menghitung jumlah iddesa signifikan per idbs
        countsC = ABC_ge20.groupby(IDBS)[IDDESA].nunique().reset_index(name="n_iddesa_ge20")
    
        # Detail persentase iddesa signifikan per idbs
        detail_tblC = ABC_ge20.assign(pct=lambda d: (d["pct_of_bs"]*100).round(2)) \
            .sort_values([IDBS, "pct"], ascending=[True, False]) \
            .groupby(IDBS) \
            .apply(lambda d: "; ".join(f"{r[IDDESA]}:{r['pct']}%" for _, r in d.iterrows())) \
            .reset_index(name="detail_iddesa_ge20")
    
        # Menyiapkan nmsls dan kdsls berdasarkan aturan Anda
        nmsls1 = (ABC_ge20_g1
                      .assign(pct=lambda d: (d["pct_of_bs"]*100).round(2))
                      .sort_values(["idbs", "pct"], ascending=[True, False])
                      .groupby(IDBS)
                      .apply(lambda d: " , ".join(
                          f"{r['nmsls']} - [{r[IDSUBSLS][-2:]}]" if r[IDSUBSLS][-2:] != '00' else f"{r['nmsls']}"
                          for _, r in d.iterrows()))
                      .reset_index(name="nmsls_g1"))
    
        kdsls1 = (ABC_ge20_g1
                      .assign(pct=lambda d: (d["pct_of_bs"]*100).round(2))
                      .sort_values(["idbs", IDSUBSLS], ascending=[True, False])
                      .groupby(IDBS)
                      .apply(lambda d: " , ".join(f"{r[IDSUBSLS][-6:]}" for _, r in d.iterrows()))
                      .reset_index(name="kdsls_g1"))
    
        kdsls2 = (ABC_ge20_g2
                      .assign(pct=lambda d: (d["pct_of_bs"]*100).round(2))
                      .sort_values(["idbs", IDSLS], ascending=[True, False])
                      .groupby(IDBS)
                      .apply(lambda d: " , ".join(f"{r[IDSLS][-4:]}" for _, r in d.iterrows()))
                      .reset_index(name="kdsls_g2"))
    
        nmsls2 = (ABC_ge20_g2
                      .assign(pct=lambda d: (d["pct_of_bs"]*100).round(2))
                      .sort_values(["idbs", "pct"], ascending=[True, False])
                      .groupby(IDBS)
                      .apply(lambda d: " , ".join(f"{r['nmsls']}" for _, r in d.iterrows()))
                      .reset_index(name="nmsls_g2"))
    
        # Gabungkan tabel hasil
        result = nmsls1 \
            .merge(kdsls1, on="idbs", how="outer") \
            .merge(kdsls2, on="idbs", how="outer") \
            .merge(nmsls2, on="idbs", how="outer")
    
        # Dominant iddesa (persen terbesar per idbs)
        idx = ABC_grouped.groupby(IDBS)["pct_of_bs"].idxmax()
        dominant = (ABC_grouped.loc[idx, [IDBS, IDDESA, "pct_of_bs"]]
                    .rename(columns={IDDESA: "dominant_iddesa", "pct_of_bs": "pct_dominant"}))
    
        # Persentase untuk encoded_iddesa
        pct_encoded = (ABC_grouped[[IDBS, IDDESA, "pct_of_bs"]]
                       .merge(bs_m[[IDBS, "encoded_iddesa"]], on=IDBS)
                       .query(f"{IDDESA} == encoded_iddesa")
                       .rename(columns={"pct_of_bs": "pct_encoded"})
                       [[IDBS, "pct_encoded"]])
    
        # Gabungkan ke BS
        bs_flag = (bs_m.merge(countsC, on=IDBS, how="left")
                       .merge(detail_tblC, on=IDBS, how="left")
                       .merge(dominant, on=IDBS, how="left")
                       .merge(pct_encoded, on=IDBS, how="left")
                       .merge(result, on=IDBS, how="left"))
    
        # Perbaiki kolom NaN dan formatkan persentase
        bs_flag["n_iddesa_ge20"] = bs_flag["n_iddesa_ge20"].fillna(0).astype(int)
        bs_flag["pct_encoded"] = (bs_flag["pct_encoded"]*100).round(2).fillna(0)
        bs_flag["pct_dominant"] = (bs_flag["pct_dominant"]*100).round(2).fillna(0)
        bs_flag["mismatch_dominant"] = (bs_flag["encoded_iddesa"] != bs_flag["dominant_iddesa"])
    
        # Tampilkan hasil di Streamlit
        hasil_analisis = bs_flag[bs_flag["mismatch_dominant"]==True]
    
        export_cols = [
            "idbs", "nmprov", "nmkab", "nmkec", "nmdesa", "area_bs", "encoded_iddesa", 
            "n_iddesa_ge20", "detail_iddesa_ge20", "dominant_iddesa", "pct_dominant", 
            "pct_encoded", "nmsls_g1", "kdsls_g1", "kdsls_g2", "nmsls_g2", "mismatch_dominant"
        ]
        
    
        # Pastikan kolom yang tidak ada tidak menyebabkan error
        export_cols = [c for c in export_cols if c in hasil_analisis.columns]
    
        export_df = hasil_analisis[export_cols].copy()
        export_df1 = bs_flag[export_cols].copy()

        # ------ SIMPAN 3 LAYER KE GPKG MENGGUNAKAN FILE TEMPORARY ------
        bs_4326 = bs.to_crs(4326)
        sls_4326 = sls.to_crs(4326)
        bs_flag_4326 = bs_flag.to_crs(4326)
        
        # Nama layer mengikuti nama file upload
        bs_layer_name = bs_file.name.replace(".geojson", "")
        sls_layer_name = sls_file.name.replace(".geojson", "")
        
        # Buat file sementara
        with tempfile.TemporaryDirectory() as tmpdir:
            gpkg_path = os.path.join(tmpdir, "output.gpkg")
        
            # 1. simpan bs_flag
            bs_flag_4326.to_file(gpkg_path, layer="bs_flag", driver="GPKG")
        
            # 2. append layer BS raw
            bs_4326.to_file(gpkg_path, layer=bs_layer_name, driver="GPKG", mode="a")
        
            # 3. append layer SLS raw
            sls_4326.to_file(gpkg_path, layer=sls_layer_name, driver="GPKG", mode="a")
        
            # Baca kembali ke BytesIO untuk download
            with open(gpkg_path, "rb") as f:
                gpkg_bytes = f.read()

        st.session_state["export_df"] = export_df
        st.session_state["export_df1"] = export_df1
        st.session_state["gpkg_bytes"] = gpkg_bytes
    
    st.markdown(f"#### Jumlah IDBS yang Bermasalah: {len(st.session_state["export_df"])}")
    st.dataframe(st.session_state["export_df"])

    # Pilih file output untuk unduhan
    st.download_button(
        label="Unduh Hasil CSV",
        data=st.session_state["export_df1"].to_csv(index=False),
        file_name="idbs_summary.csv",
        mime="text/csv"
    )

    import io

    # Unduh file Excel
    excel_data = io.BytesIO()
    with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
        st.session_state["export_df1"].to_excel(writer, index=False, sheet_name="Sheet1")
    excel_data.seek(0)  # Kembalikan pointer file ke awal

    # Tampilkan tombol download untuk Excel
    st.download_button(
        label="Unduh Hasil Excel",
        data=excel_data,
        file_name="idbs_summary.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
    
    # Tampilkan tombol download
    st.download_button(
        label="Unduh Hasil GeoPackage (EPSG:4326)",
        data=st.session_state["gpkg_bytes"],
        file_name="bs_vs_sls_epsg4326.gpkg",
        mime="application/geopackage+sqlite3"
    )

    
    st.success("Analisis selesai dan hasil telah disimpan!")
    st.success("Sebaiknya jangan gunakan file bs_flag sebagai file BS perbaikan, gunakan hanya sebagai alat bantu pengecekan.")

else:
    st.warning("Unggah file GeoJSON untuk BS dan SLS terlebih dahulu.")
