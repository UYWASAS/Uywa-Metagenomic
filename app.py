# ======================== BLOQUE 1: IMPORTS Y UTILIDADES ========================
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

from modules.diversity import diversity_tab
from modules.stats import stats_tab
from modules.taxonomy import taxonomy_tab
from modules.utils import load_table, safe_float, clean_state

# ======================== BLOQUE 2: ESTILO Y LOGO ========================
st.set_page_config(page_title="Microbiota 16S - UYWA", layout="wide")
st.markdown("""
    <style>
    html, body, .stApp, .main, .block-container {
        background: linear-gradient(120deg, #f3f6fa 0%, #e3ecf7 100%) !important;
        background-color: #f3f6fa !important;
        font-size: 14px !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    h1, h2, h3, h4, h5, h6 {
        font-size: 1.2em !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    .stMarkdown, .stText, .stDataFrame, .stTable, .stPlotlyChart, .stSelectbox, .stMultiSelect, .stNumberInput, .stTextInput, .stButton, .stFileUploader {
        font-size: 15px !important;
    }
    section[data-testid="stSidebar"] {
        background: #19345c !important;
        color: #fff !important;
    }
    section[data-testid="stSidebar"] * {
        color: #fff !important;
    }
    .block-container {
        background: transparent !important;
    }
    .stFileUploader, .stMultiSelect, .stSelectbox, .stNumberInput, .stTextInput {
        background-color: #f4f8fa !important;
        border-radius: 6px !important;
        border: none !important;
        box-shadow: none !important;
    }
    </style>
""", unsafe_allow_html=True)

with st.sidebar:
    st.image("assets/logo.png", width=110)
    st.markdown(
        """
        <div style='text-align: center; margin-bottom:10px;'>
            <div style='font-size:28px;font-family:Montserrat,Arial;color:#fff; margin-top: 10px;letter-spacing:1px; font-weight:700; line-height:1.1;'>
                UYWA-<br>MICROBIOTA<sup>®</sup>
            </div>
            <div style='font-size:14px;color:#fff; margin-top: 5px; font-family:Montserrat,Arial; line-height: 1.1;'>
                Análisis Interactivo 16S
            </div>
            <hr style='border-top:1px solid #2e4771; margin: 18px 0;'>
            <div style='font-size:13px;color:#fff; margin-top: 8px;'>
                <b>Contacto:</b> uywasas@gmail.com<br>
                Derechos reservados © 2025
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )
    # Opcional: Botón para cerrar sesión
    if st.button("Cerrar sesión"):
        st.session_state.clear()
        st.experimental_rerun()

# ======================== BLOQUE 3: LOGIN CON ARCHIVO AUTH.PY ROBUSTO ========================
# ======================== BLOQUE 3: LOGIN ========================
from auth import USERS_DB  # <-- IMPORTA TU ARCHIVO AUTH.PY AQUÍ

def login():
@@ -94,12 +98,8 @@
USER_KEY = f"uywa_mbio_{st.session_state['usuario']}"
st.markdown(f"<div style='text-align:right; font-size:13px;'>👤 Usuario: <b>{st.session_state['usuario']}</b></div>", unsafe_allow_html=True)

# ======================== BLOQUE 4: UTILIDADES DE SESIÓN ========================
# Ya se importó safe_float y clean_state desde utils.py

# ======================== BLOQUE 5: TITULO, UPLOADERS Y TABS PRINCIPALES ========================
st.title("Gestión y Análisis de Microbiota 16S")

tabs = st.tabs(["Carga de Archivos", "Diversidad", "Análisis Estadístico", "Visualización Taxonómica"])

# ======================== BLOQUE 6: CARGA DE ARCHIVOS EN PESTAÑA 0 ========================
@@ -112,37 +112,37 @@
    # Muestra información básica si los archivos están cargados
    if otus_file:
        df = load_table(otus_file)
        st.success(f"OTUs/ASVs: {df.shape[0]} muestras x {df.shape[1]} OTUs")
        st.success(f"OTUs/ASVs: {df.shape[0]} filas x {df.shape[1]} columnas")
    if taxonomy_file:
        df = load_table(taxonomy_file)
        st.success(f"Taxonomía: {df.shape[0]} filas x {df.shape[1]} columnas")
    if metadata_file:
        df = load_table(metadata_file)
        st.success(f"Metadata: {df.shape[0]} muestras x {df.shape[1]} variables")

    # Guarda en sesión para otras pestañas
    if otus_file: st.session_state["otus_file"] = otus_file
    if taxonomy_file: st.session_state["taxonomy_file"] = taxonomy_file
    if metadata_file: st.session_state["metadata_file"] = metadata_file

# ======================== BLOQUE 7: LLAMADA A CADA MÓDULO ========================
with tabs[1]:
    diversity_tab(
        st.session_state.get("otus_file"),
        st.session_state.get("taxonomy_file"),
        st.session_state.get("metadata_file"),
    )

with tabs[2]:
    stats_tab(
        st.session_state.get("otus_file"),
        st.session_state.get("taxonomy_file"),
        st.session_state.get("metadata_file"),
    )

with tabs[3]:
    taxonomy_tab(
        st.session_state.get("otus_file"),
        st.session_state.get("taxonomy_file"),
        st.session_state.get("metadata_file"),
    )
