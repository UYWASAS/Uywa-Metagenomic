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
            <div style='font-size:32px;font-family:Montserrat,Arial;color:#fff; margin-top: 10px;letter-spacing:1px; font-weight:700; line-height:1.1;'>
                UYWA-<br>MICROBIOTA<sup>®</sup>
            </div>
            <div style='font-size:16px;color:#fff; margin-top: 5px; font-family:Montserrat,Arial; line-height: 1.1;'>
                Análisis Interactivo 16S
            </div>
            <hr style='border-top:1px solid #2e4771; margin: 18px 0;'>
            <div style='font-size:14px;color:#fff; margin-top: 8px;'>
                <b>Contacto:</b> uywasas@gmail.com<br>
                Derechos reservados © 2025
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )

# ======================== BLOQUE 3: LOGIN CON ARCHIVO AUTH.PY ROBUSTO ========================
from auth import USERS_DB  # <-- IMPORTA TU ARCHIVO AUTH.PY AQUÍ

def login():
    st.title("Iniciar sesión")
    username = st.text_input("Usuario", key="usuario_login")
    password = st.text_input("Contraseña", type="password", key="password_login")
    login_btn = st.button("Entrar", key="entrar_login")
    if login_btn:
        user = USERS_DB.get(username.strip().lower())
        if user and user["password"] == password:
            st.session_state["logged_in"] = True
            st.session_state["usuario"] = username.strip()
            st.session_state["user"] = user
            st.success(f"Bienvenido, {user['name']}!")
            st.rerun()
        else:
            st.error("Usuario o contraseña incorrectos.")
    if "logged_in" not in st.session_state or not st.session_state["logged_in"]:
        st.stop()

if "logged_in" not in st.session_state or not st.session_state["logged_in"]:
    login()

USER_KEY = f"uywa_mbio_{st.session_state['usuario']}"
st.markdown(f"<div style='text-align:right'>👤 Usuario: <b>{st.session_state['usuario']}</b></div>", unsafe_allow_html=True)

# ======================== BLOQUE 4: UTILIDADES DE SESIÓN ========================
# Ya se importó safe_float y clean_state desde utils.py

# ======================== BLOQUE 5: TITULO, UPLOADERS Y TABS PRINCIPALES ========================
st.title("Gestión y Análisis de Microbiota 16S")

with st.sidebar:
    st.header("Carga de Archivos")
    otus_file = st.file_uploader("Tabla OTUs/ASVs (csv/tsv)", type=["csv", "tsv"], key="otus_upload")
    taxonomy_file = st.file_uploader("Taxonomía (csv/tsv)", type=["csv", "tsv"], key="tax_upload")
    metadata_file = st.file_uploader("Metadata (csv/tsv/xlsx)", type=["csv", "tsv", "xlsx"], key="meta_upload")

tabs = st.tabs(["Diversidad", "Análisis Estadístico", "Visualización Taxonómica"])

# ======================== BLOQUE 6: LLAMADA A CADA MÓDULO ========================
with tabs[0]:
    diversity_tab(otus_file, taxonomy_file, metadata_file)

with tabs[1]:
    stats_tab(otus_file, taxonomy_file, metadata_file)

with tabs[2]:
    taxonomy_tab(otus_file, taxonomy_file, metadata_file)
