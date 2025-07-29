import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from scipy.stats import kruskal, f_oneway
from modules.utils import load_table
import numpy as np

def rarefaction_curve(sample, steps=10):
    """
    Crea una curva de rarefacción para una muestra dada (serie de abundancias).
    """
    nseqs = int(sample.sum())
    depths = np.linspace(10, nseqs, steps, dtype=int)
    values = []
    # Expandir la muestra: cada OTU tantas veces como su abundancia
    population = []
    for otu, count in sample.items():
        population.extend([otu] * int(count))
    population = np.array(population)
    for d in depths:
        if d > len(population) or d == 0:
            values.append(np.nan)
            continue
        outs = []
        for _ in range(10):
            subsample = np.random.choice(population, d, replace=False)
            outs.append(len(np.unique(subsample)))
        values.append(np.mean(outs))
    return depths, values

def diversity_tab(otus_file, taxonomy_file, metadata_file):
    st.header("Análisis de Diversidad Alfa y Beta")
    if not otus_file or not metadata_file:
        st.warning("Por favor, sube la tabla de OTUs/ASVs y la metadata en la pestaña de carga.")
        return

    otus = load_table(otus_file)
    metadata = load_table(metadata_file)
    if otus is None or metadata is None:
        st.error("No se pudo cargar los archivos correctamente.")
        return

    # =================== DIVERSIDAD ALFA ===================
    st.subheader("Diversidad Alfa")
    alpha_metrics = {
        "shannon": "Shannon",
        "simpson": "Simpson",
        "chao1": "Chao1",
        "observed_otus": "OTUs Observados"
    }
    alpha_res = {}
    for m in alpha_metrics:
        try:
            res = alpha_diversity(m, otus.values, ids=otus.index)
            alpha_res[alpha_metrics[m]] = res
        except Exception:
            alpha_res[alpha_metrics[m]] = [None] * len(otus.index)
    alpha_df = pd.DataFrame(alpha_res, index=otus.index)
    alpha_df = alpha_df.join(metadata)

    group_col = st.selectbox("Variable de grupo para comparar", metadata.columns, index=0)
    for idx, m in enumerate(alpha_metrics.values()):
        fig = px.box(alpha_df, x=group_col, y=m, color=group_col, points="all", title=f"Índice {m}")
        st.plotly_chart(fig, use_container_width=True)
        # Estadística rápida
        groups = [alpha_df[alpha_df[group_col]==g][m].dropna() for g in alpha_df[group_col].unique()]
        if len(groups) > 1:
            kw = kruskal(*groups)
            st.caption(f"Kruskal-Wallis p={kw.pvalue:.3g}")

    # =================== CURVAS DE RAREFACCIÓN ===================
    st.subheader("Curvas de Rarefacción (experimental)")
    sample_sel = st.selectbox("Muestra para rarefacción", otus.index)
    if sample_sel:
        sample = otus.loc[sample_sel]
        depths, values = rarefaction_curve(sample)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=depths, y=values, mode="lines+markers"))
        fig.update_layout(title=f"Rarefacción: {sample_sel}", xaxis_title="Profundidad", yaxis_title="OTUs Observados")
        st.plotly_chart(fig, use_container_width=True)

    # =================== DIVERSIDAD BETA ===================
    st.subheader("Diversidad Beta (Bray-Curtis PCoA)")
    try:
        dist = beta_diversity("braycurtis", otus.values, ids=otus.index)
        coords = pcoa(dist).samples
        coords = coords.join(metadata)
        fig = px.scatter(coords, x="PC1", y="PC2", color=group_col, symbol=group_col, title="PCoA Bray-Curtis")
        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.warning(f"No se pudo calcular PCoA Bray-Curtis: {e}")
