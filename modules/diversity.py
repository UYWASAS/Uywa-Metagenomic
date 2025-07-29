import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from scipy.stats import kruskal
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
        if count > 0:
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

def plot_ellipse(fig, group_coords, color, name):
    # Calcula elipse de confianza (2D) para el grupo, solo con >=3 puntos
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    x = group_coords["PC1"].values
    y = group_coords["PC2"].values
    if len(x) < 3:
        return fig
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(vals)
    ellipse = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height, angle=theta)
    # Obtener puntos de la elipse (usamos matplotlib para generar los vértices)
    transf = plt.gca().transData
    ellipse.set_transform(transf)
    verts = ellipse.get_verts()
    fig.add_trace(go.Scatter(
        x=verts[:,0], y=verts[:,1],
        mode="lines",
        line=dict(color=color, width=2),
        name=f"{name} ellipse",
        showlegend=False,
        hoverinfo="skip"
    ))
    return fig

def diversity_tab(otus_file, taxonomy_file, metadata_file):
    st.header("Análisis de Diversidad Alfa y Beta")
    if not otus_file or not metadata_file:
        st.warning("Por favor, sube la tabla de OTUs/ASVs y la metadata en la pestaña de carga.")
        return

    otus = load_table(otus_file, index_col="OTU")
    metadata = load_table(metadata_file, index_col="SampleID")
    if otus is None or metadata is None:
        st.error("No se pudo cargar los archivos correctamente.")
        return

    # Intersección de muestras presentes en ambos archivos
    common_samples = [s for s in otus.columns if s in metadata.index]
    if not common_samples:
        st.error("No hay coincidencias entre los nombres de muestra en la tabla OTU y la metadata.")
        return
    otus = otus[common_samples]
    metadata = metadata.loc[common_samples]

    # =================== DIVERSIDAD ALFA ===================
    st.subheader("Diversidad Alfa")
    alpha_metrics = {
        "shannon": "Shannon",
        "simpson": "Simpson",
        "chao1": "Chao1",
        "observed_otus": "OTUs Observados"
    }
    alpha_df = pd.DataFrame(index=common_samples)
    otus_T = otus.T  # Filas = muestras, columnas = OTUs
    for m in alpha_metrics:
        try:
            alpha_df[alpha_metrics[m]] = alpha_diversity(m, otus_T.values, ids=otus_T.index)
        except Exception:
            alpha_df[alpha_metrics[m]] = np.nan
    alpha_df = alpha_df.join(metadata)

    # Elige variable de grupo (sólo aquellas con más de un valor y menos que el total)
    group_vars = [col for col in metadata.columns if 1 < metadata[col].nunique() < len(metadata)]
    group_col = st.selectbox("Variable de grupo para comparar", group_vars, index=0 if group_vars else None)
    if group_col:
        for m in alpha_metrics.values():
            fig = px.box(alpha_df, x=group_col, y=m, color=group_col, points="all", title=f"Índice {m}")
            st.plotly_chart(fig, use_container_width=True)
            # Estadística rápida
            groups = [alpha_df[alpha_df[group_col]==g][m].dropna() for g in alpha_df[group_col].unique()]
            if all(len(g) > 0 for g in groups) and len(groups) > 1:
                kw = kruskal(*groups)
                st.caption(f"Kruskal-Wallis p={kw.pvalue:.3g}")
            else:
                st.caption("Kruskal-Wallis p=nan")

    # =================== CURVAS DE RAREFACCIÓN ===================
    st.subheader("Curvas de Rarefacción (experimental)")
    sample_sel = st.selectbox("Muestra para rarefacción", common_samples)
    if sample_sel:
        sample = otus[sample_sel]
        sample = pd.to_numeric(sample, errors="coerce").fillna(0)
        if sample.sum() == 0:
            st.error("La muestra seleccionada está vacía. Prueba otra muestra.")
        else:
            depths, values = rarefaction_curve(sample)
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=depths, y=values, mode="lines+markers"))
            fig.update_layout(title=f"Rarefacción: {sample_sel}", xaxis_title="Profundidad", yaxis_title="OTUs Observados")
            st.plotly_chart(fig, use_container_width=True)

    # =================== DIVERSIDAD BETA (Bray-Curtis PCoA avanzado) ===================
    st.subheader("Diversidad Beta (Bray-Curtis PCoA avanzado)")
    try:
        otus_T = otus.T
        dist = beta_diversity("braycurtis", otus_T.values, ids=otus_T.index)
        coords = pcoa(dist).samples
        coords = coords.join(metadata, how="left")
        # Selección de variables para color y símbolo
        meta_cols = [col for col in metadata.columns if metadata[col].nunique() < len(metadata)]
        color_var = st.selectbox("Variable para color", meta_cols, index=0)
        symbol_var = st.selectbox("Variable para símbolo", meta_cols, index=1 if len(meta_cols) > 1 else 0)

        # Construye el scatter plot interactivo con elipses
        fig = go.Figure()
        groups = coords.groupby(color_var)
        # Paleta de colores
        palette = px.colors.qualitative.Dark24
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups.groups.keys())}
        # Mapear símbolos de Plotly
        symbol_types = ['circle', 'triangle-up', 'square', 'star', 'diamond', 'cross', 'x', 'triangle-down']
        symbol_vals = {}
        for i, val in enumerate(coords[symbol_var].unique()):
            symbol_vals[val] = symbol_types[i % len(symbol_types)]
        for i, (g, df) in enumerate(groups):
            fig.add_trace(go.Scatter(
                x=df["PC1"], y=df["PC2"],
                mode="markers",
                name=str(g),
                marker=dict(
                    color=color_map[g],
                    symbol=[symbol_vals[v] for v in df[symbol_var]],
                    size=14,
                    line=dict(width=1, color="black")
                ),
                customdata=df.index,
                hovertemplate="Sample: %{customdata}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>"+color_var+": "+str(g)
            ))
            # Dibuja la elipse para el grupo
            try:
                fig = plot_ellipse(fig, df, color_map[g], g)
            except Exception:
                pass
        fig.update_layout(
            xaxis_title="PCoA 1",
            yaxis_title="PCoA 2",
            legend_title=color_var,
            title="PCoA Bray-Curtis con elipses de grupo",
            width=800,
            height=600
        )
        st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.warning(f"No se pudo calcular PCoA Bray-Curtis: {e}")
