import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, DistanceMatrix
from skbio.stats.ordination import nmds
from scipy.stats import kruskal, f_oneway
from modules.utils import load_table
import numpy as np

def rarefaction_curve(sample, steps=10):
    nseqs = int(sample.sum())
    depths = np.linspace(10, nseqs, steps, dtype=int)
    values = []
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

def get_ellipse(x, y, n_std=2.0, num_points=100):
    """Devuelve puntos para una elipse de dispersión tipo confidence ellipse."""
    if len(x) < 3:
        return None, None
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    t = np.linspace(0, 2 * np.pi, num_points)
    ellipse = np.array([width/2 * np.cos(t) , height/2 * np.sin(t)])
    R = np.array([[np.cos(np.radians(theta)), -np.sin(np.radians(theta))],
                  [np.sin(np.radians(theta)),  np.cos(np.radians(theta))]])
    ellipse_rot = R @ ellipse
    x0, y0 = np.mean(x), np.mean(y)
    return ellipse_rot[0] + x0, ellipse_rot[1] + y0

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

    # Selección interactiva de variables para color y símbolo
    cat_vars = [col for col in metadata.columns if 1 < metadata[col].nunique() < len(metadata)]
    color_var = st.selectbox("Variable para color en boxplot", cat_vars, index=0 if cat_vars else None, key="alpha_color")
    symbol_var = st.selectbox("Variable para símbolo en boxplot", cat_vars, index=1 if len(cat_vars)>1 else 0, key="alpha_symbol")

    # Boxplot + puntos + interacción de variables + elipses en 2D (scatter de índices)
    for metric in alpha_metrics.values():
        st.markdown(f"**{metric}**")
        fig = px.box(
            alpha_df, x=color_var, y=metric, color=color_var, points="all",
            symbol=symbol_var if symbol_var != color_var else None,
            title=f"{metric} por grupo"
        )
        # Overlay: scatter (para mostrar interacción/color/símbolo)
        fig_scatter = px.scatter(
            alpha_df, x=symbol_var, y=metric, color=color_var, symbol=symbol_var if symbol_var != color_var else None,
            title=f"{metric}: interacción {color_var} y {symbol_var}"
        )
        st.plotly_chart(fig, use_container_width=True)
        st.plotly_chart(fig_scatter, use_container_width=True)

        # Estadística: ANOVA o Kruskal-Wallis
        groups = [alpha_df[alpha_df[color_var]==g][metric].dropna() for g in alpha_df[color_var].unique()]
        if all(len(g) > 1 for g in groups) and len(groups) > 1:
            try:
                fval, pval = f_oneway(*groups)
                st.caption(f"ANOVA: p = {pval:.3g}")
            except Exception:
                kw = kruskal(*groups)
                st.caption(f"Kruskal-Wallis: p = {kw.pvalue:.3g}")
        else:
            st.caption("No hay replicación suficiente para ANOVA/Kruskal-Wallis.")

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

    # =================== DIVERSIDAD BETA (NMDS + elipses) ===================
    st.subheader("Diversidad Beta (NMDS Bray-Curtis + Elipses)")
    try:
        otus_T = otus.T
        dist = beta_diversity("braycurtis", otus_T.values, ids=otus_T.index)
        # NMDS (2D)
        nmds_res = nmds(distance_matrix=dist, n_components=2, random_state=42)
        coords = pd.DataFrame(nmds_res.samples.values, index=otus_T.index, columns=["NMDS1", "NMDS2"])
        coords = coords.join(metadata, how="left")
        # Selección de variables para color y símbolo
        color_var_beta = st.selectbox("Variable para color NMDS", cat_vars, index=0, key="beta_color")
        symbol_var_beta = st.selectbox("Variable para símbolo NMDS", cat_vars, index=1 if len(cat_vars)>1 else 0, key="beta_symbol")

        fig = go.Figure()
        palette = px.colors.qualitative.Dark24
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(coords[color_var_beta].unique())}
        # Símbolos para grupos
        symbol_types = ['circle', 'triangle-up', 'square', 'star', 'diamond', 'cross', 'x', 'triangle-down']
        symbol_vals = {}
        for i, val in enumerate(coords[symbol_var_beta].unique()):
            symbol_vals[val] = symbol_types[i % len(symbol_types)]
        # Dibuja puntos y elipses por grupo
        for i, (group, group_df) in enumerate(coords.groupby(color_var_beta)):
            fig.add_trace(go.Scatter(
                x=group_df["NMDS1"], y=group_df["NMDS2"],
                mode="markers",
                name=str(group),
                marker=dict(
                    color=color_map[group],
                    symbol=[symbol_vals[v] for v in group_df[symbol_var_beta]],
                    size=14,
                    line=dict(width=1, color="black")
                ),
                customdata=group_df.index,
                hovertemplate="Sample: %{customdata}<br>NMDS1: %{x:.2f}<br>NMDS2: %{y:.2f}<br>"+color_var_beta+": "+str(group)
            ))
            # Elipse para el grupo
            ex, ey = get_ellipse(group_df["NMDS1"].values, group_df["NMDS2"].values)
            if ex is not None:
                fig.add_trace(go.Scatter(
                    x=ex, y=ey,
                    mode="lines",
                    line=dict(color=color_map[group], width=2),
                    name=f"{group} ellipse",
                    showlegend=False,
                    hoverinfo="skip"
                ))
        fig.update_layout(
            xaxis_title="NMDS1",
            yaxis_title="NMDS2",
            legend_title=color_var_beta,
            title="NMDS Bray-Curtis con elipses de grupo",
            width=800,
            height=600
        )
        st.plotly_chart(fig, use_container_width=True)

        # PERMANOVA
        if color_var_beta:
            try:
                dm = DistanceMatrix(dist.data, ids=dist.ids)
                permanova_res = permanova(dm, grouping=coords[color_var_beta], permutations=999)
                st.caption(f"PERMANOVA (adonis) para {color_var_beta}: p = {permanova_res['p-value']:.3g}, pseudo-F = {permanova_res['test statistic']:.3g}, R2 = {permanova_res['test statistic']/permanova_res['denominator']:.3g}")
            except Exception as e:
                st.caption(f"No se pudo calcular PERMANOVA: {e}")

    except Exception as e:
        st.warning(f"No se pudo calcular NMDS Bray-Curtis: {e}")
