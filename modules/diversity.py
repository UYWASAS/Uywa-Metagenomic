import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova, DistanceMatrix
from scipy.stats import kruskal, f_oneway
from sklearn.manifold import MDS
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
    ellipse = np.array([width/2 * np.cos(t), height/2 * np.sin(t)])
    R = np.array([[np.cos(np.radians(theta)), -np.sin(np.radians(theta))],
                  [np.sin(np.radians(theta)),  np.cos(np.radians(theta))]])
    ellipse_rot = R @ ellipse
    x0, y0 = np.mean(x), np.mean(y)
    return ellipse_rot[0] + x0, ellipse_rot[1] + y0

def plot_alpha_index_tabbed(alpha_df, cat_vars, alpha_metrics):
    color_var = st.selectbox("Variable de agrupación", cat_vars, index=0 if cat_vars else None, key="alpha_color")
    use_interaction = st.checkbox("¿Mostrar interacción entre dos variables? (alfa diversidad)", value=False)
    symbol_var = None
    if use_interaction:
        symbol_var = st.selectbox("Variable para interacción (símbolo)", cat_vars, index=1 if len(cat_vars) > 1 else 0, key="alpha_symbol")
        if symbol_var == color_var:
            st.info("Selecciona dos variables diferentes para la interacción.")

    tabs = st.tabs(list(alpha_metrics.values()))
    for i, metric in enumerate(alpha_metrics.values()):
        with tabs[i]:
            st.markdown(f"**{metric}**")
            fig = px.box(
                alpha_df, x=color_var, y=metric, color=color_var, points="all",
                title=f"{metric} por grupo"
            )
            st.plotly_chart(fig, use_container_width=True)
            if use_interaction and symbol_var and symbol_var != color_var:
                fig_scatter = px.scatter(
                    alpha_df, x=symbol_var, y=metric, color=color_var, symbol=symbol_var,
                    title=f"{metric}: interacción {color_var} y {symbol_var}"
                )
                st.plotly_chart(fig_scatter, use_container_width=True)
            groups = [alpha_df[alpha_df[color_var] == g][metric].dropna() for g in alpha_df[color_var].unique()]
            if all(len(g) > 1 for g in groups) and len(groups) > 1:
                try:
                    fval, pval = f_oneway(*groups)
                    st.caption(f"ANOVA: p = {pval:.3g}")
                except Exception:
                    kw = kruskal(*groups)
                    st.caption(f"Kruskal-Wallis: p = {kw.pvalue:.3g}")
            else:
                st.caption("No hay replicación suficiente para ANOVA/Kruskal-Wallis.")

def plot_beta_diversity(coords, metadata, dist, cat_vars_beta):
    color_var_beta = st.selectbox("Variable para color NMDS", cat_vars_beta, index=0, key="beta_color")
    use_interaction_beta = st.checkbox("¿Mostrar interacción entre dos variables? (beta diversidad)", value=False)
    symbol_var_beta = None
    if use_interaction_beta:
        symbol_var_beta = st.selectbox("Variable para símbolo NMDS", cat_vars_beta, index=1 if len(cat_vars_beta) > 1 else 0, key="beta_symbol")
        if symbol_var_beta == color_var_beta:
            st.info("Selecciona dos variables diferentes para la interacción.")
    else:
        symbol_var_beta = color_var_beta

    fig = go.Figure()
    palette = px.colors.qualitative.Dark24
    color_map = {g: palette[i % len(palette)] for i, g in enumerate(coords[color_var_beta].unique())}
    symbol_types = ['circle', 'triangle-up', 'square', 'star', 'diamond', 'cross', 'x', 'triangle-down']

    if use_interaction_beta and symbol_var_beta and symbol_var_beta != color_var_beta:
        symbol_vals = {}
        for i, val in enumerate(coords[symbol_var_beta].unique()):
            symbol_vals[val] = symbol_types[i % len(symbol_types)]
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
                hovertemplate="Sample: %{customdata}<br>NMDS1: %{x:.2f}<br>NMDS2: %{y:.2f}<br>" + color_var_beta + ": " + str(group)
            ))
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
    else:
        for i, (group, group_df) in enumerate(coords.groupby(color_var_beta)):
            fig.add_trace(go.Scatter(
                x=group_df["NMDS1"], y=group_df["NMDS2"],
                mode="markers",
                name=str(group),
                marker=dict(
                    color=color_map[group],
                    size=14,
                    line=dict(width=1, color="black")
                ),
                customdata=group_df.index,
                hovertemplate="Sample: %{customdata}<br>NMDS1: %{x:.2f}<br>NMDS2: %{y:.2f}<br>" + color_var_beta + ": " + str(group)
            ))
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

    # PERMANOVA debug y limpieza
    if color_var_beta and color_var_beta in metadata.columns:
        try:
            grouping = metadata.loc[coords.index, color_var_beta]
            group_counts = grouping.value_counts()
            st.write("PERMANOVA grouping (antes de limpiar):", grouping)
            st.write("Grouping unique values (antes):", grouping.unique())
            grouping = grouping.astype(str).str.replace(r"[^a-zA-Z0-9_\-]", "_", regex=True)
            st.write("PERMANOVA grouping (después de limpiar):", grouping)
            st.write("Grouping unique values (después):", grouping.unique())
            st.write("Grouping as list:", list(grouping))

            # NUEVO: Debug de índices y matriz de distancias
            st.write("IDs en dist:", dist.ids)
            st.write("Índice en grouping:", list(grouping.index))
            st.write("Shape de dist.data:", dist.data.shape)
            st.write("dist.data:", dist.data)

            # Prueba manual con subconjunto reducido
            try:
                test_ids = dist.ids[:5]
                test_dm = DistanceMatrix(dist.data[:5, :5], ids=test_ids)
                test_grouping = grouping.loc[test_ids]
                st.write("Test grouping:", test_grouping)
                test_result = permanova(test_dm, grouping=test_grouping, permutations=99)
                st.write("Test PERMANOVA:", test_result)
            except Exception as e:
                st.write("Error en test PERMANOVA:", e)

            # PERMANOVA real (solo si hay >1 grupo y cada grupo >1 muestra)
            if grouping.nunique() > 1 and all(group_counts > 1):
                dm = DistanceMatrix(dist.data, ids=dist.ids)
                permanova_res = permanova(dm, grouping=grouping, permutations=999)
                st.write("PERMANOVA result (type):", type(permanova_res))
                st.write("PERMANOVA result (value):", permanova_res)
                pval = stat = r2 = None
                if hasattr(permanova_res, "iloc") and hasattr(permanova_res, "columns"):
                    row = permanova_res.iloc[0]
                    pval = row.get('p-value', row.get('p_value', None)) if hasattr(row, 'get') else row['p-value'] if 'p-value' in row else None
                    stat = row.get('test statistic', row.get('statistic', None)) if hasattr(row, 'get') else row['test statistic'] if 'test statistic' in row else (row['statistic'] if 'statistic' in row else None)
                    r2 = row.get('r2', None) if hasattr(row, 'get') else row['r2'] if 'r2' in row else None
                elif isinstance(permanova_res, dict):
                    pval = permanova_res.get('p-value', permanova_res.get('p_value'))
                    stat = permanova_res.get('test statistic', permanova_res.get('statistic'))
                    r2 = permanova_res.get('r2')
                else:
                    try:
                        pval = getattr(permanova_res, 'p_value', None)
                        stat = getattr(permanova_res, 'test_statistic', getattr(permanova_res, 'statistic', None))
                        r2 = getattr(permanova_res, 'r2', None)
                    except Exception:
                        pass
                msg = f"PERMANOVA (adonis) para {color_var_beta}: "
                msg += f"p = {pval}, " if pval is not None else "p = None, "
                msg += f"pseudo-F = {stat}" if stat is not None else "pseudo-F = None"
                msg += f", R2 = {r2}" if r2 is not None else ""
                st.caption(msg)
            else:
                st.caption("PERMANOVA requiere al menos 2 grupos con más de 1 muestra cada uno.")
        except Exception as e:
            st.caption(f"No se pudo calcular PERMANOVA: {e}")

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
    otus_T = otus.T
    for m in alpha_metrics:
        try:
            alpha_df[alpha_metrics[m]] = alpha_diversity(m, otus_T.values, ids=otus_T.index)
        except Exception:
            alpha_df[alpha_metrics[m]] = np.nan
    alpha_df = alpha_df.join(metadata)
    cat_vars = [col for col in metadata.columns if 1 < metadata[col].nunique() < len(metadata)]
    plot_alpha_index_tabbed(alpha_df, cat_vars, alpha_metrics)

    # =================== DIVERSIDAD BETA (NMDS + elipses) ===================
    st.subheader("Diversidad Beta (NMDS Bray-Curtis + Elipses)")
    try:
        otus_T = otus.T
        dist = beta_diversity("braycurtis", otus_T.values, ids=otus_T.index)
        mds = MDS(n_components=2, metric=False, dissimilarity='precomputed', random_state=42, n_init=10, max_iter=300)
        nmds_coords = mds.fit_transform(dist.data)
        coords = pd.DataFrame(nmds_coords, index=otus_T.index, columns=["NMDS1", "NMDS2"])
        coords = coords.join(metadata, how="left")
        cat_vars_beta = [col for col in metadata.columns if 1 < metadata[col].nunique() < len(metadata)]
        plot_beta_diversity(coords, metadata, dist, cat_vars_beta)
    except Exception as e:
        st.warning(f"No se pudo calcular NMDS Bray-Curtis: {e}")

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
