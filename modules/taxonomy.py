import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from modules.utils import load_table

def taxonomy_tab(otus_file, taxonomy_file, metadata_file):
    st.header("Visualización Taxonómica")

    otus = load_table(otus_file, index_col=0)
    taxonomy = load_table(taxonomy_file, index_col=0)
    metadata = load_table(metadata_file) if metadata_file else None

    if otus is None or taxonomy is None:
        st.warning("Carga archivos para visualizar taxonomía.")
        return

    def clean_otu_id(x):
        s = str(x).strip().upper()
        if s.endswith('.0'):
            s = s[:-2]
        return s

    otus.index = otus.index.map(clean_otu_id)
    taxonomy.index = taxonomy.index.map(clean_otu_id)

    missing_in_tax = set(otus.index) - set(taxonomy.index)
    missing_in_otus = set(taxonomy.index) - set(otus.index)

    st.write(f"Primeros 5 OTUs en matriz: {list(otus.index[:5])}")
    st.write(f"Primeros 5 OTUs en taxonomía: {list(taxonomy.index[:5])}")
    st.write(f"OTUs en matriz y NO en taxonomía (max 5): {list(missing_in_tax)[:5]}")
    st.write(f"OTUs en taxonomía y NO en matriz (max 5): {list(missing_in_otus)[:5]}")

    if len(set(otus.index) & set(taxonomy.index)) == 0:
        st.error("No hay coincidencias entre los OTU IDs de la matriz y la tabla de taxonomía.")
        st.stop()

    taxonomy.columns = [str(c).strip().capitalize() for c in taxonomy.columns]
    tax_levels = [col for col in taxonomy.columns if taxonomy[col].nunique(dropna=True) > 1]
    if not tax_levels:
        st.warning("No se detectaron niveles taxonómicos múltiples en el archivo de taxonomía.")
        return

    cat_vars = []
    meta_df = None
    if metadata is not None:
        meta_df = metadata.reset_index()
        meta_df.columns = [str(c).strip() for c in meta_df.columns]
        if "Sampleid" not in [c.lower() for c in meta_df.columns]:
            if "index" in meta_df.columns:
                meta_df = meta_df.rename(columns={"index": "SampleID"})
        cat_vars = [col for col in meta_df.columns if 1 < meta_df[col].nunique() < len(meta_df)]

    tabs = st.tabs(tax_levels)
    for i, nivel in enumerate(tax_levels):
        with tabs[i]:
            st.subheader(f"Barplot apilado por {nivel} (top 10 + Otros)")
            color_var = None
            if cat_vars:
                color_var = st.selectbox("Variable de agrupación", cat_vars, index=0, key=f"tax_color_{nivel}")

            otus_tax = otus.join(taxonomy, how="inner")
            if otus_tax.empty:
                st.warning("La unión OTU-taxonomía no produjo datos. Revisa que los IDs coincidan exactamente.")
                continue

            num_cols = otus.columns
            if nivel not in otus_tax.columns:
                st.warning(f"No se encontró el nivel '{nivel}' en la tabla de taxonomía.")
                continue

            tax_sum = otus_tax.groupby(nivel)[num_cols].sum()
            if tax_sum.empty or tax_sum.shape[0] == 0:
                st.warning(f"No se encontraron datos agrupados por el nivel '{nivel}'.")
                continue

            top_taxa = tax_sum.sum(axis=1).sort_values(ascending=False).head(10).index
            tax_sum_top = tax_sum.loc[top_taxa].copy()
            other_taxa = tax_sum.drop(top_taxa, errors="ignore")
            if not other_taxa.empty:
                tax_sum_top.loc["Otros"] = other_taxa.sum()

            tax_sum_top = tax_sum_top.T
            tax_sum_pct = tax_sum_top.div(tax_sum_top.sum(axis=1), axis=0) * 100
            tax_sum_pct.index.name = "Muestra"
            plot_df = tax_sum_pct.reset_index().melt(id_vars="Muestra", var_name=nivel, value_name="Porcentaje")

            if plot_df["Porcentaje"].isnull().all() or (plot_df["Porcentaje"].sum() == 0):
                st.warning(f"No hay datos para graficar en el nivel '{nivel}'.")
                continue

            if meta_df is not None and color_var:
                plot_df["Muestra"] = plot_df["Muestra"].astype(str)
                meta_df = meta_df.copy()
                for col in meta_df.columns:
                    meta_df[col] = meta_df[col].astype(str)
                # Encuentra la columna de ID de muestra
                id_col = None
                for col in meta_df.columns:
                    if set(plot_df["Muestra"]).issubset(set(meta_df[col])):
                        id_col = col
                        break
                if id_col is None:
                    st.error("No se encuentra la columna de ID de muestra en la metadata para hacer merge.")
                    return
                meta_df = meta_df.drop_duplicates(subset=[id_col])

                # --- SELECCIÓN MÚLTIPLE DE GRUPOS ---
                unique_groups = meta_df[color_var].dropna().unique()
                selected_groups = st.multiselect(
                    f"Selecciona uno o varios valores de '{color_var}' para comparar:",
                    unique_groups,
                    default=list(unique_groups),  # por defecto todos
                    key=f"multiselect_group_{nivel}"
                )

                muestras_en_grupo = meta_df.loc[meta_df[color_var].isin(selected_groups), id_col]
                plot_df_group = plot_df[plot_df["Muestra"].isin(muestras_en_grupo)]
                # Haz el merge para obtener la columna de grupo en plot_df_group
                plot_df_group = plot_df_group.merge(meta_df[[id_col, color_var]], left_on="Muestra", right_on=id_col, how="left")
                plot_df_group = plot_df_group[plot_df_group["Porcentaje"] > 0]

                # Ordenar por grupo y muestra (simple y robusto)
                plot_df_group = plot_df_group.sort_values([color_var, "Muestra"])

                # Gráfico manteniendo solo color por taxón
                fig = px.bar(
                    plot_df_group,
                    x="Muestra", y="Porcentaje", color=nivel,
                    title=f"Abundancia relativa por {nivel} en {color_var}: {', '.join([str(g) for g in selected_groups])}",
                    labels={"Porcentaje": "% abundancia relativa", color_var: color_var}
                )
                fig.update_layout(barmode="stack", xaxis_title="Muestra", yaxis_title="% abundancia relativa")

                # --- ANOTACIONES Y LÍNEAS DIVISORIAS PARA GRUPOS ---
                muestras_x = plot_df_group[["Muestra", color_var]].drop_duplicates().reset_index(drop=True)
                for group in selected_groups:
                    muestras_grupo = muestras_x[muestras_x[color_var] == group]["Muestra"]
                    if muestras_grupo.empty:
                        continue
                    x0 = muestras_grupo.index[0]
                    x1 = muestras_grupo.index[-1]
                    # Línea divisoria antes del primer elemento del grupo (excepto el primero)
                    if x0 != 0:
                        fig.add_vline(x=x0-0.5, line_width=1, line_dash="dot", line_color="grey")
                    # Anotación centrada sobre las muestras del grupo (encima del gráfico)
                    xpos = (x0 + x1) / 2
                    fig.add_annotation(
                        x=xpos, y=105,
                        text=f"{color_var}: {group}",
                        showarrow=False,
                        font=dict(size=13, color="black"),
                        xref="x", yref="y"
                    )

                st.plotly_chart(fig, use_container_width=True)
            else:
                plot_df = plot_df[plot_df["Porcentaje"] > 0]
                fig = px.bar(
                    plot_df,
                    x="Muestra", y="Porcentaje", color=nivel,
                    title=f"Abundancia relativa por {nivel} (Top 10 + Otros)",
                    labels={"Porcentaje": "% abundancia relativa"}
                )
                fig.update_layout(barmode="stack", xaxis_title="Muestra", yaxis_title="% abundancia relativa")
                st.plotly_chart(fig, use_container_width=True)
