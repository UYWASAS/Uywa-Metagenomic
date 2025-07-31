import streamlit as st
import pandas as pd
import plotly.express as px
from modules.utils import load_table

def taxonomy_tab(otus_file, taxonomy_file, metadata_file):
    st.header("Visualización Taxonómica")
    otus = load_table(otus_file)
    taxonomy = load_table(taxonomy_file)
    metadata = load_table(metadata_file)
    if otus is None or taxonomy is None:
        st.warning("Carga archivos para visualizar taxonomía.")
        return

    # Detecta niveles taxonómicos válidos: más de un valor único
    tax_levels = [col for col in taxonomy.columns if taxonomy[col].nunique() > 1]
    if not tax_levels:
        st.warning("No se detectaron niveles taxonómicos múltiples en el archivo de taxonomía.")
        return

    # Detecta variables categóricas de la metadata
    cat_vars = []
    if metadata is not None:
        meta_df = metadata.reset_index()
        if "SampleID" not in meta_df.columns:
            if "index" in meta_df.columns:
                meta_df = meta_df.rename(columns={"index": "SampleID"})
        cat_vars = [col for col in meta_df.columns if 1 < meta_df[col].nunique() < len(meta_df)]

    tabs = st.tabs(tax_levels)
    for i, nivel in enumerate(tax_levels):
        with tabs[i]:
            st.subheader(f"Barplot apilado por {nivel} (top 10 + Otros)")

            # Selección de variable de agrupación y posible interacción como en diversidad
            color_var = None
            symbol_var = None
            use_interaction = False
            if cat_vars:
                color_var = st.selectbox("Variable de agrupación", cat_vars, index=0, key=f"tax_color_{nivel}")
                use_interaction = st.checkbox("¿Mostrar interacción entre dos variables?", value=False, key=f"tax_inter_{nivel}")
                if use_interaction:
                    symbol_var = st.selectbox("Variable para interacción (símbolo)", cat_vars, index=1 if len(cat_vars) > 1 else 0, key=f"tax_symbol_{nivel}")
                    if symbol_var == color_var:
                        st.info("Selecciona dos variables diferentes para la interacción.")

            # Procesamiento de los datos taxonómicos
            otus_tax = otus.T.join(taxonomy[nivel])
            tax_sum = otus_tax.groupby(nivel).sum().T
            top_taxa = tax_sum.sum().sort_values(ascending=False).head(10).index
            tax_sum_top = tax_sum[top_taxa]
            other_cols = [col for col in tax_sum.columns if col not in top_taxa]
            if other_cols:
                tax_sum_top["Otros"] = tax_sum[other_cols].sum(axis=1)
            tax_sum_pct = tax_sum_top.div(tax_sum_top.sum(axis=1), axis=0) * 100
            tax_sum_pct.index.name = "Muestra"
            plot_df = tax_sum_pct.reset_index().melt(id_vars="Muestra", var_name=nivel, value_name="Porcentaje")

            # Añadir metadata para agrupación/interacción (merge robusto)
            if metadata is not None and color_var:
                meta_df = metadata.reset_index()
                # Encuentra la columna de ID de muestra para hacer merge
                id_col = None
                for col in meta_df.columns:
                    if set(plot_df["Muestra"]).issubset(set(meta_df[col])):
                        id_col = col
                        break
                if id_col is None:
                    st.error("No se encuentra la columna de ID de muestra en la metadata.")
                    continue
                plot_df = plot_df.merge(meta_df, left_on="Muestra", right_on=id_col, how="left")
                if use_interaction and symbol_var and symbol_var != color_var:
                    fig = px.bar(
                        plot_df,
                        x="Muestra", y="Porcentaje", color=nivel,
                        facet_col=color_var,
                        pattern_shape=symbol_var,
                        title=f"Abundancia relativa por {nivel} agrupado por {color_var} y {symbol_var}",
                        labels={"Porcentaje": "% abundancia relativa"}
                    )
                else:
                    fig = px.bar(
                        plot_df,
                        x="Muestra", y="Porcentaje", color=nivel,
                        facet_col=color_var,
                        title=f"Abundancia relativa por {nivel} agrupado por {color_var}",
                        labels={"Porcentaje": "% abundancia relativa"}
                    )
            else:
                fig = px.bar(
                    plot_df,
                    x="Muestra", y="Porcentaje", color=nivel,
                    title=f"Abundancia relativa por {nivel} (Top 10 + Otros)",
                    labels={"Porcentaje": "% abundancia relativa"}
                )

            fig.update_layout(barmode="stack", xaxis_title="Muestra", yaxis_title="% abundancia relativa")
            st.plotly_chart(fig, use_container_width=True)
