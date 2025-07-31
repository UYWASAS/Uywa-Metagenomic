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

    # Normaliza nombres de columna de taxonomía
    taxonomy.columns = [str(c).strip().capitalize() for c in taxonomy.columns]

    # Detecta niveles taxonómicos válidos
    tax_levels = [col for col in taxonomy.columns if taxonomy[col].nunique(dropna=True) > 1]
    if not tax_levels:
        st.warning("No se detectaron niveles taxonómicos múltiples en el archivo de taxonomía.")
        return

    # Variables categóricas de metadata
    cat_vars = []
    meta_df = None
    if metadata is not None:
        meta_df = metadata.reset_index()
        if "Sampleid" not in [c.lower() for c in meta_df.columns]:
            if "index" in meta_df.columns:
                meta_df = meta_df.rename(columns={"index": "SampleID"})
        cat_vars = [col for col in meta_df.columns if 1 < meta_df[col].nunique() < len(meta_df)]

    tabs = st.tabs(tax_levels)
    for i, nivel in enumerate(tax_levels):
        with tabs[i]:
            st.subheader(f"Barplot apilado por {nivel} (top 10 + Otros)")

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

            # DEPURACIÓN: muestra las dimensiones y los índices/columnas originales
            st.write("Dimensiones OTUs:", otus.shape)
            st.write("Primeras columnas OTUs:", otus.columns[:5])
            st.write("Dimensiones Taxonomía:", taxonomy.shape)
            st.write("Primeros índices Taxonomía:", taxonomy.index[:5])
            st.write("Primeras filas Taxonomía:", taxonomy.head())

            # Forzar tipo str en índices y columnas
            taxonomy.index = taxonomy.index.astype(str)
            otus.columns = otus.columns.astype(str)

            # Solo usar OTUs presentes en ambos
            comunes = [otu for otu in otus.columns if otu in taxonomy.index]
            if len(comunes) == 0:
                st.error("No hay coincidencias entre OTU IDs en la matriz y la tabla de taxonomía.")
                return

            otus_T = otus[comunes].T.copy()
            tax_for_otus = taxonomy.loc[otus_T.index]

            # DEPURACIÓN: muestra después del join
            st.write("Dimensiones OTUs_T (tras filtrar):", otus_T.shape)
            st.write("Dimensiones Tax_for_otus (tras join):", tax_for_otus.shape)

            otus_tax = otus_T.copy()
            otus_tax[nivel] = tax_for_otus[nivel].values

            # DEPURACIÓN: muestra después de añadir columna nivel
            st.write("Dimensiones otus_tax:", otus_tax.shape)
            st.dataframe(otus_tax.head())

            # Agrupa por ese nivel taxonómico y suma
            tax_sum = otus_tax.groupby(nivel).sum()
            st.write("tax_sum shape:", tax_sum.shape)
            st.dataframe(tax_sum.head())

            # Solo deja top 10 + Otros
            top_taxa = tax_sum.sum().sort_values(ascending=False).head(10).index
            tax_sum_top = tax_sum.loc[top_taxa]
            other_cols = [col for col in tax_sum.index if col not in top_taxa]
            if other_cols:
                sum_otros = tax_sum.loc[other_cols].sum()
                tax_sum_top.loc["Otros"] = sum_otros
            tax_sum_top = tax_sum_top.T
            # Normaliza a porcentaje por muestra
            tax_sum_pct = tax_sum_top.div(tax_sum_top.sum(axis=1), axis=0) * 100
            tax_sum_pct.index.name = "Muestra"
            st.write("tax_sum_pct shape:", tax_sum_pct.shape)
            st.dataframe(tax_sum_pct.head())

            plot_df = tax_sum_pct.reset_index().melt(id_vars="Muestra", var_name=nivel, value_name="Porcentaje")

            # DEPURACIÓN: muestra la tabla final antes del plot
            st.dataframe(plot_df.head(20))
            st.write(plot_df["Porcentaje"].describe())

            if plot_df["Porcentaje"].isnull().all() or (plot_df["Porcentaje"].sum() == 0):
                st.warning(f"No hay datos para graficar en el nivel '{nivel}'.")
                continue

            if meta_df is not None and color_var:
                id_col = None
                for col in meta_df.columns:
                    if set(plot_df["Muestra"]).issubset(set(meta_df[col].astype(str))):
                        id_col = col
                        break
                if id_col is None:
                    st.error("No se encuentra la columna de ID de muestra en la metadata para hacer merge.")
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
