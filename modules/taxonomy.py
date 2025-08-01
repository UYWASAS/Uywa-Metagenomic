import streamlit as st
import pandas as pd
import plotly.express as px
from modules.utils import load_table

def taxonomy_tab(otus_file, taxonomy_file, metadata_file):
    st.header("Visualización Taxonómica")

    # Lee archivos con index_col=0 para que OTU sea el índice
    otus = load_table(otus_file, index_col=0)
    taxonomy = load_table(taxonomy_file, index_col=0)
    metadata = load_table(metadata_file) if metadata_file else None

    if otus is None or taxonomy is None:
        st.warning("Carga archivos para visualizar taxonomía.")
        return

    # Normaliza y muestra los índices para debug
    def clean_otu_id(x):
        s = str(x).strip().upper()
        if s.endswith('.0'):  # Excel puede guardar OTU_1 como OTU_1.0
            s = s[:-2]
        return s

    otus.index = otus.index.map(clean_otu_id)
    taxonomy.index = taxonomy.index.map(clean_otu_id)

    # Debug: Muestra los IDs que no coinciden
    missing_in_tax = set(otus.index) - set(taxonomy.index)
    missing_in_otus = set(taxonomy.index) - set(otus.index)

    st.write(f"Primeros 5 OTUs en matriz: {list(otus.index[:5])}")
    st.write(f"Primeros 5 OTUs en taxonomía: {list(taxonomy.index[:5])}")
    st.write(f"OTUs en matriz y NO en taxonomía (max 5): {list(missing_in_tax)[:5]}")
    st.write(f"OTUs en taxonomía y NO en matriz (max 5): {list(missing_in_otus)[:5]}")

    # Verifica que haya al menos algún OTU compartido
    if len(set(otus.index) & set(taxonomy.index)) == 0:
        st.error("No hay coincidencias entre los OTU IDs de la matriz y la tabla de taxonomía.")
        st.stop()

    # Normaliza nombres de columna de taxonomía (por si acaso)
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
        # Intenta encontrar el nombre correcto de la columna de IDs de muestra
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
            symbol_var = None
            use_interaction = False
            if cat_vars:
                color_var = st.selectbox("Variable de agrupación", cat_vars, index=0, key=f"tax_color_{nivel}")
                use_interaction = st.checkbox("¿Mostrar interacción entre dos variables?", value=False, key=f"tax_inter_{nivel}")
                if use_interaction:
                    symbol_var = st.selectbox("Variable para interacción (símbolo)", cat_vars, index=1 if len(cat_vars) > 1 else 0, key=f"tax_symbol_{nivel}")
                    if symbol_var == color_var:
                        st.info("Selecciona dos variables diferentes para la interacción.")

            # Unir OTUs con taxonomía (por índice)
            otus_tax = otus.join(taxonomy, how="inner")
            if otus_tax.empty:
                st.warning("La unión OTU-taxonomía no produjo datos. Revisa que los IDs coincidan exactamente.")
                continue

            # Agrupar por nivel taxonómico y sumar SOLO columnas numéricas (abundancias)
            num_cols = otus.columns  # Se asume que todas las columnas en otus son numéricas
            if nivel not in otus_tax.columns:
                st.warning(f"No se encontró el nivel '{nivel}' en la tabla de taxonomía.")
                continue

            tax_sum = otus_tax.groupby(nivel)[num_cols].sum()
            if tax_sum.empty or tax_sum.shape[0] == 0:
                st.warning(f"No se encontraron datos agrupados por el nivel '{nivel}'.")
                continue

            # Selecciona top 10 taxones + Otros
            top_taxa = tax_sum.sum(axis=1).sort_values(ascending=False).head(10).index
            tax_sum_top = tax_sum.loc[top_taxa].copy()
            other_taxa = tax_sum.drop(top_taxa, errors="ignore")
            if not other_taxa.empty:
                tax_sum_top.loc["Otros"] = other_taxa.sum()

            # Normaliza a porcentaje por muestra
            tax_sum_top = tax_sum_top.T  # columnas: taxones, filas: muestras
            tax_sum_pct = tax_sum_top.div(tax_sum_top.sum(axis=1), axis=0) * 100
            tax_sum_pct.index.name = "Muestra"
            plot_df = tax_sum_pct.reset_index().melt(id_vars="Muestra", var_name=nivel, value_name="Porcentaje")

            if plot_df["Porcentaje"].isnull().all() or (plot_df["Porcentaje"].sum() == 0):
                st.warning(f"No hay datos para graficar en el nivel '{nivel}'.")
                continue

            # Añadir metadata para agrupación/interacción, si aplica
            if meta_df is not None and color_var:
                # Asegura que las columnas relevantes sean string
                plot_df["Muestra"] = plot_df["Muestra"].astype(str)
                meta_df = meta_df.copy()
                for col in meta_df.columns:
                    meta_df[col] = meta_df[col].astype(str)
                # Encuentra la columna de ID de muestra en metadata
                id_col = None
                for col in meta_df.columns:
                    if set(plot_df["Muestra"]) <= set(meta_df[col]):
                        id_col = col
                        break
                if id_col is None:
                    st.error("No se encuentra la columna de ID de muestra en la metadata para hacer merge.")
                    continue
                # Elimina duplicados en el id_col de metadata (obligatorio para merge 1:1)
                meta_df = meta_df.drop_duplicates(subset=[id_col])
                # Antes del merge, verifica unicidad en ambos lados
                if meta_df[id_col].duplicated().any():
                    st.error(f"La columna '{id_col}' en metadata tiene valores duplicados. Esto causará duplicación en el gráfico.")
                    continue
                if plot_df["Muestra"].duplicated().any():
                    st.error("La columna 'Muestra' en datos a graficar tiene duplicados. Esto causará duplicación en el gráfico.")
                    continue
                # Realiza el merge 1:1, comprobando que no cambie el número de filas
                rows_before = len(plot_df)
                plot_df = plot_df.merge(meta_df, left_on="Muestra", right_on=id_col, how="left")
                if len(plot_df) != rows_before:
                    st.error("Se generaron filas duplicadas tras el merge. Verifica la unicidad y correspondencia exacta de los IDs de muestra en metadata.")
                    continue
                # Si hay interacción entre dos variables
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
