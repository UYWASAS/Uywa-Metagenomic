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

    # Detecta niveles taxonómicos (Phylum, Class, Order, Family, Genus, Species...)
    tax_levels = [col for col in taxonomy.columns if taxonomy[col].nunique() > 1]
    if not tax_levels:
        st.warning("No se detectaron niveles taxonómicos múltiples en el archivo de taxonomía.")
        return

    nivel = st.selectbox("Nivel taxonómico", tax_levels, index=tax_levels.index("Phylum") if "Phylum" in tax_levels else 0)

    st.subheader(f"Barplot apilado por {nivel} (top 10)")
    if nivel in taxonomy.columns:
        # Suma por nivel taxonómico
        otus_tax = otus.T.join(taxonomy[nivel])
        tax_sum = otus_tax.groupby(nivel).sum().T
        top_taxa = tax_sum.sum().sort_values(ascending=False).head(10).index
        tax_sum_top = tax_sum[top_taxa]
        # Agrupa en "Otros" el resto
        other_cols = [col for col in tax_sum.columns if col not in top_taxa]
        if other_cols:
            tax_sum_top["Otros"] = tax_sum[other_cols].sum(axis=1)
        # Normaliza a porcentaje por muestra
        tax_sum_pct = tax_sum_top.div(tax_sum_top.sum(axis=1), axis=0) * 100
        # Prepara formato largo para plotly
        plot_df = tax_sum_pct.reset_index().melt(id_vars="index", var_name=nivel, value_name="Porcentaje")
        plot_df = plot_df.rename(columns={"index": "Muestra"})

        fig = px.bar(
            plot_df,
            x="Muestra", y="Porcentaje", color=nivel,
            title=f"Abundancia relativa por {nivel} (Top 10 + Otros)",
            labels={"Porcentaje":"% abundancia relativa"}
        )
        fig.update_layout(barmode="stack", xaxis_title="Muestra", yaxis_title="% abundancia relativa")
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning(f"No se encuentra la columna '{nivel}' en el archivo de taxonomía.")
