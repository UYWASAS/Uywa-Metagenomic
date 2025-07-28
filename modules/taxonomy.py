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

    st.subheader("Barplot apilado por Filo (top 10)")
    if "Phylum" in taxonomy.columns:
        # Suma por filo
        otus_tax = otus.T.join(taxonomy["Phylum"])
        phylum_sum = otus_tax.groupby("Phylum").sum().T
        top_phyla = phylum_sum.sum().sort_values(ascending=False).head(10).index
        phylum_sum = phylum_sum[top_phyla]
        phylum_sum_norm = phylum_sum.div(phylum_sum.sum(axis=1), axis=0)
        fig = px.bar(phylum_sum_norm, title="Abundancia Relativa por Filo", labels={"value":"Abundancia Relativa","index":"Muestra"})
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("No se encuentra la columna 'Phylum' en el archivo de taxonomía.")
