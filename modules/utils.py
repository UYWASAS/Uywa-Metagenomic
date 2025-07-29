import pandas as pd

def load_table(file, index_col=None, sep=None):
    """
    Carga un archivo de tabla (.csv, .tsv, .xlsx) y configura el índice si se indica.
    Ahora es robusto a variantes comunes de nombres de columna y espacios.
    - file: archivo cargado (st.file_uploader)
    - index_col: nombre de la columna a usar como índice (ej: "OTU", "SampleID")
    - sep: separador opcional (por defecto autodetecta por extensión)
    """
    if file is None:
        return None
    filename = file.name.lower()
    if filename.endswith(".csv"):
        df = pd.read_csv(file, sep=sep if sep else ",")
    elif filename.endswith(".tsv") or filename.endswith(".txt"):
        df = pd.read_csv(file, sep=sep if sep else "\t")
    elif filename.endswith(".xlsx"):
        df = pd.read_excel(file)
    else:
        raise ValueError("Formato de archivo no soportado.")
    # Normaliza columnas: quita espacios y pone en minúsculas para la búsqueda
    colnames_raw = list(df.columns)
    colnames_norm = [str(col).strip().replace(" ", "").lower() for col in df.columns]
    if index_col is not None:
        idx_norm = index_col.strip().replace(" ", "").lower()
        # Busca una columna que coincida con el index_col normalizado
        if idx_norm in colnames_norm:
            true_col = colnames_raw[colnames_norm.index(idx_norm)]
            df.set_index(true_col, inplace=True)
        else:
            raise ValueError(
                f"La columna '{index_col}' no se encuentra en el archivo {filename}\n"
                f"Las columnas encontradas son: {colnames_raw}"
            )
    df.columns = [str(col).strip() for col in df.columns]
    df.index = df.index.map(lambda x: str(x).strip())
    return df

def safe_float(val, default=0.0):
    try:
        if isinstance(val, str):
            val = val.replace(",", ".")
        return float(val)
    except Exception:
        return default

def clean_state(keys_prefix, valid_names):
    import streamlit as st
    for key in list(st.session_state.keys()):
        for prefix in keys_prefix:
            if key.startswith(prefix):
                found = False
                for n in valid_names:
                    if key.endswith(f"{n}_incl_input") or key.endswith(f"{n}_input"):
                        found = True
                        break
                if not found:
                    del st.session_state[key]
