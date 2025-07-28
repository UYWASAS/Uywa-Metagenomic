import pandas as pd

def load_table(file):
    if file is None:
        return None
    if file.name.endswith("csv"):
        return pd.read_csv(file, index_col=0)
    elif file.name.endswith("tsv"):
        return pd.read_csv(file, index_col=0, sep="\t")
    else:
        return pd.read_excel(file, index_col=0)

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
