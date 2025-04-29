import streamlit as st
import scanpy as sc

def run():
    st.title("📊 Visualization")

    if "adata" not in st.session_state:
        st.warning("No data found. Please run the analysis first.")
        return

    adata = st.session_state["adata"]

    st.subheader("UMAP Plot")
    obs_keys = list(adata.obs_keys())
    color_by = st.selectbox("Color by", obs_keys, index=0)

    if st.button("Show UMAP"):
        with st.spinner("Rendering UMAP..."):
            fig = sc.pl.umap(adata, color=[color_by], return_fig=True)
            st.pyplot(fig)

    st.subheader("Metadata Table")
    st.dataframe(adata.obs.head())
