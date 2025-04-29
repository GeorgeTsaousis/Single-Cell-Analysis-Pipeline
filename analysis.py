import streamlit as st
import scanpy as sc

def run():
    st.title("🔬 Analysis Pipeline")

    data_file = st.file_uploader("Upload .h5ad file", type=["h5ad"])
    if data_file:
        adata = sc.read_h5ad(data_file)
        st.success("Data loaded successfully!")

        st.subheader("Data Summary")
        st.write("Shape:", adata.shape)
        st.dataframe(adata.obs.head())

        st.subheader("Preprocessing Parameters")
        min_genes = st.slider("Minimum genes per cell", 0, 1000, 600)
        min_cells = st.slider("Minimum cells per gene", 0, 50, 3)
        target_sum = st.number_input("Normalization target sum", value=1e4)

        if st.button("Run Preprocessing"):
            with st.spinner("Processing..."):
                sc.pp.filter_cells(adata, min_genes=min_genes)
                sc.pp.filter_genes(adata, min_cells=min_cells)
                adata = adata[:, [g for g in adata.var_names if not str(g).startswith(tuple(['ERCC', 'MT-', 'mt-']))]]
                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)
                sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                adata.raw = adata
                adata = adata[:, adata.var.highly_variable]
                sc.pp.scale(adata, max_value=10)
                sc.pp.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
            st.success("Preprocessing complete!")
            st.session_state["adata"] = adata
    else:
        st.info("Please upload a .h5ad file to begin analysis.")
