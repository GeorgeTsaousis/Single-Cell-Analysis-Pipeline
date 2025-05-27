import streamlit as st
import scanpy as sc

def run():
    st.title("Analysis Pipeline")

    data_file = st.file_uploader("Upload .h5ad file", type=["h5ad"])
    if data_file:
        adata = sc.read_h5ad(data_file)
        st.success("Data loaded successfully!")

        st.subheader("Data Summary")
        st.write("Shape:", adata.shape)
        st.dataframe(adata.obs.head())

        # Clear the anndata set from ERCC and MT
        adata = adata[:, [g for g in adata.var_names if not str(g).startswith(('ERCC', 'MT-', 'mt-'))]]

        # Raw now has not ERCC and MT to intefere to the analysis bellow
        adata.raw = adata.copy()

        # Save to session in order to be used in all tabs
        st.session_state["adata"] = adata

        # Preprocessing Parameters 
        st.subheader("Preprocessing Parameters")
        min_genes = st.slider("Minimum genes per cell", 0, 1000, 600)
        min_cells = st.slider("Minimum cells per gene", 0, 50, 3)
        target_sum = st.number_input("Normalization target sum", value=1e4)

        use_harmony = st.checkbox("Use Harmony Integration (Harmony)", value=False)

        if st.button("Run Preprocessing"):
            adata = st.session_state["adata"]

            with st.spinner("Processing..."):
                sc.pp.filter_cells(adata, min_genes=min_genes)
                sc.pp.filter_genes(adata, min_cells=min_cells)

                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)
                sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                adata = adata[:, adata.var.highly_variable]
                sc.pp.scale(adata, max_value=10)
                sc.pp.pca(adata)

                if use_harmony and "batch" in adata.obs.columns:
                    from scanpy.external.pp import harmony_integrate
                    harmony_integrate(adata, 'batch')
                    st.success("Harmony integration completed.")
                elif use_harmony:
                    st.warning("No 'batch' column found in .obs — skipping Harmony.")

                sc.pp.neighbors(adata)
                sc.tl.umap(adata)

                st.success("Preprocessing complete!")
                st.session_state["adata"] = adata

        # Differential Expression Section 
        st.subheader("Differential Expression Analysis")

        # filtering collumns
        adata = st.session_state["adata"]
        groupable_cols = [col for col in adata.obs.columns if adata.obs[col].dtype.name in ['category', 'object']]

        groupby_col = st.selectbox("Group by (for DE analysis)", groupable_cols)

        if st.button("Run Differential Expression Analysis"):
            if groupby_col:
                with st.spinner(f"Running DE analysis using '{groupby_col}'..."):
                    sc.tl.rank_genes_groups(adata, groupby_col, method='wilcoxon')
                st.success("Differential Expression Analysis completed.")
                st.session_state["adata"] = adata
            else:
                st.error("Please select a valid group column for DEA.")
    else:
        st.info("Please upload a .h5ad file to begin analysis.")