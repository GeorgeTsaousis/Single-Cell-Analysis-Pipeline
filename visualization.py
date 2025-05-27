import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.express as px

def run():
    st.title("Visualization")

    if "adata" not in st.session_state:
        st.warning("No data found. Please upload and process data in the Analysis section first.")
        return

    adata = st.session_state["adata"]

    tab1, tab2 = st.tabs(["UMAP", "Differential Expression"])

    with tab1:
        st.subheader("UMAP")
        if "X_umap" in adata.obsm:
            color_by = st.selectbox("Color UMAP by", adata.obs.columns, key="umap_color")
            sc.pl.umap(adata, color=color_by, show=False)
            st.pyplot(plt.gcf())
        else:
            st.info("UMAP not computed yet.")

    with tab2:
        st.subheader("Differential Expression Results")
        if "rank_genes_groups" in adata.uns:
            top_n = st.slider("Top N genes to show", 5, 50, 25)
            sc.pl.rank_genes_groups(adata, n_genes=top_n, sharey=False, show=False)
            st.pyplot(plt.gcf())

            # Volcano Plot
            st.subheader("Volcano Plot")
            group = st.selectbox("Select group", list(adata.uns["rank_genes_groups"]["names"].dtype.names))

            df = pd.DataFrame({
                'names': adata.uns["rank_genes_groups"]["names"][group],
                'pvals': adata.uns["rank_genes_groups"]["pvals"][group],
                'logfoldchanges': adata.uns["rank_genes_groups"]["logfoldchanges"][group]
            })

            df['-log10(pval)'] = -np.log10(df['pvals'] + 1e-10)
            df['significant'] = (df['pvals'] < 0.05) & (abs(df['logfoldchanges']) > 1)

            fig = px.scatter(
                df,
                x="logfoldchanges",
                y="-log10(pval)",
                hover_name="names",
                color="significant",
                color_discrete_map={True: "red", False: "grey"},
                title="Volcano Plot"
            )
            fig.update_layout(xaxis_title="Log2 Fold Change", yaxis_title="-log10(p-value)")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No DE results found. Run Differential Expression in Analysis section.")