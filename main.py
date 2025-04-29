import streamlit as st
import analysis
import visualization
import team_info

st.set_page_config(page_title="scRNA-seq Analyzer", layout="wide")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Analysis", "Visualization", "Team Info"])

# Load corresponding page
if page == "Analysis":
    analysis.run()
elif page == "Visualization":
    visualization.run()
elif page == "Team Info":
    team_info.run()
