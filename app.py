# Importing required libraries
import streamlit as st
from network_analysis import calculate_centrality_measure, visualize_network
from data_processing import fetch_data, process_data, create_networkx_graph

# Title and description
st.title("Protein Interaction Network Analysis")
st.write("This web app uses NetworkX Centrality Measures to visualize and analyze a protein interaction network.")

# Fetching Dataset 
protein_list = ['TPH1', 'COMT', 'SLC18A2', 'HTR1B', 'HTR2C']
df = fetch_data(protein_list)
interactions = process_data(df)
G = create_networkx_graph(interactions)

# Sidebar for user input
st.sidebar.header("Select Centrality Measure")
centrality_measure = st.sidebar.radio("Select Centrality Measure:", [
    "Degree Centrality", "Eigenvector Centrality", "Closeness Centrality",
    "Information Centrality", "Betweenness Centrality",
    "Current Flow Betweeenness Centrality", "Communicability Betweenness Centrality",
    "Load Centrality", "Subgraph Centrality", "Harmonic Centrality", "Second Order Centrality"
])

# Calculate Centrality Measure
centrality_scores = calculate_centrality_measure(G, centrality_measure)

# Display Centrality Measure
st.header(f"{centrality_measure} Scores")
st.write("Centrality scores for the entire protein interaction network:")
for protein, score in centrality_scores.items():
    st.write(f"{protein}: {score:.2f}")

# Visualization of Centrality Measure
visualize_network(G, centrality_scores, centrality_measure)

# Results and Conclusion for Centrality Measure
st.header(f"Results and Conclusion for {centrality_measure}")
st.write("Centrality measures for the entire protein interaction network have been calculated and visualized.")
st.write(f"This analysis provides insights into the network structure and the importance of individual proteins based on {centrality_measure}.")
