# Importing required libraries
import networkx as nx
import streamlit as st
import matplotlib.pyplot as plt

def calculate_centrality_measure(G, centrality_measure):
    if centrality_measure == "Degree Centrality":
        centrality_scores = nx.degree_centrality(G)
    elif centrality_measure == "Eigenvector Centrality":
        centrality_scores = nx.eigenvector_centrality(G)
    elif centrality_measure == "Closeness Centrality":
        centrality_scores = nx.closeness_centrality(G)
    elif centrality_measure == "Information Centrality":
        centrality_scores = nx.information_centrality(G)
    elif centrality_measure == "Betweenness Centrality":
        centrality_scores = nx.betweenness_centrality(G)
    elif centrality_measure == "Current Flow Betweeenness Centrality":
        centrality_scores = nx.current_flow_betweenness_centrality(G)
    elif centrality_measure == "Communicability Betweenness Centrality":
        centrality_scores = nx.communicability_betweenness_centrality(G)
    elif centrality_measure == "Load Centrality":
        centrality_scores = nx.load_centrality(G)
    elif centrality_measure == "Subgraph Centrality":
        centrality_scores = nx.subgraph_centrality(G)
    elif centrality_measure == "Harmonic Centrality":
        centrality_scores = nx.harmonic_centrality(G)
    elif centrality_measure == "Second Order Centrality":
        centrality_scores = nx.second_order_centrality(G)
    else:
        st.error("Invalid Centrality Measure selected.")
        st.stop()

    return centrality_scores

def visualize_network(G, centrality_scores, centrality_measure):
    colors_centrality = [centrality_scores[node] for node in G.nodes()]

    # Draw the graph using spring layout and node colors based on Centrality Measure
    plt.figure(figsize=(12, 12))
    pos_centrality = nx.spring_layout(G, k=0.5, iterations=50)
    nx.draw_networkx_edges(G, pos_centrality, alpha=0.2)
    nx.draw_networkx_nodes(G, pos_centrality, node_size=100, node_color=colors_centrality, cmap='coolwarm')
    nx.draw_networkx_labels(G, pos_centrality, font_size=10, font_family='sans-serif')
    plt.axis('off')
    plt.title(f'Protein Interaction Graph with {centrality_measure}')
    
    # Display the plot in the Streamlit app
    st.pyplot()
