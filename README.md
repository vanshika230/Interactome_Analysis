# Interactome_Analysis
The Protein-Protein Interaction Network Bipartite Centrality Analysis using NetworkX Notebook aims to analyze and visualize protein-protein interaction networks using bipartite centrality measures. 

Protein interaction data plays a crucial role in comprehending the intricate relationships between gene-encoded biomolecules, enabling us to gain insight into cellular function and predict therapeutic possibilities. Understanding the structure and function of these interaction networks is essential for advancing our knowledge of complex biological systems. PPINs are complex and highly interconnected, with thousands of nodes and edges. Understanding the structure and function of PPINs is a challenging task due to their complexity. Centrality measures can be used to identify the most important nodes in a PPIN and provide insights into network structure and function.
Several databases store protein interaction data, with STRING being one of the most prominent. At present, STRING encompasses over 3.1 billion interactions, involving 20 million proteins across 5,000 organisms, thereby offering a comprehensive resource for exploring protein-protein interactions. We will work on selected human proteins. 

Centrality measures are quantitative measures that capture different aspects of network structure and importance. The most commonly used centrality measures are degree centrality, betweenness centrality, closeness centrality, and eigenvector centrality. Each of these measures provides different insights into network structure and can be used to identify different types of important nodes.This analysis provides insight into the network structure and helps to identify key players in protein-protein interaction networks. You can find more measures provided by NetworkX here :- https://networkx.org/documentation/latest/reference/algorithms/centrality.html 

NetworkX is a popular Python library for network analysis and visualization, offering a range of algorithms and tools for exploring complex networks. The notebook employs NetworkX to construct and analyze protein-protein interaction networks, which are represented as graphs with nodes representing proteins and edges representing interactions between them.

Further work can be extended on this dataset with finding minimum spanning trees, applying other centrality measures such as pagerank algorithm and Leidenalg. We can also use this example to learn isomorphism. The possibilities of analyzing interactomes with NetworkX are endless.

# Features
Dynamic Centrality Analysis: Users can select from a variety of centrality measures, including Degree Centrality, Eigenvector Centrality, Closeness Centrality, Information Centrality, Betweenness Centrality, and more.

Graph Visualization: The results of the centrality analysis are visually represented through an interactive graph, providing an intuitive understanding of protein interactions and their significance.

Customized Results and Conclusions: Detailed insights and conclusions are provided based on the selected centrality measure, offering users a deeper understanding of the network structure.

# Installation
        To run the Protein Interaction Analysis project locally, follow these steps:

# Clone the Repository:

        open Git bash
        git clone https://github.com/your-username/protein-interaction-analysis.git
        cd protein-interaction-analysis

# Install Dependencies:
        pip install -r requirements.txt

# Run the Web App:


        streamlit run src/app.py

The web app should now be accessible in your web browser at http://localhost:8501.

# Usage
- Select Centrality Measure:
- Choose from various centrality measures in the sidebar.

# Explore Results:
- View centrality scores for the entire protein interaction network, along with a visual representation of the graph.

- Customized Conclusions:
Read detailed results and conclusions specific to the selected centrality measure.

# Directory Structure
/protein_interaction_analysis
    |-- src
    |    |-- __init__.py
    |    |-- app.py
    |    |-- network_analysis.py
    |    |-- data_processing.py
    |
    |-- tests
    |    |-- test_network_analysis.py
    |
    |-- requirements.txt
    |-- README.md
src/: Contains the main source code files.
tests/: Houses unit tests for the project.
requirements.txt: Lists project dependencies.

# Dependencies
NetworkX 2.7.3
Requests 2.26.0
Pandas 1.3.3
NumPy 1.21.2
Matplotlib 3.4.3
Streamlit 1.14.0

# Contributing
If you would like to contribute to the project, please follow the guidelines outlined in CONTRIBUTING.md.

# License
This project is licensed under the MIT License.

# Acknowledgments
The project relies on open-source libraries and tools, including NetworkX, Streamlit, and others. See ACKNOWLEDGMENTS.md for details.