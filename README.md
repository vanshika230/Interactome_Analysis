# Interactome_Analysis
The Protein-Protein Interaction Network Bipartite Centrality Analysis using NetworkX Notebook aims to analyze and visualize protein-protein interaction networks using bipartite centrality measures. 

Protein interaction data plays a crucial role in comprehending the intricate relationships between gene-encoded biomolecules, enabling us to gain insight into cellular function and predict therapeutic possibilities. Understanding the structure and function of these interaction networks is essential for advancing our knowledge of complex biological systems. PPINs are complex and highly interconnected, with thousands of nodes and edges. Understanding the structure and function of PPINs is a challenging task due to their complexity. Centrality measures can be used to identify the most important nodes in a PPIN and provide insights into network structure and function.
Several databases store protein interaction data, with STRING being one of the most prominent. At present, STRING encompasses over 3.1 billion interactions, involving 20 million proteins across 5,000 organisms, thereby offering a comprehensive resource for exploring protein-protein interactions. We will work on selected human proteins. 

Centrality measures are quantitative measures that capture different aspects of network structure and importance. The most commonly used centrality measures are degree centrality, betweenness centrality, closeness centrality, and eigenvector centrality. Each of these measures provides different insights into network structure and can be used to identify different types of important nodes.This analysis provides insight into the network structure and helps to identify key players in protein-protein interaction networks. You can find more measures provided by NetworkX here :- https://networkx.org/documentation/latest/reference/algorithms/centrality.html 

NetworkX is a popular Python library for network analysis and visualization, offering a range of algorithms and tools for exploring complex networks. The notebook employs NetworkX to construct and analyze protein-protein interaction networks, which are represented as graphs with nodes representing proteins and edges representing interactions between them.

Further work can be extended on this dataset with finding minimum spanning trees, applying other centrality measures such as pagerank algorithm and Leidenalg. We can also use this example to learn isomorphism. The possibilities of analyzing interactomes with NetworkX are endless.

# NOTES
Please note that the ppin.md and jupyter notebook was created 8 months ago as the NetworkX Guides are created as Jupyter Markdown Files so this project is old. This project was deployed in streamlit for the purpose of easy usage and for demonstrations purposes in the open source community.The code has been added from a private repository here!

# Features
Dynamic Centrality Analysis: Users can select from a variety of centrality measures, including Degree Centrality, Eigenvector Centrality, Closeness Centrality, Information Centrality, Betweenness Centrality, and more.

Graph Visualization: The results of the centrality analysis are visually represented through an interactive graph, providing an intuitive understanding of protein interactions and their significance.

Customized Results and Conclusions: Detailed insights and conclusions are provided based on the selected centrality measure, offering users a deeper understanding of the network structure.
# The Web App 
![Front_Page1](https://github.com/vanshika230/Interactome_Analysis/assets/74042272/2685f0d0-5754-4838-ab79-a1731b081e0d)
![Front_Page2](https://github.com/vanshika230/Interactome_Analysis/assets/74042272/d9072059-ba12-475e-aaf5-b8abf03377fe)
![Results](https://github.com/vanshika230/Interactome_Analysis/assets/74042272/0dd836cb-111d-4688-9149-eca5c75acf3c)


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
                    
Here's a brief explanation of each file:

/protein_interaction_analysis/src/__init__.py: An empty file that makes the src directory a Python package.

/protein_interaction_analysis/src/app.py: Contains the Streamlit web application code.

/protein_interaction_analysis/src/network_analysis.py: Contains functions related to network analysis, including fetching data, creating a NetworkX graph, and calculating centrality measures.

/protein_interaction_analysis/src/data_processing.py: Contains functions related to data processing, such as handling API responses and converting data into a pandas DataFrame.

/protein_interaction_analysis/tests/test_network_analysis.py: Contains test cases for the functions defined in network_analysis.py.

/protein_interaction_analysis/requirements.txt: A file containing the required Python packages and their versions.

# Dependencies

NetworkX 2.7.3

Requests 2.26.0

Pandas 1.3.3

NumPy 1.21.2

Matplotlib 3.4.3

Streamlit 1.14.0

# Contributing
If you would like to contribute to the project, please follow the guidelines outlined in CONTRIBUTING.md.

# Acknowledgments
The project relies on open-source libraries and tools, including NetworkX, Streamlit, and others. See ACKNOWLEDGMENTS.md for details.
