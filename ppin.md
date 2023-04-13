---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Introduction
Protein interaction data plays a crucial role in comprehending the intricate relationships between gene-encoded biomolecules, enabling us to gain insight into cellular function and predict therapeutic possibilities. Understanding the structure and function of these interaction networks is essential for advancing our knowledge of complex biological systems. PPINs are complex and highly interconnected, with thousands of nodes and edges. Understanding the structure and function of PPINs is a challenging task due to their complexity. Centrality measures can be used to identify the most important nodes in a PPIN and provide insights into network structure and function.The Protein-Protein Interaction Network Bipartite Centrality Analysis using NetworkX Notebook aims to analyze and visualize protein-protein interaction networks using bipartite centrality measures.
Several databases store protein interaction data, with STRING being one of the most prominent. At present, STRING encompasses over 3.1 billion interactions, involving 20 million proteins across 5,000 organisms, thereby offering a comprehensive resource for exploring protein-protein interactions. We will work on selected human proteins.
Centrality measures are quantitative measures that capture different aspects of network structure and importance. The most commonly used centrality measures are degree centrality, betweenness centrality, closeness centrality, and eigenvector centrality. Each of these measures provides different insights into network structure and can be used to identify different types of important nodes.This analysis provides insight into the network structure and helps to identify key players in protein-protein interaction networks. You can find more measures provided by NetworkX here :- https://networkx.org/documentation/latest/reference/algorithms/centrality.html
NetworkX is a popular Python library for network analysis and visualization, offering a range of algorithms and tools for exploring complex networks. The notebook employs NetworkX to construct and analyze protein-protein interaction networks, which are represented as graphs with nodes representing proteins and edges representing interactions between them.

# Importing the required libraries and packages 

```python
import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
```

## About the Dataset
The STRING database is a biological database and web resource of known and predicted protein interactions, including both direct physical interactions and indirect functional associations. The interactions are derived from multiple sources, including high-throughput experiments, co-expression analysis, and text mining of the scientific literature. The database covers a wide range of organisms, including bacteria, archaea, and eukaryotes.
The STRING database provides users with several tools for visualizing and analyzing protein interaction networks. Users can search for individual proteins and view their known interactions, or upload their own list of proteins to generate a custom interaction network. The database also provides tools for analyzing the functional enrichment of proteins in a network, predicting protein interactions based on sequence similarity, and comparing different protein networks.
Dataset Source :- https://string-db.org/cgi/about?footer_active_subpage=references

# Fetching the Dataset
We will use the STRING database API to retrieve a protein interaction network for a list of human proteins. This code uses the STRING database API to retrieve a protein interaction network for a list of human proteins.We make a list of human protein names.We use the join() method to combine the protein names into a string separated by %0d (a URL-encoded newline character). The resulting string is stored in the proteins variable. A URL is constructed using the proteins string and a species identifier for human (9606). The requests module is used to send a GET request to the URL, and the response is stored in the r variable. The list of proteins was selected from a larger set of human proteins related to serotonin. 9606 is the species number for humans. You can checkout more methods to fetch data from String database here :- https://string-db.org/cgi/help?sessionId=bbzQfhU4eZmR


```python
# Define a list of protein IDs to query the STRING database
protein_list = ['TPH1','COMT','SLC18A2','HTR1B','HTR2C','HTR2A','MAOA',
            'TPH2','HTR1A','HTR7','SLC6A4','GABBR2','POMC','GNAI3',
            'NPY','ADCY1','PDYN','GRM2','GRM3','GABBR1']
proteins = '%0d'.join(protein_list) # Concatenate the protein IDs into a single string using the '%0d' separator
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606' # Construct the URL for the STRING database API endpoint, specifying the protein IDs and species code
r = requests.get(url)
```

```python
lines = r.text.split('\n') # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines] # split each line into its components based on tabs
# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0]) # Send a GET request to the API endpoint and store the response in a variable
```

```python
df.head(n=10)
```

```python
# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df[['preferredName_A', 'preferredName_B', 'score']] 
```

```python
interactions[:10] 
```

# Construction of Graph using NetworkX
Now we will construct a protein interaction graph using the NetworkX library.We create a new empty graph G is created using the NetworkX Graph class, with a name attribute set to 'Protein Interaction Graph'.
The code iterates through each row of the interactions array. For each row, it extracts the two proteins involved in the interaction (a and b) and the weight of the interaction (w). The weight is stored as a float, and represents the strength of the interaction, with higher scores indicating weaker interactions.
Overall, this code constructs a protein interaction graph using the NetworkX library, where nodes represent proteins and edges represent interactions between proteins. The edges are weighted based on the strength of the interaction, with higher weights indicating weaker interactions.





```python
G=nx.Graph(name='Protein Interaction Graph')
interactions = np.array(interactions)
for i in range(len(interactions)):
    interaction = interactions[i]
    a = interaction[0] # protein a node which will be preferredName_A
    b = interaction[1] # protein b node which will be preferredName_B
    w = float(interaction[2]) # score as weighted edge where high scores = low weight
    G.add_weighted_edges_from([(a,b,w)]) # add weighted edge to graph
```

```python
print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())
```

# Visualizing Network

```python
pos = nx.spring_layout(G) # position the nodes using the spring layout
plt.figure(figsize=(11, 11), facecolor='white')
nx.draw_networkx(G, node_color='blue', edge_color='grey')
plt.axis('off')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/network.jpg))

# Finding Centrality Measures for the Protein-Protein Network
Finding centrality measures for a protein-protein interaction (PPI) network dataset can give various insights on the protein-protein network.
1) Identifying highly connected nodes: Centrality measures can help identify highly connected nodes in the PPI network, which are often key regulators or hubs of the network. 

2) Predicting essential genes: We can identify essential genes in the PPI network, which are genes whose deletion or knockout would result in a lethal phenotype. Highly central genes are often essential genes, as they are required for maintaining the integrity of the network.

3) Identifying disease-associated genes: Centrality measures can help identify genes that are associated with disease or other phenotypic traits. Disease-associated genes are often highly central in PPI networks, as they are involved in key pathways or processes that are disrupted in disease.

4) Identifying important proteins: Centrality measures can be used to identify proteins that are most important in a PPI network. For example, proteins with high degree centrality are likely to be involved in many interactions and may play important roles in regulating cellular processes. Identifying these important proteins can help researchers better understand the underlying biology of a system.

5) Characterizing network structure: Centrality measures can also be used to characterize the overall structure of a PPI network. For example, networks with high clustering coefficient and low betweenness centrality are indicative of modular networks, where proteins tend to interact within well-defined modules. Understanding the network structure can provide insights into the functional organization of a biological system.

Overall, finding centrality measures for a PPI network dataset can provide valuable insights into the organization and function of the network, as well as its role in disease and other biological processes. These insights can help guide future research and drug development efforts.


# Degree Centrality 
The degree centrality for a node v is the fraction of nodes it is connected to. We find the degree centrality of our network to identify the most connected nodes.

```python
dc=nx.degree_centrality(G)
dc
```

## Visualization of Degree Centrality 
The nodes with higher degree centrality are colored with warmer colors (e.g., red), while nodes with lower degree centrality are colored with cooler colors (e.g., blue). We can observe that the visualization correlates with the conclusions drawn above.

```python
# Create a list of node colors based on their Degree centrality
colors = [dc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Degree centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Degree Centrality')
plt.show()

```
![DegreeCentrality](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/DegreeCentrality.jpg))

## Results
From the above calculations we can infer that HTR1A, HTR2A, and SLC6A4 have the highest degree centrality, indicating that they are highly connected to other proteins in the network and are likely to be important hubs in the network, and their disruption may have significant effects on the overall network structure and function.'HTR7','TPH2','ADCY1' have the lowest degree centrality indicating that they aren't well connected to the network. 
# Eigen Vector Centrality
Eigenvector centrality measures the influence of a node in the network based on the connections it has to other highly connected nodes in the network. The higher the eigenvector centrality of a node, the more influential it is in the network.

```python
ec = nx.eigenvector_centrality(G)
```

```python
sorted((v, f"{c:0.2f}") for v, c in ec.items())
```

## Lets Visualize using Eigenvector Centrality
The nodes with higher eigenvector centrality are colored with warmer colors (e.g., red), while nodes with lower eigenvector centrality are colored with cooler colors (e.g., blue).  We can observe that the visualization correlates with the conclusions drawn above.

```python
# Create a list of node colors based on their eigenvector centrality
colors = [ec[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their eigenvector centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Eigenvector Centrality')
plt.show()
```
![EigenVector](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/EigenVector.jpg))

## Results
The results suggest that HTR1A, SLC6A4, and COMT are the most influential nodes in the network, as they have the highest eigenvector centrality scores. This indicates that they are connected to other highly connected nodes in the network and may play important roles in regulating the interactions among other proteins in the network.These proteins may be good candidates for further investigation as potential targets for drug development or as biomarkers for certain diseases. On the other hand, ADCY1 has the lowest eigenvector centrality score, indicating that it may be less influential in the network.


# Closeness Centrality
Closeness centrality of a node u is the reciprocal of the average shortest path distance to u over all n-1 reachable nodes. A high closeness centrality score indicates that a node is close to many other nodes and can easily interact with them.

```python
cc= nx.closeness_centrality(G)
cc
```

## Lets Visualize using Closeness Centrality
The nodes with higher Closeness centrality are colored with warmer colors (e.g., red), while nodes with lower closeness centrality are colored with cooler colors (e.g., blue). We can observe that the visualization correlates with the conclusions drawn above.

```python
# Create a list of node colors based on their closeness centrality
colors = [cc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their closeness centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Closeness Centrality')
plt.show()
```
![Closeness](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/Closeness.jpg))
## Results
The nodes with the highest closeness centrality values are HTR1A with a value of 0.8636, followed by TPH1, COMT, and HTR2A with values of 0.76. These nodes are, therefore, the most central nodes in the graph according to this metric.

On the other hand, the nodes with the lowest closeness centrality values are ADCY1 with a value of 0.4872, followed by GABBR2 and HTR7 with values of 0.5938 and 0.5758 respectively. These nodes are, therefore, the least central nodes in the graph according to this metric.
# Information Centrality 
The information centrality is a measure that quantifies the influence of a node in a network based on the amount of information that passes through it. In a PPIN, this could be interpreted as the ability of a protein to transmit or receive signals from other proteins in the network.

```python
ic=nx.information_centrality(G)
ic
```

## Lets Visualize using Information Centrality
The nodes with higher information centrality are colored with warmer colors (e.g., red), while nodes with lower information centrality are colored with cooler colors (e.g., blue).We find that visualization correlates with the results above. 

```python
# Create a list of node colors based on their information centrality
colors = [ic[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their information centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Information Centrality')
plt.show()
```
![InformationCentrality](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/InformationCentrality.jpg))
## Results
TPH1, COMT, and HTR1A have relatively high information centrality values, suggesting that they are important in the flow of information in the PPIN. On the other hand, TPH2, ADCY1,GABBR2 has a low information centrality value, indicating that it may have less influence in the flow of information in the network.
# Betweenness Centrality 
Betweenness centrality is a measure of the importance of a node in a network based on the number of shortest paths that pass through it.Therefore, a node with high betweenness centrality is important in maintaining the connectivity of the network. Nodes with high betweenness centrality may be bottlenecks or bridges that connect different parts of the network. In contrast, nodes with low betweenness centrality are less important for maintaining the overall connectivity of the network.

```python
bc= nx.betweenness_centrality(G)
bc
```

## Lets Visualize using Betweenness Centrality
The nodes with higher betweenness centrality are colored with warmer colors (e.g., red), while nodes with lower betweenness centrality are colored with cooler colors (e.g., blue).The visualization also correlates with the results above. 

```python
# Create a list of node colors based on their Betweenness centrality
colors = [bc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Betweenness centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Betweenness Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/betweenness.jpg))
## Result
In the given result, some nodes have relatively high betweenness centrality, such as 'HTR1A', 'HTR1B', 'GRM2', 'GNAI3', and 'HTR2A', which suggests that these nodes are important for maintaining the connectivity of the network. On the other hand, some nodes have low betweenness centrality, such as 'ADCY1','GABBR2' which suggests that these nodes are less important for maintaining the overall connectivity of the network.
# Current Flow Betweeenness Centrality
In network analysis, current flow betweenness centrality is a measure of a node's importance in a graph based on the current flow between all pairs of nodes. The current flow betweenness centrality of a node is calculated as the sum of the current flows between all pairs of nodes passing through that node.The nodes with higher current flow betweenness centrality values are considered to be more important in the graph, as they play a crucial role in the flow of current between nodes in the network.

```python
cfb=nx.current_flow_betweenness_centrality(G)
cfb
```

## Lets Visualize using Current Flow Betweenness Centrality
The nodes with higher current flow betweenness centrality are colored with warmer colors (e.g., red), while nodes with lower current flow betweenness centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their current flow betweenness centrality
colors = [cfb[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their current flow betweenness centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Current Flow Betweenness Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/currentflowbetween.jpg))

## Results
In this specific result, the nodes with the highest current flow betweenness centrality are HTR1A, HTR1B, GRM2, and HTR2A, which indicates that these nodes are particularly important for maintaining the connectivity and information flow in the network. These nodes may be potential targets for interventions to modulate the network's function or prevent the spread of information through the network.

On the other hand, nodes with low current flow betweenness centrality, such as ADCY1,HTR7 may have a smaller impact on the network's overall connectivity and function.




# Communicability Betweenness Centrality
The score for each node indicates the extent to which that node mediates communication between other nodes in the network. A higher score indicates that the node is more important for communication between other nodes.

```python
cbc=nx.communicability_betweenness_centrality(G)
cbc
```

## Lets Visualize using Communicability Betweenness Centrality
The nodes with higher  Communicability Betweenness Centrality are colored with warmer colors (e.g., red), while nodes with lower Communicability Betweenness Centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their Communicability Betweenness Centrality
colors = [cbc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their  Communicability Betweenness Centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Communicability Betweenness Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/communicability.jpg))
## Results
The node 'HTR1A' has the highest score of 0.676, indicating that it is more important for communication between other nodes than any other node in the network. On the other hand, the node 'GABBR2' has a score of 0.119, which suggests that it is a relatively peripheral node in terms of communication flow in the network.
# Load Centrality 
Load centrality measures the importance of a node based on the number of shortest paths that pass through it. In this case, the higher the load centrality value of a node, the more important it is in maintaining the network's connectivity.

```python
lc=nx.load_centrality(G)
lc
```

## Lets Visualize using Load Centrality
The nodes with higher Load centrality are colored with warmer colors (e.g., red), while nodes with lower Load centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their Load Centrality
colors = [lc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Load Centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Load Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/load.jpg))
## Results
The nodes "HTR1A" and "HTR1B" have the highest load centrality values of 0.086 and 0.073, respectively, suggesting that these nodes are critical for the network's connectivity. Other nodes with relatively high load centrality values include "GRM2," "GNAI3," "HTR2A," and "NPY," indicating that these nodes are also important in maintaining network connectivity.

On the other hand, some nodes have lower load centrality values, such as "ADCY1" with a value of 0.0, which suggests that these nodes may be less important in maintaining network connectivity.




# Subgraph Centrality 
Subgraph centrality measures the importance of a node based on its participation in all the subgraphs of the network. In other words, it reflects how many closed paths, or cycles, a node participates in.

```python
sc= nx.subgraph_centrality(G)
sc
```

## Lets Visualize using Subgraph Centrality
The nodes with higher Subgraph Centrality are colored with warmer colors (e.g., red), while nodes with lower Subgraph Centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their Subgraph Centrality
colors = [sc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Subgraph Centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Subgraph Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/subgraph.jpg))
## Results
The nodes "HTR1A," "HTR2A," and "COMT" have the highest subgraph centrality values of 7753, 6855, and 6585, respectively, indicating that these nodes are the most important in terms of their participation in the subgraphs of the network. Other nodes with relatively high subgraph centrality values include "SLC6A4," "TPH1," "NPY," "HTR2C," and "TPH2," suggesting that these nodes are also important in participating in the subgraphs of the network.

On the other hand, some nodes have lower subgraph centrality values, such as "ADCY1" with a value of 87.3, which suggests that these nodes participate in fewer subgraphs and may be less important in terms of their contribution to the overall network structure.

# Estrada Index 
The Estrada index is a measure of the graph complexity that takes into account the contribution of all nodes and their interactions. It is calculated by exponentiating the adjacency matrix of the network and summing all of its elements.

A higher Estrada index value suggests a higher level of complexity and a greater number of interactions between nodes in the network. In other words, the more interactions and interdependencies between nodes, the higher the Estrada index value.

```python
nx.estrada_index(G)
```

## Results
The Estrada index value of 80826.4 indicates that the ppin network is relatively complex and has many interactions between its nodes. This suggests that the network may be biologically important, as biological systems are often characterized by complex interactions between their components.
# Harmonic Centrality
The harmonic centrality measures the sum of the inverse of the shortest path lengths between a node and all other nodes in the network. A higher harmonic centrality value suggests that a node is more centrally located within the network and has more influence over the flow of information within the network.

```python
hc=nx.harmonic_centrality(G)
hc
```

## Lets Visualize using Harmonic Centrality
The nodes with higher harmonic centrality are colored with warmer colors (e.g., red), while nodes with lower Harmonic Centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their Harmonic Centrality
colors = [hc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Harmonic Centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Harmonic Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/harmonic.jpg))
## Results
The nodes with the highest harmonic centrality values are HTR1A (17.5), COMT (16.0), and TPH1 (16.0), indicating that these nodes are more central and influential in the network. The nodes with the lowest harmonic centrality values are ADCY1 (10.0), GABBR2 (12.5), and HTR7 (12.3), indicating that these nodes are less central and influential in the network.
# Second Order Centrality 
The second-order centrality measures the influence of a node based on its indirect connections, which means how many nodes are reachable through a node's immediate neighbors.

```python
soc=nx.second_order_centrality(G)
soc
```

## Lets Visualize using Second Order Centrality
The nodes with higher Second Order centrality are colored with warmer colors (e.g., red), while nodes with lower Second Order Centrality are colored with cooler colors (e.g., blue).

```python
# Create a list of node colors based on their Second Order Centrality
colors = [soc[node] for node in G.nodes()]

# Draw the graph using spring layout and node colors based on their Second Order Centrality
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors, cmap='coolwarm')
nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
plt.axis('off')
plt.title('Protein Interaction Graph with Second Order Centrality')
plt.show()
```
![Network](https://github.com/vanshika230/Interactome_Analysis/blob/main/Images/secondorder.png)
## Results
The node with the highest second-order centrality is ADCY1 with a value of 73.88, which means that it can reach a lot of nodes indirectly through its neighbors. The next highest values are for HTR7 and GABBR2, which also have high values, indicating their indirect influence on the network. On the other hand, HTR1A and HTR2A have low second-order centrality values, indicating their limited indirect influence on the network.
# Voterank 
The nx.voterank(G) function returns a list of nodes ranked by their influence in the network using a voter model approach. In this case, the nodes are ordered from most to least influential.

```python
nx.voterank(G)
```

## Results
The most influential node in the network is HTR1A, followed by SLC6A4 and HTR2A. These nodes have the highest degree and centrality measures in the network. Other highly influential nodes include TPH1, GNAI3, COMT, NPY, GRM2, MAOA, and HTR1B.


## Conclusion 
In conclusion, the use of NetworkX and centrality measures provides a powerful and flexible framework for analyzing protein-protein interaction networks (PPINs). NetworkX offers a range of functions and tools for constructing and visualizing networks, computing different centrality measures, and analyzing network properties and structures.

Through the use of NetworkX and centrality measures, we can gain a deeper understanding of the complex organization and function of PPINs. For instance, by analyzing the degree, betweenness, closeness, and eigenvector centrality measures of nodes in a PPIN, we can identify hub proteins, bottlenecks, and other types of important nodes that play critical roles in maintaining the network's stability and function.

Moreover, NetworkX allows us to easily visualize and interpret the results of our analyses, making it easier to communicate our findings to other researchers and stakeholders. With these capabilities, NetworkX and centrality measures offer a powerful and valuable toolset for studying PPINs and other types of complex networks.

We have looked through various centrality measures applied on String Dataset using NetworkX library. Now it is ypur turn to play around with it. Choose different set of proteins available on the String Website and make your own networks. 
