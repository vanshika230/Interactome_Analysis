# Importing required libraries
import networkx as nx
import pytest
from network_analysis import calculate_centrality_measure, visualize_network
from data_processing import fetch_data, process_data, create_networkx_graph

# Sample data for testing
protein_list = ['TPH1', 'COMT', 'SLC18A2', 'HTR1B', 'HTR2C']

def test_calculate_centrality_measure():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1)])
    
    degree_centrality = calculate_centrality_measure(G, "Degree Centrality")
    assert degree_centrality == {1: 0.6666666666666666, 2: 1.0, 3: 0.6666666666666666}

    eigenvector_centrality = calculate_centrality_measure(G, "Eigenvector Centrality")
    assert eigenvector_centrality == {1: 0.453182105296043, 2: 0.7663649885037259, 3: 0.453182105296043}

def test_visualize_network():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1)])
    centrality_scores = {1: 0.6666666666666666, 2: 1.0, 3: 0.6666666666666666}
    visualize_network(G, centrality_scores, "Degree Centrality")

def test_fetch_data():
    # Mock the requests.get function for testing purposes
    class MockResponse:
        @staticmethod
        def text():
            return 'ID1\tID2\tScore\n1\t2\t0.8\n2\t3\t0.7\n3\t1\t0.9\n'
    
    with monkeypatch.context() as m:
        m.setattr(requests, 'get', lambda *args, **kwargs: MockResponse)
        df = fetch_data(protein_list)
    
    expected_columns = ['ID1', 'ID2', 'Score']
    assert list(df.columns) == expected_columns
    assert len(df) == 3

def test_process_data():
    # Assuming df is a DataFrame with 'ID1', 'ID2', and 'Score' columns
    df = pd.DataFrame({'ID1': [1, 2, 3], 'ID2': [2, 3, 1], 'Score': [0.8, 0.7, 0.9]})
    interactions = process_data(df)
    expected_columns = ['ID1', 'ID2', 'Score']
    assert list(interactions.columns) == expected_columns
    assert len(interactions) == 3

def test_create_networkx_graph():
    interactions = pd.DataFrame({'ID1': [1, 2, 3], 'ID2': [2, 3, 1], 'Score': [0.8, 0.7, 0.9]})
    G = create_networkx_graph(interactions)
    assert isinstance(G, nx.Graph)
    assert len(G.nodes) == 3
    assert len(G.edges) == 3