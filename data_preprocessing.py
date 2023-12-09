# Importing required libraries
import requests
import pandas as pd
import numpy as np

def fetch_data(protein_list):
    proteins = '%0d'.join(protein_list)
    url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606'
    r = requests.get(url)
    lines = r.text.split('\n')
    data = [l.split('\t') for l in lines]
    df = pd.DataFrame(data[1:-1], columns=data[0])
    return df

def process_data(df):
    interactions = df[['preferredName_A', 'preferredName_B', 'score']]
    return interactions

def create_networkx_graph(interactions):
    G = nx.Graph(name='Protein Interaction Graph')
    interactions = np.array(interactions)
    for i in range(len(interactions)):
        interaction = interactions[i]
        a = interaction[0]
        b = interaction[1]
        w = float(interaction[2])
        G.add_weighted_edges_from([(a, b, w)])
    return G
