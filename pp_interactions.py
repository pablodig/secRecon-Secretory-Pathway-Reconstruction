import requests
import json
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm

def fetch_string_interactions(gene_list):
    """
    Fetch protein-protein interactions from STRING database for a list of genes.
    
    Parameters:
    - gene_list (list): List of gene names to fetch interactions for.
    
    Returns:
    - G (networkx.Graph): A graph where nodes are proteins and edges represent interactions.
    """
    # Initialize an empty graph
    G = nx.Graph()
    
    # Fetch interactions for each gene in the list
    for gene in gene_list:
        url = f"https://string-db.org/api/json/network?identifiers={gene}"
        response = requests.get(url)
        
        if response.status_code == 200:
            interactions = json.loads(response.text)
            
            for interaction in interactions:
                protein1 = interaction['preferredName_A']
                protein2 = interaction['preferredName_B']
                score = interaction['score']
                
                # Add nodes and edges to the graph
                G.add_node(protein1)
                G.add_node(protein2)
                G.add_edge(protein1, protein2, weight=score)
        else:
            print(f'No interaction found for {gene}')
                
    return G

def fetch_biogrid_interactions(gene_list):
    """
    Fetch protein-protein interactions from BioGRID database for a list of genes.
    
    Parameters:
    - gene_list (list): List of gene names to fetch interactions for.
    
    Returns:
    - G (networkx.Graph): A graph where nodes are proteins and edges represent interactions.
    """
    # Initialize an empty graph
    G = nx.Graph()
    
    # Fetch interactions for each gene in the list
    url = f"https://webservice.thebiogrid.org/interactions/?accesskey=0659585a790acc50563e75a55f72b15f&format=json"
    response = requests.get(url)
    
    if response.status_code == 200:
        interactions = response.json()
        
        for gene in gene_list:
            for interaction_id, interaction in interactions.items():
                # Filter interactions to only include those involving organism ID 9606
                if interaction['ORGANISM_A'] == 9606 and interaction['ORGANISM_B'] == 9606:
                    # Further filter to only include interactions involving the gene of interest
                    if interaction['OFFICIAL_SYMBOL_A'] == gene or interaction['OFFICIAL_SYMBOL_B'] == gene:
                        protein1 = interaction['OFFICIAL_SYMBOL_A']
                        protein2 = interaction['OFFICIAL_SYMBOL_B']
                        score = interaction.get('SCORE', 1)  # Use 1 as default if no score is provided
                        
                        # Add nodes and edges to the graph
                        G.add_node(protein1)
                        G.add_node(protein2)
                        G.add_edge(protein1, protein2, weight=score)
    else:
        print(f'No interaction found for {gene}')
                
    return G

def visualize_network(G, gene_list=None, node_size=400, labels_size=10, filename=None):
    """
    Visualize a protein-protein interaction network using matplotlib.
    
    Parameters:
    - G (networkx.Graph): The graph to visualize.
    - gene_list (list, optional): List of gene names for the title.
    - node_size (int, optional): Size of the nodes.
    - labels_size (int, optional): Font size for labels.
    - filename (str, optional): If provided, save the plot to this filename.
    """
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(15, 15))
    
    # Add a title
    if gene_list:
        title = f"Protein-Protein Interaction Network for {', '.join(gene_list)}"
    else:
        title = "Protein-Protein Interaction Network"
    plt.title(title)
    
    # Draw the network
    pos = nx.spring_layout(G, k=0.3, iterations=50)
    
    # Scale edge widths
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    min_width = 0.1
    max_width = 2.0
    epsilon = 1e-10
    edge_weights = [min_width + (w - min(edge_weights)) * (max_width - min_width) / (max(edge_weights) - min(edge_weights) + epsilon) for w in edge_weights]
    
    # Draw edges with scaled widths
    nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color='lightgray', ax=ax)
    
    # Draw nodes with a black border
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=node_size, ax=ax, linewidths=0.5, edgecolors='black',)
    
    # Draw labels shifted a bit to the upper right
    label_pos = {key: [value[0]+0.01, value[1]+0.01] for key, value in pos.items()}
    
    # Draw labels
    nx.draw_networkx_labels(G, label_pos, font_weight='bold', font_size=labels_size, font_color='black', ax=ax)
    
    if filename:
        plt.savefig(filename)
    
    plt.show()


    import requests

def get_gene_symbol_uniprot(uniprot_id):
    """
    Given a UniProt ID, fetches the corresponding gene symbol.

    Parameters:
    - uniprot_id (str): A UniProt ID for which the gene symbol is to be fetched.

    Returns:
    - str: The gene symbol corresponding to the UniProt ID.
    """
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.txt'
    response = requests.get(url)
    
    if response.status_code == 200:
        for line in response.text.split('\n'):
            if line.startswith('GN'):
                gene_info = line.split()
                for info in gene_info:
                    if info.startswith("Name="):
                        gene_name = info.split("=")[1]
                        return gene_name.rstrip(';')
    else:
        print(f"Failed to fetch data for UniProt ID: {uniprot_id}")
        return None


def get_uniprot_id(gene_symbols):
    """
    Given a list of gene symbols, fetches the corresponding UniProt IDs.

    Parameters:
    - gene_symbols (list): A list of gene symbols for which the UniProt IDs are to be fetched.

    Returns:
    - dict: A dictionary mapping gene symbols to their corresponding UniProt IDs.
    """
    uniprot_ids = {}
    
    for gene_symbol in gene_symbols:
        url = f'https://rest.uniprot.org/uniprotkb/search?size=1&query={gene_symbol}&fields=accession%2Cgene_names'
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            data = json.loads(response.text)
            if data['results']:
                uniprot_id = data['results'][0]['primaryAccession']
                uniprot_ids[gene_symbol] = uniprot_id
        else:
            print(f"Failed to fetch data for gene symbol: {gene_symbol}")
    
    return uniprot_ids