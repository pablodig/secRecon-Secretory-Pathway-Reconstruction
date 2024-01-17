from Bio import Entrez
import time
import requests
import os
from Bio import SeqIO       
import shutil
import pandas as pd

def read_gene_ID(gene_EntrezID):
    gene_fetch_tries = 0
    gene_fetch_found = True
    record = []
    while gene_fetch_tries < 3 and gene_fetch_found:
        try:
            search_results = Entrez.esearch(db="gene", term = gene_EntrezID)            
            record = Entrez.read(search_results)
            gene_fetch_found = False
        except Exception as e:
            time.sleep(3)
            gene_fetch_tries += 1
    if gene_fetch_tries > 2:        
        print(f"EntrezID of {gene_EntrezID} could not be read; An error occurred")

    return record