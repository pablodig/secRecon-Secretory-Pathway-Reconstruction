from Bio import Entrez
import time
import requests
import os
from Bio import SeqIO       
import shutil
import pandas as pd
import requests
from bs4 import BeautifulSoup
import json


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

    if int(record["Count"]) > 0:

        # Get the Gene ID of the first result
        gene_id = record["IdList"][0]

        # Fetch the gene information, maximum 3 tries
        fetch_tries = 0
        fetch_found = True
        while fetch_tries < 3 and fetch_found:
            gene_fetch_tries = 0
            gene_fetch_found = True
            gene_info = []
            while gene_fetch_tries < 3 and gene_fetch_found:
                try:
                    gene_info = Entrez.efetch(db="gene", id=gene_id, retmode="xml")    
                    gene_fetch_found = False
                except Exception as e:
                    time.sleep(3)
                    gene_fetch_tries += 1
                    print(f"Gene ID {gene_id}; An error occurred: {e}")
            if gene_fetch_tries > 2:
                print(f"{gene_EntrezID}. Gene ID {gene_id}; An error occurred")

            if gene_info == []:
                KeyError("sda")
            
            try:
                gene_record = Entrez.read(gene_info)
                fetch_found = False
            except Exception as e:
                time.sleep(3)
                fetch_tries += 1                
    return gene_record

def description_fetch_and_try(human_ortholog_EntrezID):
    human_fetch_tries = 0
    human_fetch_found = True
    gene_description = []

    while human_fetch_tries < 3 and human_fetch_found:
        try:
            search_results_human = Entrez.esearch(db="gene", term=human_ortholog_EntrezID)     
            record_human = Entrez.read(search_results_human)
            gene_id_human = record_human["IdList"][0]
            gene_info_human = Entrez.efetch(db="gene", id=gene_id_human, retmode="xml")
            gene_record_human = Entrez.read(gene_info_human)
            gene_description = gene_record_human[0]['Entrezgene_summary']
            human_fetch_found = False
        except Exception as e:
            time.sleep(3)
            human_fetch_tries += 1
            # print(f"Human_Description {human_ortholog_EntrezID}; An error occurred: {e}")     
    if human_fetch_tries > 2:
        print(f"Human_Description {human_ortholog_EntrezID}; An error occurred")


    return gene_description



def get_gene_ids(gene_entrezID, target_tax_id):
    url = f'https://www.ncbi.nlm.nih.gov/gene/{gene_entrezID}/ortholog/'
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        script_tag = soup.find('script', string=lambda text: 'var appData =' in text)
        script_content = script_tag.text if script_tag else ""
        gene_data_start = script_content.find('appData.genes =') + len('appData.genes =')
        gene_data_end = script_content.find(';', gene_data_start)
        gene_data_json = script_content[gene_data_start:gene_data_end]
        genes_data = json.loads(gene_data_json)
        for gene_info in genes_data:
            if gene_info.get('tax_id') == int(target_tax_id):
                gene_id_orth = (gene_info.get('gene_id'))
        return gene_id_orth