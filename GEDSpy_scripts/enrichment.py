#Requirements import
import urllib.request 
import re
import pandas as pd
from tqdm import tqdm
import json
import zipfile
import requests
import gzip
import shutil
import numpy as np
from scipy import stats
import os
from bs4 import BeautifulSoup
import warnings
import copy
import networkx as nx
from pyvis.network import Network
import webbrowser
import tkinter as tk
from collections import Counter
from datetime import datetime
import pkg_resources
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram
import math
import gdown
import sqlite3
import pandas as pd
from scipy.stats import fisher_exact, binom_test
from matplotlib.gridspec import GridSpec
from matplotlib import rc
rc("svg", fonttype='path')







   

pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore")

#Geta data directory



#Geta data directory

class PathMetadata:
    def __init__(self):
        
        def get_package_directory():
            return pkg_resources.resource_filename(__name__, '')
        
        self._cwd = get_package_directory()
        self.path_inside = os.path.join(self._cwd, 'data')
        self.path_in_inside = os.path.join(self._cwd, 'data', 'in_use')
        self.path_tmp = os.path.join(self._cwd, 'data', 'tmp')
        
        os.makedirs(self.path_inside, exist_ok=True)
        os.makedirs(self.path_in_inside, exist_ok=True)
        os.makedirs(self.path_tmp, exist_ok=True)

    def __repr__(self):
        return (f"PathMetadata(\n"
                f"  path_inside='{self.path_inside}',\n"
                f"  path_in_inside='{self.path_in_inside}',\n"
                f"  path_tmp='{self.path_tmp}'\n"
                f")")



class GetDataRaw(PathMetadata):
    
    
    def __init__(self):
       super().__init__()
    
    def get_raw_REACTOME(self):
        
        """
        This method gets the REACTOME data.
    
        Returns:
           dict: REACTOME data
        """
        
        try:
            with open(os.path.join(self.path_inside, 'reactome_jbio.json'), 'r') as json_file:
                reactome_jbio = json.load(json_file)  
        
            return reactome_jbio
        
        except:
            print("Something went wrong. Check the function input data and try again!")



    def get_raw_REF_GEN(self):
        
        """
        This method gets the REF_GEN which is the combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes for scientific use.
         
        Returns:
           dict: Combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes
        """
        
        try:
            
            with open(os.path.join(self.path_inside, 'gene_dictionary_jbio.json'), 'r') as json_file:
                gene_dictionary = (json.load(json_file))
        
            return gene_dictionary
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_raw_REF_GEN_RNA_SEQ(self):
        
        """
        This method gets the tissue-specific RNA-SEQ data including:
            -human_tissue_expression_HPA
            -human_tissue_expression_RNA_total_tissue
            -human_tissue_expression_fetal_development_circular

      
        Returns:
           dict: Tissue specific RNA-SEQ data
        """
        
        try:
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_HPA.json'), 'r') as json_file:
                human_tissue_expression_HPA = json.load(json_file)
                        
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_RNA_total_tissue.json'), 'r') as json_file:
                human_tissue_expression_RNA_total_tissue = json.load(json_file)
                
            
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_fetal_development_circular.json'), 'r') as json_file:
                human_tissue_expression_fetal_development_circular = json.load(json_file)
                
            
          

            rna_seq_list = {'tissue_expression_HPA':human_tissue_expression_HPA, 
                            'tissue_expression_RNA_total_tissue':human_tissue_expression_RNA_total_tissue, 
                            'tissue_expression_fetal_development_circular':human_tissue_expression_fetal_development_circular}
        
            return rna_seq_list
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            
            
            
    def get_raw_HPA(self):
        
        """
        This method gets the HPA (Human Protein Atlas) data.
         
        Returns:
           dict: HPA data
        """
        
        try:
            
            with open(os.path.join(self.path_inside, 'HPA_jbio.json'), 'r') as json_file:
                HPA_jbio = json.load(json_file)
                
                
        
            return HPA_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")




    def get_raw_DISEASES(self):
        
        """
        This method gets the DISEASES data.

        Returns:
           dict: DISEASES data
        """
        
        try:
            
            #load diseases
            with open(os.path.join(self.path_inside, 'diseases_jbio.json'), 'r') as json_file:
                disease_dict_jbio = json.load(json_file)
                
        
            return disease_dict_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_raw_ViMIC(self):
        
        """
        This method gets the ViMIC data.
         
        Returns:
           dict: ViMIC data
        """
        
        try:
            
            #load viral diseases
            with open(os.path.join(self.path_inside, 'viral_diseases_jbio.json'), 'r') as json_file:
                viral_dict_jbio = json.load(json_file)
            
        
            return viral_dict_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_raw_KEGG(self):
        
        """
        This method gets the KEGG data.

        Returns:
           dict: KEGG data
        """
        
        try:
            
            #load kegg
            with open(os.path.join(self.path_inside, 'kegg_jbio.json'), 'r') as json_file:
                kegg_jbio = json.load(json_file)
                 
        
            return kegg_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")




    def get_raw_GO(self):
        
        """
        This method gets the GO-TERM data.
         
        Returns:
           dict: GO-TERM data
        """
        
        try:
        
            with open(os.path.join(self.path_inside, 'goterm_jbio.json'), 'r') as json_file:
                go_term_jbio = json.load(json_file)
           
        
            return go_term_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_raw_IntAct(self):
        
        """
        This method gets the IntAct data.

        Returns:
           dict: IntAct data
        """
        
        try:
            with open(os.path.join(self.path_inside, 'IntAct_jbio.json'), 'r') as json_file:
                IntAct_dict = json.load(json_file)
                
        
            return IntAct_dict
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")


    def get_raw_STRING(self):
        
        """
        This method gets the STRING data.
         
        Returns:
           dict: STRING data
        """
        
        try:
            with open(os.path.join(self.path_inside, 'string_jbio.json'), 'r') as json_file:
                string_dict = json.load(json_file)
                
                
            return string_dict
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            
    def get_raw_CellTalk(self):
        
        """
        This method gets the CellTalk data.
         
        Returns:
           dict: CellTalk data
        """
        
        try:
            with open(os.path.join(self.path_inside, 'cell_talk_jbio.json'), 'r') as json_file:
                cell_talk = json.load(json_file)
                
                
            return cell_talk
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
    
    def get_raw_CellPhone(self):
        
        """
        This method gets the CellPhone data.
         
        Returns:
           dict: CellPhone data
        """
        
        try:
            with open(os.path.join(self.path_inside, 'cell_phone_jbio.json'), 'r') as json_file:
                cell_phone = json.load(json_file)
                
                
            return cell_phone
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            




class GetData(PathMetadata):
    
    def __init__(self):
       super().__init__()
    
    def get_REACTOME(self):
        
        """
        This method gets the REACTOME data including the id to connect with REF_GENE by id_reactome
    
        Returns:
           data_frame: REACTOME data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'reactome_jbio_dict.json'), 'r') as json_file:
                reactome_jbio = json.load(json_file)
        
            return reactome_jbio
        
        except:
            print("Something went wrong. Check the function input data and try again!")



    def get_REF_GEN(self):
        
        """
        This method gets the REF_GEN which is the combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes for scientific use.
         
        Returns:
           dict: Combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes
        """
        
        try:
            
            with open(os.path.join(self.path_in_inside, 'gene_dictionary_jbio_annotated.json'), 'r') as json_file:
                gene_dictionary = json.load(json_file)
        
            return gene_dictionary
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_REF_GEN_RNA_SEQ(self):
        
        """
        This method gets the tissue-specific RNA-SEQ data including:
            -human_tissue_expression_HPA
            -human_tissue_expression_RNA_total_tissue
            -human_tissue_expression_fetal_development_circular

      
        Returns:
           dict: Tissue specific RNA-SEQ data
        """
        
        try:
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_HPA.json'), 'r') as json_file:
                human_tissue_expression_HPA = json.load(json_file)
                        
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_RNA_total_tissue.json'), 'r') as json_file:
                human_tissue_expression_RNA_total_tissue = json.load(json_file)
                
            
            
            with open(os.path.join(self.path_in_inside, 'tissue_expression_fetal_development_circular.json'), 'r') as json_file:
                human_tissue_expression_fetal_development_circular = json.load(json_file)
                
                

            rna_seq_list = {'tissue_expression_HPA':human_tissue_expression_HPA, 
                            'tissue_expression_RNA_total_tissue':human_tissue_expression_RNA_total_tissue, 
                            'tissue_expression_fetal_development_circular':human_tissue_expression_fetal_development_circular}
        
            return rna_seq_list
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            
            
            
    def get_HPA(self):
        
        """
        This method gets the HPA (Human Protein Atlas) data including the id to connect with REF_GENE by id_HPA
         
        Returns:
           dict: HPA data
        """
        
        try:
            
            with open(os.path.join(self.path_in_inside, 'HPA_jbio_dict.json'), 'r') as json_file:
                HPA_jbio = json.load(json_file)
                
                
        
            return HPA_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")




    def get_DISEASES(self):
        
        """
        This method gets the DISEASES data including the id to connect with REF_GENE by id_diseases

        Returns:
           dict: DISEASES data
        """
        
        try:
            
            #load diseases
            with open(os.path.join(self.path_in_inside, 'disease_jbio_dict.json'), 'r') as json_file:
                disease_dict_jbio = json.load(json_file)
                
        
            return disease_dict_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_ViMIC(self):
        
        """
        This method gets the ViMIC data including the id to connect with REF_GENE by id_viral_diseases
         
        Returns:
           dict: ViMIC data
        """
        
        try:
            
            #load viral diseases
            with open(os.path.join(self.path_in_inside, 'viral_jbio_dict.json'), 'r') as json_file:
                viral_dict_jbio = json.load(json_file)
            
        
            return viral_dict_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_KEGG(self):
        
        """
        This method gets the KEGG data including the id to connect with REF_GENE by id_kegg

        Returns:
           dict: KEGG data
        """
        
        try:
            
            #load kegg
            with open(os.path.join(self.path_in_inside, 'kegg_jbio_dict.json'), 'r') as json_file:
                kegg_jbio = json.load(json_file)
                 
        
            return kegg_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")




    def get_GO(self):
        
        """
        This method gets the GO-TERM data including the id to connect with REF_GENE by id_go
         
        Returns:
           data_frame: GO-TERM data
        """
        
        try:
        
            with open(os.path.join(self.path_in_inside, 'goterm_jbio_dict.json'), 'r') as json_file:
                go_term_jbio = (json.load(json_file))
           
        
            return go_term_jbio
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")



    def get_IntAct(self):
        
        """
        This method gets the IntAct data including the id to connect with REF_GENE by id_IntAct

        Returns:
           data_frame: IntAct data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'intact_jbio_dict.json'), 'r') as json_file:
                IntAct_dict = json.load(json_file)
                
        
            return IntAct_dict
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")


    def get_STRING(self):
        
        """
        This method gets the STRING data including the id to connect with REF_GENE by id_string
         
        Returns:
           data_frame: STRING data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'string_jbio_dict.json'), 'r') as json_file:
                string_dict = json.load(json_file)
                
                
            return string_dict
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            
            
    def get_CellTalk(self):
        
        """
        This method gets the CellTalk data after adjustment.
         
        Returns:
           dict: CellTalk data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'cell_talk_jbio.json'), 'r') as json_file:
                cell_talk = json.load(json_file)
                
                
            return cell_talk
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
    
    def get_CellPhone(self):
        
        """
        This method gets the CellPhone data after adjustment.
         
        Returns:
           dict: CellPhone data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'cell_phone_jbio.json'), 'r') as json_file:
                cell_phone = json.load(json_file)
                
                
            return cell_phone
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            
     
    def get_interactions(self):
        
        """
        This method gets the CellPhone & CellTalk data including the id to connect with REF_GENE by id_go.
         
        Returns:
           dict: CellPhone data
        """
        
        try:
            with open(os.path.join(self.path_in_inside, 'cell_int_jbio.json'), 'r') as json_file:
                cell_int = json.load(json_file)
                
                
            return cell_int
        
        except:
            print('\n')
            print("Something went wrong. Check the function input data and try again!")
            
            

            
           
 
    
    

class Enrichment(GetData):
    
    def __init__(self):
        self.features_list = None
        self.genome = None
        self.found_genes = None
        self.species_study = ['Mus musculus', 'Homo sapiens', 'Rattus norvegicus']
        self.species_ids = None
        self.species_genes = ['Mus musculus', 'Homo sapiens', 'Rattus norvegicus']
        self.ids = None
        self.HPA = None
        self.STRING = None
        self.GO_TERM = None
        self.REACTOME = None
        self.CellCon = None
        self.IntAct = None
        self.ViMIC = None
        self.Diseases = None
        self.KEGG = None
        self.mapper = None
        self.RNA_SEQ = None
        
        super().__init__()
        
        
        
    @property
    def show_founded_features(self):
        return self.found_genes['found_genes']
        
        
    @property
    def show_non_founded_features(self):
        return self.found_genes['not_found']



    @property
    def get_gene_info(self):
        
        if isinstance(self.genome, pd.DataFrame):
            
            return self.genome.to_dict(orient = 'list')
        
        else:
            raise ValueError('\nLack of Genome data...')



    @property
    def get_RNA_SEQ(self):
        
        if self.RNA_SEQ != None:
            
            return self.RNA_SEQ
        
        else:
            raise ValueError('\nLack of enrichment of RNA-SEQ data...')



    @property
    def get_HPA(self):
        
        if self.HPA != None:
            
            return self.HPA
        
        else:
            raise ValueError('\nLack of enrichment of HPA data...')



    @property
    def get_STRING(self):
        
        if self.STRING != None:
            
            return self.STRING
        
        else:
            raise ValueError('\nLack of enrichment of STRING data...')



    @property
    def get_GO_TERM(self):
        
        if self.GO_TERM != None:
            
            return self.GO_TERM
        
        else:
            raise ValueError('\nLack of enrichment of GO-TERM data...')



    @property
    def get_REACTOME(self):
        
        if self.REACTOME != None:
            
            return self.REACTOME
        
        else:
            raise ValueError('\nLack of enrichment of REACTOME data...')

     
  
    @property
    def get_CellCon(self):
        
        if self.CellCon != None:
            
            return self.CellCon
        
        else:
            raise ValueError('\nLack of enrichment of CellCon data...')

     
    @property
    def get_IntAct(self):
        
        if self.IntAct != None:
            
            return self.IntAct
        
        else:
            raise ValueError('\nLack of enrichment of IntAct data...')

     
        
    @property
    def get_KEGG(self):
        
        if self.KEGG != None:
            
            return self.KEGG
           
        else:
            raise ValueError('\nLack of enrichment of KEGG data...')

        
    @property
    def get_DISEASES(self):
        
        if self.Diseases != None:
            
            return self.Diseases
        
        else:
            raise ValueError('\nLack of enrichment of Diseases data...')

        
        
    @property
    def get_ViMIC(self):
        
        if self.ViMIC != None:
            
            return self.ViMIC
        
        else:
            raise ValueError('\nLack of enrichment of ViMIC data...')

            
            
    @property
    def get_columns_names(self):
        conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
        
        cursor = conn.cursor()

        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        
        conn.close()
        
        for table in tables:
            print(table[0])
          
            
    @property
    def get_results(self):
        
        results = {}

        try:
            results['gene_info'] = self.get_gene_info
        except:
            pass
        
        
        try:
            results['HPA'] = self.get_HPA
        except:
            pass
              
        try:    
            results['STRING'] = self.get_STRING
        except:
            pass
           
        
        try:    
            results['GO-TERM'] = self.get_GO_TERM
        except:
            pass
              
        
        try:
            results['REACTOME']  = self.get_REACTOME
        except:
            pass
              
            
        try:
            results['CellConnections'] = self.get_CellCon
        except:
            pass
              
             
        try:
            results['IntAct'] = self.get_IntAct
        except:
            pass
            
        
        try:
            results['KEGG'] = self.get_KEGG
        except:
            pass
        
        
        try:
            results['DISEASES'] = self.get_DISEASES
        except: 
            pass
              
            
        try:
            results['ViMIC'] = self.get_ViMIC
        except:
            pass
        
        
        try:
            results['RNA-SEQ'] = self.get_RNA_SEQ   
        except:
            pass
        
        
        results['species'] = {}
        results['species']['species_genes'] = self.species_genes   
        results['species']['species_study'] = self.species_study  
            
        return results
            
    
    
    def deserialize_data(self, value):
        try:
            return json.loads(value)
        except (json.JSONDecodeError, TypeError):
            return value
                
    
    def load_genome(self):
        
        print('\n')
        print('Metadata loading...')
    
        conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
        
        query = "SELECT * FROM RefGenome;"
        
        df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
        
        conn.close()

        self.genome = df
    
    
            
    def spec_dic(self):
        
        species = self.species_genes 
        
        if isinstance(species, list) and len(species) == 1 and 'Homo sapiens' in species:
            tf_h = ['Homo sapiens' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 1 and 'Mus musculus' in species:
            tf_h = ['Mus musculus' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 1 and 'Rattus norvegicus' in species:
            tf_h = ['Rattus norvegicus' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 2 and 'Homo sapiens' in species and 'Mus musculus' in species:
            tf_h = ['Homo sapiens' in x or 'Mus musculus' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 2 and 'Homo sapiens' in species and 'Rattus norvegicus' in species:
            tf_h = ['Homo sapiens' in x or 'Rattus norvegicus' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 2 and 'Mus musculus' in species and 'Rattus norvegicus' in species:
            tf_h = ['Rattus norvegicus' in x or 'Mus musculus' in x for x in self.genome['species']]
        elif isinstance(species, list) and len(species) == 3 and 'Mus musculus' in species and 'Rattus norvegicus' in species and 'Homo sapiens' in species:
            tf_h = [True for x in self.genome['species']]  # Selects all species
        else:
            raise ValueError("Invalid species specified.")
            
        ids = [value for value, flag in zip(self.genome['sid'], tf_h) if flag]
        
        self.species_ids = ids
    

    
    def find_fetures_id(self):
        
        genome_dict = pd.DataFrame(self.genome)
        genome_dict = genome_dict.reset_index(drop = True)
        gene = []
        not_found = []
        ids = []
        for fet in tqdm(self.features_list):
            idg = genome_dict['sid'][genome_dict['possible_names'].apply(lambda x: fet.upper() in x)]
            if len(idg) > 0:
                ids.append(int(idg))
                gene.append(fet)
            else:
                not_found.append(fet)
                
        self.found_genes = {'found_genes':gene, 'found_ids':ids, 'not_found':not_found}
    
    
    def v1_is_in_v2(self):
        v1 = self.found_genes['found_ids']
        v2 = self.species_ids
        v3 = [x for x in v1 if x in v2]
        self.ids =  v3
        
    
    
    def return_dictionary(self):
        genome_dict = pd.DataFrame(self.genome)
        genome_dict = genome_dict[genome_dict['sid'].isin(self.ids)]
        genome_dict = genome_dict.reset_index(drop = True)
        
        self.genome = genome_dict

    
    
    def add_found_names(self):
     
        mapper = dict(zip(self.found_genes['found_ids'], self.found_genes['found_genes']))
        gene_dictionary = pd.DataFrame(self.genome)
        gene_dictionary['found_names'] = gene_dictionary['sid']
        gene_dictionary['found_names'] = gene_dictionary['found_names'].map(mapper)
        
        self.genome = gene_dictionary
        
        
    def enriche_RNA_SEQ(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:


            genes_hs = [str(x) for x in self.genome['gene_Homo_sapiens'] if x == x] + ['tissue']
    
            df = self.get_REF_GEN_RNA_SEQ()
            
            for k in df.keys():
                tmp = df[k]
                tmp = {key: tmp[key] for key in tmp if key in genes_hs}
                if len(tmp.keys()) == 1:
                    tmp = {}
                
                df[k] = tmp

            
            self.RNA_SEQ = df
            

        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
        
    def enriche_specificiti(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            hpa_datasets = [
                "HPA_RNA_tissue",
                "HPA_RNA_single_cell",
                "HPA_RNA_cancer",
                "HPA_RNA_brain",
                "HPA_RNA_blood",
                "HPA_RNA_blood_lineage",
                "HPA_RNA_cell_line",
                "HPA_RNA_mouse_brain_region",
                "HPA_subcellular_location",
                "HPA_blood_markers"
            ]
            
            
            ids = [int(x) for x in self.genome['id_HPA'] if x == x]
            
            full_dict = {}
            for k in hpa_datasets:
            
                query = f"SELECT * FROM {k} WHERE id IN ({', '.join(map(str, ids))});"
    
                df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
                
                
                df = pd.merge(df, self.genome[['id_HPA', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_HPA')
                df = df.drop(['id_HPA'], axis = 1)
                                        
    
                full_dict[k] = df.to_dict(orient = 'list')
                
                
            conn.close()
    
            self.HPA = full_dict
        
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
            
        
    def enriche_KEGG(self):

        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_KEGG'] if x == x]
            
            query = f"SELECT * FROM KEGG WHERE id IN ({', '.join(map(str, ids))});"
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_KEGG', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_KEGG')
            df = df.drop(['id_KEGG'], axis = 1)
                                    
            conn.close()
    
            self.KEGG = df.to_dict(orient = 'list')
            
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
             
            
        
            
    def enriche_GOTERM(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:
        
            final_dict = {}
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_GO'] if x == x]
            
            
             
            species = self.species_study 
            species = ', '.join(map(lambda x: f"'{x}'", species))
            
            
            query = f"""SELECT * FROM GO_gene_info 
            WHERE id IN ({', '.join(map(str, ids))})
              AND species IN ({species});
            """
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_GO', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_GO')
            df = df.drop(['id_GO'], axis = 1)
                                    
            final_dict['gene_info'] = df.to_dict(orient = 'list')
                
            ids = ', '.join(map(lambda x: f"'{x}'", df['GO_id']))
            hids = list(df['GO_id'])

            query = f"SELECT * FROM GO_go_names WHERE GO_id IN ({ids});"
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            final_dict['go_names'] = df.to_dict(orient = 'list')
            
            query = f"SELECT * FROM GO_hierarchy WHERE GO_id IN ({ids});"
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            for cl in ['is_a_ids', 'part_of_ids', 'has_part_ids', 'regulates_ids', 'negatively_regulates_ids', 'positively_regulates_ids']:
                tmp_val = [x if x in hids else None for x in df[cl]]
                df[cl] = tmp_val
    
            final_dict['hierarchy'] = df.to_dict(orient = 'list')
            
            del df
            
            conn.close()
    
    
            self.GO_TERM = final_dict
            
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
     
        
           
           
            
        
    def enriche_REACTOME(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_reactome'] if x == x]
            
            query = f"SELECT * FROM REACTOME WHERE id IN ({', '.join(map(str, ids))});"
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_reactome', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_reactome')
            df = df.drop(['id_reactome'], axis = 1)
            
            conn.close()
            
            self.REACTOME = df.to_dict(orient = 'list')
       
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
    
    
    def enriche_DISEASES(self):
                
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_diseases'] if x == x]
            
            query = f"SELECT * FROM disease WHERE id IN ({', '.join(map(str, ids))});"
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_diseases', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_diseases')
            df = df.drop(['id_diseases'], axis = 1)
            
    
            conn.close()
    
            self.Diseases = df.to_dict(orient = 'list')
            
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
        
    def enriche_ViMIC(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_viral_diseases'] if x == x]
            
            query = f"SELECT * FROM ViMIC WHERE id IN ({', '.join(map(str, ids))});"
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_viral_diseases', 'found_names']], how = 'left', left_on = 'id', right_on = 'id_viral_diseases')
            df = df.drop(['id_viral_diseases'], axis = 1)
            
            conn.close()
    
            self.ViMIC = df.to_dict(orient = 'list')
        
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
        
    def enriche_IntAct(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            final_dict = {}
    
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_IntAct'] if x == x]
            
             
            species = self.species_study 
            species = ', '.join(map(lambda x: f"'{x}'", species))
            
            
            
            query = f"""
            SELECT * 
            FROM IntAct_gene_product 
            WHERE id_1 IN ({', '.join(map(str, ids))}) 
              AND id_2 IN ({', '.join(map(str, ids))})
              AND species IN ({species});
            """
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_IntAct', 'found_names']], how = 'left', left_on = 'id_1', right_on = 'id_IntAct')
            df = df.drop(['id_IntAct'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_1'})
    
            
            df = pd.merge(df, self.genome[['id_IntAct', 'found_names']], how = 'left', left_on = 'id_2', right_on = 'id_IntAct')
            df = df.drop(['id_IntAct'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_2'})
            
            final_dict['gene_products'] = df.to_dict(orient = 'list')
            
            
            query = f"""
            SELECT * 
            FROM IntAct_non_gene_product 
            WHERE id_1 IN ({', '.join(map(str, ids))}) 
              OR id_2 IN ({', '.join(map(str, ids))})
              AND species IN ({species});
            """
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            final_dict['non_gene_products'] = df.to_dict(orient = 'list')
    
            conn.close()
    
            self.IntAct = final_dict
            
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
 
    def enriche_STRING(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:
        
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_STRING'] if x == x]
            
            
            
            species = self.species_study 
            species = ', '.join(map(lambda x: f"'{x}'", species))

            
            query = f"""
            SELECT * 
            FROM STRING 
            WHERE protein1 IN ({', '.join(map(str, ids))}) 
              AND protein2 IN ({', '.join(map(str, ids))})
              AND species IN ({species});
            """      
    
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_STRING', 'found_names']], how = 'left', left_on = 'protein1', right_on = 'id_STRING')
            df = df.drop(['id_STRING'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_1'})
    
            
            df = pd.merge(df, self.genome[['id_STRING', 'found_names']], how = 'left', left_on = 'protein2', right_on = 'id_STRING')
            df = df.drop(['id_STRING'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_2'})
    
            conn.close()
    
            self.STRING = df.to_dict(orient = 'list')
          
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
            
        
    def enriche_CellCon(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:
            
            
            tmp = {}
            ###################################################################
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_cell_int'] if x == x]
            
            
             
            species = self.species_study 
            species = ', '.join(map(lambda x: f"'{x}'", species))
            
            
            query = f"""
            SELECT * 
            FROM CellInteractions 
            WHERE protein_id_2 IN ({', '.join(map(str, ids))})
              AND Species IN ({species});
            """  
     
                
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_1', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_1'})
    
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_2', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_2'})
    
            tmp['interactor1'] = df.to_dict(orient = 'list')
            ###################################################################
            
            
            query = f"""
            SELECT * 
            FROM CellInteractions 
            WHERE protein_id_1 IN ({', '.join(map(str, ids))}) 
              AND Species IN ({species});
            """  
     
                
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_1', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_1'})
    
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_2', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_2'})
            
            tmp['interactor2'] = df.to_dict(orient = 'list')
            
            
            ###################################################################
            
            
            query = f"""
            SELECT * 
            FROM CellInteractions 
            WHERE protein_id_1 IN ({', '.join(map(str, ids))}) 
              AND protein_id_2 IN ({', '.join(map(str, ids))})
              AND Species IN ({species});
            """  
     
                
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_1', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_1'})
    
            
            df = pd.merge(df, self.genome[['id_cell_int', 'found_names']], how = 'left', left_on = 'protein_id_2', right_on = 'id_cell_int')
            df = df.drop(['id_cell_int'], axis = 1)
            df = df.rename(columns={'found_names': 'found_names_2'})
    
            tmp['mutual'] = df.to_dict(orient = 'list')

            
            conn.close()
    
    
            self.CellCon = tmp
            
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
              
            
        
        
    def select_features(self, features_list:list):
        

            self.features_list = features_list
            
            self.load_genome()
            self.spec_dic()
            
            print('\nFeatures selection...')

            self.find_fetures_id()
            self.v1_is_in_v2()
            self.return_dictionary()
            self.add_found_names()
            
            nf = self.show_non_founded_features
            
            if len(nf) > 0:
                print('\nSome features were not found in Database:')
                for n in nf:
                    print(f' --> {n}')
            
    
              
              
        
    def full_enrichment(self):
        
        if isinstance(self.genome, pd.DataFrame) and 'found_names' in self.genome.columns:

            print('\nSpecificity enrichment...')

            self.enriche_specificiti()
            
            print('\nEnrichment with KEGG information...')

            self.enriche_KEGG()
            
            print('\nEnrichment with GO-TERM information...')

            self.enriche_GOTERM()
            
            print('\nEnrichment with REACTOME information...')

            self.enriche_REACTOME()
            
            print('\nEnrichment with DISEASES information...')

            self.enriche_DISEASES()
            
            print('\nEnrichment with ViMIC information...')

            self.enriche_ViMIC()
            
            print('\nEnrichment with IntAct information...') 

            self.enriche_IntAct()
            
            print('\nEnrichment with STRING information...')

            self.enriche_STRING()
            
            print('\nEnrichment with CellConnections information...')

            self.enriche_CellCon()
            
            print('\nEnrichment with tissue specific RNA-SEQ information...')

            self.enriche_RNA_SEQ()
            
        
        else:
            raise ValueError('\nSelect features to enriche first! Use select_features() method...')
       
        
        
        
        
class Analysis(Enrichment):
    
    def __init__(self, input_data:dict):
        
        self.input_data = input_data
        self.KEGG_stat = None
        self.GO_stat = None
        self.REACTOME_stat = None
        self.DISEASE_stat = None
        self.ViMIC_stat = None
        self.features_interactions = None
        self.specificity_stat = None
        self.KEGG_net = None
        self.GO_net = None
        self.REACTOME_net = None
        # self.ViMIC_net = None
        # self.DISEASES_net = None
        self.network_stat = {'test':'BIN', 'adj':'BF', 'p_val':0.05}
        self.go_grade = 1
        self.occ = None
        self.interaction_strength = 900
        self.interaction_source = ['Affinomics',
                                'Alzheimers',
                                'BioCreative',
                                'Cancer',
                                'Cardiac',
                                'Chromatin',
                                'Coronavirus',
                                'Diabetes',
                                "Huntington's",
                                'IBD',
                                'Neurodegeneration',
                                'Parkinsons',
                                'STRING']
        
        
        
        super().__init__()

      

    @property
    def interactions_metadata(self):
        
        """
        This method returns current interactions parameters:
            
            -interaction_strength : int - value of enrichment strenght for STRING data:
                *900 - very high probabylity of interaction; 
                *700 - medium probabylity of interaction, 
                *400 - low probabylity of interaction, 
                *<400 - very low probabylity of interaction
                
                
            -interaction_source : list - list of sources for interaction estimation:
                *STRING: ['STRING']
                *IntAct: ['Affinomics', 'Alzheimers','BioCreative', 'Cancer',
                          'Cardiac', 'Chromatin', 'Coronavirus', 'Diabetes',
                          "Huntington's", 'IBD', 'Neurodegeneration', 'Parkinsons']
                  
        
        Returns:
            dict : {'interaction_strength' : int,
                    'interaction_source' : list} 
           
        """
        
        return {'interaction_strength':self.interaction_strength,
                'interaction_source':self.interaction_source} 
    
    
    
    def set_interaction_strength(self, value:int):
        
        """
        This method sets self.interaction_strength parameter.
        
        The 'interaction_strength' value is used for enrichment strenght of STRING data:
             *900 - very high probabylity of interaction; 
             *700 - medium probabylity of interaction, 
             *400 - low probabylity of interaction, 
             *<400 - very low probabylity of interaction
        
        """
        
        if value < 1000 and value > 0:
            self.interaction_strength = value
            
        else:
            
            raise ValueError('\nValue should be integer between 0 and 1000')
    
    def set_go_grade(self, grade:int):
        
        """
        This method sets self.go_grade parameter.
        
        The 'go_grade' value is used for GO-TERM data gradation [1-4].
             
        
        """
        
        if grade < 5 and grade > 0:
            self.go_grade = grade
            
        else:
            
            raise ValueError('\nValue should be integer between 0 and 1000')
    
    
    def set_interaction_source(self, sources_list:list):
        
        """
        This method sets self.interaction_source parameter.
        
        The 'interaction_source' value is list of sources for interaction estimation:
            
            *STRING / IntAct: ['STRING', 'Affinomics', 'Alzheimers','BioCreative', 
                               'Cancer', 'Cardiac', 'Chromatin', 'Coronavirus', 
                               'Diabetes', "Huntington's", 'IBD', 'Neurodegeneration', 
                               'Parkinsons']
        
        """
        
        values = ['STRING', 'Affinomics', 'Alzheimers','BioCreative', 
                           'Cancer', 'Cardiac', 'Chromatin', 'Coronavirus', 
                           'Diabetes', "Huntington's", 'IBD', 'Neurodegeneration', 
                           'Parkinsons']
        
        if all(value in values for value in sources_list):
            self.interaction_source = sources_list
            
        else:
            
            raise ValueError('\nValues in list should be included in STRING, Affinomics, Alzheimers, BioCreative, Cancer, Cardiac, Chromatin, Coronavirus, Diabetes, Huntington`s, IBD, Neurodegeneration, Parkinsons')
         
         
         
    @property
    def networks_metadata(self):
        
        """
        This method returns current networks creation parameters:
            
            -test : str - test type for enrichment overrepresentation analysis. 
                Available test:                    
                    *BIN - binomial test
                    *FISH - Fisher's exact test
                    
            -adj : str | None - p_value correction.
                Available correction:
                    *BF - Bonferroni correction
                    *BH - Benjamini-Hochberg correction
                    *None - lack of correction
                    
            -p_val : float - threshold for p-value in network creation
        
        
        Returns:
            dict : {'test': str, 
                    'adj': str | None, 
                    'p_val': float}
           
        """
        
        return self.network_stat
    
    
    def set_p_value(self, value:float):
        
        if value > 0 and value < 1:
            self.network_stat['p_val'] = value
            
        else:
            
            raise ValueError('\nValue should be float between 0 and 1')
            
            
    def set_test(self, test):
        
        if test.upper() == 'FISH' or test.upper() == 'BIN':
            self.network_stat['test'] = test.upper()
            
        else:
            
            raise ValueError('\nTest should be included in BIN or FISH')
            
            
            
    def set_correction(self, correction):
        
        if correction == None:
            self.network_stat['adj'] = None
            
        else:
            
            if correction.upper() == 'BF' or correction.upper() == 'BH':
                self.network_stat['adj'] = correction.upper()
                
            else:
                
                raise ValueError('\nTest should be included in BF or BH')
               
            
   
    
    @property
    def get_KEGG_statistics(self):
        
        if self.KEGG_stat != None:
            
            return self.KEGG_stat
        
        else:
            raise ValueError('\nNo data to return...')
            
            
    @property
    def get_REACTOME_statistics(self):
        
        if self.REACTOME_stat != None:
            
            return self.REACTOME_stat
        
        else:
            raise ValueError('\nNo data to return...')
        
        
    @property
    def get_GO_statistics(self):
        
        if self.GO_stat != None:
            
            return self.GO_stat
        
        else:
            raise ValueError('\nNo data to return...')
    
    
    @property
    def get_DISEASE_statistics(self):
        
        if self.DISEASE_stat != None:
            
            return self.DISEASE_stat
        
        else:
            raise ValueError('\nNo data to return...')
    
    
    @property
    def get_ViMIC_statistics(self):
        
        if self.ViMIC_stat != None:
            
            return self.ViMIC_stat
        
        else:
            raise ValueError('\nNo data to return...')
    
    
    @property
    def get_features_interactions_statistics(self):
        
        if self.features_interactions != None:
            
            return self.features_interactions
        
        else:
            raise ValueError('\nNo data to return...')
    
    
    @property
    def get_specificity_statistics(self):
        
        if self.specificity_stat != None:
            
            return self.specificity_stat
        
        else:
            raise ValueError('\nNo data to return...')
    
    
    @property
    def get_KEGG_network(self):
        
        if self.KEGG_net != None:
            
            return self.KEGG_net
        
        else:
            raise ValueError('\nNo data to return...')
            
        
    @property
    def get_REACTOME_network(self):
        
        if self.REACTOME_net != None:
            
            return self.REACTOME_net
        
        else:
            raise ValueError('\nNo data to return...')
            
         
    @property
    def get_GO_network(self):
        
        if self.GO_net != None:
            
            return self.GO_net
        
        else:
            raise ValueError('\nNo data to return...')
            
    
          
    # @property
    # def get_DISEASES_network(self):
        
    #     if self.DISEASES_net != None:
            
    #         return self.DISEASES_net
        
    #     else:
    #         raise ValueError('\nNo data to return...')
       
            
    # @property
    # def get_ViMIC_network(self):
        
    #     if self.ViMIC_net != None:
            
    #         return self.ViMIC_net
        
    #     else:
    #         raise ValueError('\nNo data to return...')
            
    
    
    def map_interactions_flat(self, row, mapping):
        interactions = [typ for col in row if col in mapping for typ in mapping[col]]
        return list(set(interactions))  
    
    
    def select_test(self, test, adj):
        try:
            test_string = ''
    
            if adj != None and adj.upper() in ['BF','BH']:
                test_string = test_string + 'adj_pval_'
            else:
                test_string = test_string + 'pval_'
            
            
            if test != None and test.upper() == 'BIN':
                test_string = test_string + 'BIN'
            elif test != None and test.upper() == 'FISH':
                test_string = test_string + 'FISH'
            else:
                test_string = test_string + 'BIN'
                
            
            if adj != None and adj.upper() == 'BF':
                test_string = test_string + '-BF'
            elif adj != None and adj.upper() == 'BH':
                test_string = test_string + '-BH'
            else:
                test_string = test_string + ''
            
            return test_string
        except:
            print('\n')
            print('Provided wrong test input!')
            
        
    

    def run_enrichment_tests(self, N, K, n, k):
        """
        Oblicza test Fishera i test dwumianowy dla wzbogacenia genów.
    
        Parametry:
        N (int): Liczba wszystkich genów w tle.
        K (int): Liczba genów związanych z danym terminem.
        n (int): Liczba genów w analizowanym zbiorze.
        k (int): Liczba genów w analizowanym zbiorze przypisanych do termu.
    
        Zwraca:
        dict: Wyniki testów (p-value i odds ratio dla Fishera oraz p-value dla dwumianowego).
        """
        fisher_table = [[k, K - k],
                        [n - k, N - K - (n - k)]]
        
        fisher_odds_ratio, fisher_p_value = fisher_exact(fisher_table, alternative='greater')
        
        p_background = K / N
        
        binomial_p_value = binom_test(k, n, p_background, alternative='greater')
        
        return {
            "fisher_p_value": fisher_p_value,
            "fisher_odds_ratio": fisher_odds_ratio,
            "binomial_p_value": binomial_p_value
        }



    def create_full_conections(self, go_data, grade = 1):
        
        go_data = go_data[['GO_id', 'is_a_ids', 'part_of_ids',
                'has_part_ids', 'regulates_ids', 'negatively_regulates_ids',
                'positively_regulates_ids']]

        
        for i in ['is_a_ids', 'part_of_ids', 'has_part_ids', 'regulates_ids', 
                          'negatively_regulates_ids', 'positively_regulates_ids']:

            go_data = go_data.explode(i)





        go_wide = pd.DataFrame()

        init_list = list(set(go_data['GO_id'][go_data['is_a_ids'].isin([None])]))
        full_list = list(set(go_data['GO_id'][~go_data['is_a_ids'].isin([None])]))



        for i in ['is_a_ids', 'part_of_ids', 'has_part_ids', 'regulates_ids', 
                          'negatively_regulates_ids', 'positively_regulates_ids']:
        
            go_tmp = go_data[['GO_id', i]][go_data[i].isin(init_list)]
            go_tmp.columns = ['grade_1', 'parent']
            
            go_wide = pd.concat([go_wide, go_tmp])

        final_wide = go_wide.copy()
        final_wide = final_wide.drop('parent', axis = 1)

        grade_set = set(go_wide['grade_1'])
        
        grade_dict = {}
        
        n = 1
        while(True):
            
            # print(f'Round {n}')
            

            if len(full_list) == 0:
                break
            
            go_wide_tmp = pd.DataFrame()
            

            full_list = [x for x in full_list if x not in grade_set]
            
            
            for i in ['is_a_ids', 'part_of_ids', 'has_part_ids', 'regulates_ids', 
                              'negatively_regulates_ids', 'positively_regulates_ids']:
            
                go_tmp = go_data[['GO_id', i]][go_data[i].isin(grade_set)]
                go_tmp.columns = [f'grade_{n+1}', 'parent2']
                
                go_wide_tmp = pd.concat([go_wide_tmp, go_tmp])
            
            if len(go_wide_tmp.index) == 0:
                break
            
            final_wide = pd.merge(final_wide, go_wide_tmp, left_on = f'grade_{n}', right_on = 'parent2', how = 'left') 
            final_wide = final_wide.drop('parent2', axis = 1)
            
            grade_dict[f'grade_1_grade_{n+1}'] = final_wide.groupby('grade_1').agg({f'grade_{n+1}': set}).reset_index().copy()
            
            if n != 1:
                final_wide = final_wide.drop(f'grade_{n}', axis = 1).drop_duplicates()
            
     
            n += 1

            grade_set = set(go_wide_tmp[f'grade_{n}'])
            
        
        primary = None
        for n, i in enumerate(grade_dict.keys()):
            
            if n+1 == grade:
                primary = grade_dict[i].reset_index(drop = True)
            else:
                if isinstance(primary, pd.DataFrame):
                    secondary = grade_dict[i].reset_index(drop = True)

                    for inx in primary.index:
                        go_name = primary.iloc[inx,0]
                        set1_cleaned = secondary[secondary.iloc[:,0] == go_name].reset_index(drop = True).iloc[0, 1]
                        set1_cleaned = {x for x in set1_cleaned if x == x}
                        
                        set2_cleaned = primary[primary.iloc[:,0] == go_name].reset_index(drop = True).iloc[0, 1]
                        set2_cleaned = {x for x in set2_cleaned if x == x}

                        set1_cleaned = set1_cleaned | set2_cleaned
                        
                        if len(set1_cleaned) == 0:
                            set1_cleaned = {}
                        
                        primary.iloc[inx,1] = list(set1_cleaned)
                
                
        primary.columns = ['parent', 'children']
      
        primary = primary.explode('children')
        
        return primary.reset_index(drop = True)




    def GO_overrepresentation(self):
        
        if "GO-TERM" in self.input_data.keys():
            
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            go1 = pd.DataFrame(self.input_data['GO-TERM']['gene_info'])
            go2 = pd.DataFrame(self.input_data['GO-TERM']['go_names'])
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            
            query = "SELECT * FROM GO_hierarchy;"
            
            hierarchy = pd.read_sql_query(query, conn).applymap(self.deserialize_data)

            results = self.create_full_conections(hierarchy, grade = self.go_grade)
            
            for i in results.index:
                
                if results.iloc[i, 0] == results.iloc[i, 0] and results.iloc[i, 1] != results.iloc[i, 1]:
                    
                    results.iloc[i, 1] = results.iloc[i, 0]
                    
                    results.iloc[i, 0] = 'Core group'
                    
                    
    
            go3 = pd.merge(go1, go2, on = 'GO_id', how = 'left')
            
            results = results[results.iloc[:,0].isin(list(go3['GO_id']) + ['Core group'])]
            results = results[results.iloc[:,1].isin(list(go3['GO_id']))]


            
            del go1, go2
            
            
            go_out = {}
            
              
            go_out['parent'] = []
            go_out['parent_genes'] = []
            go_out['parent_pval_FISH'] = []
            go_out['parent_pval_BIN'] = []
            go_out['parent_n'] = []
            go_out['parent_pct'] = []
    
            go_out['child'] = []
            go_out['child_genes'] = []
            go_out['child_pval_FISH'] = []
            go_out['child_pval_BIN'] = []
            go_out['child_n'] = []
            go_out['child_pct'] = []
    
            
            species = self.input_data['species']['species_study']
            species = ', '.join(map(lambda x: f"'{x}'", species))
            
            
            for i in tqdm(set(results['parent'])):
                
                if i != 'Core group':
            
                    query = f"""SELECT * FROM GO_gene_info 
                            WHERE GO_id IN ('{i}')
                              AND species IN ({species});
                            """
            
                    df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
        
                    res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                K = len(set(df['id'])), 
                                                n = len(set(self.input_data['gene_info']['sid'])), 
                                                k = len(set(go3['id'][go3['GO_id'] == i])))
                    
                    for c in set(results['children'][results['parent'].isin([i])]):
                        
                        query = f"""SELECT * FROM GO_gene_info 
                                WHERE GO_id IN ('{c}')
                                  AND species IN ({species});
                                """
                        df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)

                        res2 = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                    K = len(set(df['id'][df['GO_id'] == c])), 
                                                    n = len(set(self.input_data['gene_info']['sid'])), 
                                                    k = len(set(go3['id'][go3['GO_id'] == c])))
                
                        
                                 
                        go_out['parent'].append(i)
                        go_out['parent_genes'].append(list(set(go3['found_names'][go3['GO_id'] == i])))
                        go_out['parent_pval_FISH'].append(res['fisher_p_value'])
                        go_out['parent_pval_BIN'].append(res['binomial_p_value'])
                        go_out['parent_n'].append(len(set(go3['id'][go3['GO_id'] == i])))
                        go_out['parent_pct'].append(len(set(go3['id'][go3['GO_id'] == i])) / len(set(self.input_data['gene_info']['sid'])))
        
                        go_out['child'].append(c)
                        go_out['child_genes'].append(list(set(go3['found_names'][go3['GO_id'] == c])))
                        go_out['child_pval_FISH'].append(res2['fisher_p_value'])
                        go_out['child_pval_BIN'].append(res2['binomial_p_value'])
                        go_out['child_n'].append(len(set(go3['id'][go3['GO_id'] == c])))
                        go_out['child_pct'].append(len(set(go3['id'][go3['GO_id'] == c])) / len(set(self.input_data['gene_info']['sid'])))
            
                else:
                    
                
                    
                    for c in set(results['children'][results['parent'].isin([i])]):
                        
                        query = f"""SELECT * FROM GO_gene_info 
                                WHERE GO_id IN ('{c}')
                                  AND species IN ({species});
                                """
                        df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
    
                        res2 = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                    K = len(set(df['id'][df['GO_id'] == c])), 
                                                    n = len(set(self.input_data['gene_info']['sid'])), 
                                                    k = len(set(go3['id'][go3['GO_id'] == c])))
                
                        
                                 
                        go_out['parent'].append(i)
                        go_out['parent_genes'].append(list(set(go3['found_names'][go3['GO_id'].isin(set(results['children'][results['parent'] == i]))])))
                        go_out['parent_pval_FISH'].append(0)
                        go_out['parent_pval_BIN'].append(0)
                        go_out['parent_n'].append(len(set(go3['found_names'][go3['GO_id'].isin(set(results['children'][results['parent'] == i]))])))
                        go_out['parent_pct'].append(len(set(go3['found_names'][go3['GO_id'].isin(set(results['children'][results['parent'] == i]))])) / len(set(self.input_data['gene_info']['sid'])))
        
        
                        go_out['child'].append(c)
                        go_out['child_genes'].append(list(set(go3['found_names'][go3['GO_id'] == c])))
                        go_out['child_pval_FISH'].append(res2['fisher_p_value'])
                        go_out['child_pval_BIN'].append(res2['binomial_p_value'])
                        go_out['child_n'].append(len(set(go3['id'][go3['GO_id'] == c])))
                        go_out['child_pct'].append(len(set(go3['id'][go3['GO_id'] == c])) / len(set(self.input_data['gene_info']['sid'])))
        
            
            
    
    
    
            go_out = pd.DataFrame(go_out)
            
            # parent adjustment
            go_out['parent_adj_pval_BIN-BF'] = go_out['parent_pval_BIN'] * len(go_out['parent_pval_BIN'])
            go_out['parent_adj_pval_BIN-BF'][go_out['parent_adj_pval_BIN-BF'] >= 1] = 1
            go_out['parent_adj_pval_FISH-BF'] = go_out['parent_pval_FISH'] * len(go_out['parent_pval_FISH'])
            go_out['parent_adj_pval_FISH-BF'][go_out['parent_adj_pval_FISH-BF'] >= 1] = 1
            
            go_out = go_out.sort_values(by='parent_pval_BIN',  ascending=True)
        
            n = len(go_out['parent_pval_BIN'])
        
            go_out['parent_pval_BIN-BH'] = (go_out['parent_pval_BIN'] * n) / np.arange(1, n+1)
            
            go_out = go_out.sort_values(by='parent_pval_FISH',  ascending=True)
        
            go_out['parent_adj_pval_FISH-BH'] = (go_out['parent_pval_FISH'] * n) / np.arange(1, n+1)
            
            go_out['parent_adj_pval_FISH-BH'][go_out['parent_adj_pval_FISH-BH'] >= 1] = 1
            go_out['parent_pval_BIN-BH'][go_out['parent_pval_BIN-BH'] >= 1] = 1
            
            
            # child adjustment
            go_out['child_adj_pval_BIN-BF'] = go_out['child_pval_BIN'] * len(go_out['child_pval_BIN'])
            go_out['child_adj_pval_BIN-BF'][go_out['child_adj_pval_BIN-BF'] >= 1] = 1
            go_out['child_adj_pval_FISH-BF'] = go_out['child_pval_FISH'] * len(go_out['child_pval_FISH'])
            go_out['child_adj_pval_FISH-BF'][go_out['child_adj_pval_FISH-BF'] >= 1] = 1
            
            go_out = go_out.sort_values(by='child_pval_BIN',  ascending=True)
        
            n = len(go_out['child_pval_BIN'])
        
            go_out['child_pval_BIN-BH'] = (go_out['child_pval_BIN'] * n) / np.arange(1, n+1)
            
            go_out = go_out.sort_values(by='child_pval_FISH',  ascending=True)
        
            go_out['child_adj_pval_FISH-BH'] = (go_out['child_pval_FISH'] * n) / np.arange(1, n+1)
            
            go_out['child_adj_pval_FISH-BH'][go_out['child_adj_pval_FISH-BH'] >= 1] = 1
            go_out['child_pval_BIN-BH'][go_out['child_pval_BIN-BH'] >= 1] = 1
                    
            conn.close()
    
            gn = pd.DataFrame(self.input_data['GO-TERM']['go_names'])

            name_mapping = dict(zip(gn['GO_id'], gn['name']))
            go_out['parent_name'] = go_out['parent'].map(name_mapping)
            go_out['child_name'] = go_out['child'].map(name_mapping)
            
            del gn
            
            self.GO_stat = go_out.to_dict(orient = 'list')
            
        else:
            print("\nGO enrichment analysis could not be performed due to missing GO information in the input data.")



    def KEGG_overrepresentation(self):
        
        if "KEGG" in self.input_data.keys():
            
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            kegg_out = {}
            kegg = pd.DataFrame(self.input_data['KEGG'])
            
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            kegg_out['2nd'] = []
            kegg_out['2nd_genes'] = []
            kegg_out['2nd_pval_FISH'] = []
            kegg_out['2nd_pval_BIN'] = []
            kegg_out['2nd_n'] = []
            kegg_out['2nd_pct'] = []
    
            kegg_out['3rd'] = []
            kegg_out['3rd_genes'] = []
            kegg_out['3rd_pval_FISH'] = []
            kegg_out['3rd_pval_BIN'] = []
            kegg_out['3rd_n'] = []
            kegg_out['3rd_pct'] = []
    
    
    
            for i in tqdm(set(kegg['2nd'])):
            
                query = f'SELECT * FROM KEGG WHERE "2nd" IN ("{i}");'
        
                df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
    
                res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                            K = len(set(df['id'])), 
                                            n = len(set(self.input_data['gene_info']['sid'])), 
                                            k = len(set(kegg['id'][kegg['2nd'] == i])))
                
                for c in set(kegg['3rd'][kegg['2nd'].isin([i])]):
                    
                    res2 = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                K = len(set(df['id'][df['3rd'] == c])), 
                                                n = len(set(self.input_data['gene_info']['sid'])), 
                                                k = len(set(kegg['id'][kegg['3rd'] == c])))
            
                    
                             
                    kegg_out['2nd'].append(i)
                    kegg_out['2nd_genes'].append(list(set(kegg['found_names'][kegg['2nd'] == i])))
                    kegg_out['2nd_pval_FISH'].append(res['fisher_p_value'])
                    kegg_out['2nd_pval_BIN'].append(res['binomial_p_value'])
                    kegg_out['2nd_n'].append(len(set(kegg['id'][kegg['2nd'] == i])))
                    kegg_out['2nd_pct'].append(len(set(kegg['id'][kegg['2nd'] == i])) / len(set(self.input_data['gene_info']['sid'])))
    
                    kegg_out['3rd'].append(c)
                    kegg_out['3rd_genes'].append(list(set(kegg['found_names'][kegg['3rd'] == c])))
                    kegg_out['3rd_pval_FISH'].append(res2['fisher_p_value'])
                    kegg_out['3rd_pval_BIN'].append(res2['binomial_p_value'])
                    kegg_out['3rd_n'].append(len(set(kegg['id'][kegg['3rd'] == c])))
                    kegg_out['3rd_pct'].append(len(set(kegg['id'][kegg['3rd'] == c])) / len(set(self.input_data['gene_info']['sid'])))
    
            kegg_out = pd.DataFrame(kegg_out)
            
            # 2nd adjustment
            kegg_out['2nd_adj_pval_BIN-BF'] = kegg_out['2nd_pval_BIN'] * len(kegg_out['2nd_pval_BIN'])
            kegg_out['2nd_adj_pval_BIN-BF'][kegg_out['2nd_adj_pval_BIN-BF'] >= 1] = 1
            kegg_out['2nd_adj_pval_FISH-BF'] = kegg_out['2nd_pval_FISH'] * len(kegg_out['2nd_pval_FISH'])
            kegg_out['2nd_adj_pval_FISH-BF'][kegg_out['2nd_adj_pval_FISH-BF'] >= 1] = 1
            
            kegg_out = kegg_out.sort_values(by='2nd_pval_BIN',  ascending=True)
        
            n = len(kegg_out['2nd_pval_BIN'])
        
            kegg_out['2nd_pval_BIN-BH'] = (kegg_out['2nd_pval_BIN'] * n) / np.arange(1, n+1)
            
            kegg_out = kegg_out.sort_values(by='2nd_pval_FISH',  ascending=True)
        
            kegg_out['2nd_adj_pval_FISH-BH'] = (kegg_out['2nd_pval_FISH'] * n) / np.arange(1, n+1)
            
            kegg_out['2nd_adj_pval_FISH-BH'][kegg_out['2nd_adj_pval_FISH-BH'] >= 1] = 1
            kegg_out['2nd_pval_BIN-BH'][kegg_out['2nd_pval_BIN-BH'] >= 1] = 1
            
            
            # 3rd adjustment
            kegg_out['3rd_adj_pval_BIN-BF'] = kegg_out['3rd_pval_BIN'] * len(kegg_out['3rd_pval_BIN'])
            kegg_out['3rd_adj_pval_BIN-BF'][kegg_out['3rd_adj_pval_BIN-BF'] >= 1] = 1
            kegg_out['3rd_adj_pval_FISH-BF'] = kegg_out['3rd_pval_FISH'] * len(kegg_out['3rd_pval_FISH'])
            kegg_out['3rd_adj_pval_FISH-BF'][kegg_out['3rd_adj_pval_FISH-BF'] >= 1] = 1
            
            kegg_out = kegg_out.sort_values(by='3rd_pval_BIN',  ascending=True)
        
            n = len(kegg_out['3rd_pval_BIN'])
        
            kegg_out['3rd_pval_BIN-BH'] = (kegg_out['3rd_pval_BIN'] * n) / np.arange(1, n+1)
            
            kegg_out = kegg_out.sort_values(by='3rd_pval_FISH',  ascending=True)
        
            kegg_out['3rd_adj_pval_FISH-BH'] = (kegg_out['3rd_pval_FISH'] * n) / np.arange(1, n+1)
            
            kegg_out['3rd_adj_pval_FISH-BH'][kegg_out['3rd_adj_pval_FISH-BH'] >= 1] = 1
            kegg_out['3rd_pval_BIN-BH'][kegg_out['3rd_pval_BIN-BH'] >= 1] = 1
                    
            conn.close()
            
            self.KEGG_stat = kegg_out.to_dict(orient = 'list')
            
        else:
            print("\nKEGG enrichment analysis could not be performed due to missing KEGG information in the input data.")

        
        
        
        
    def REACTOME_overrepresentation(self):
        
        if "REACTOME" in self.input_data.keys():
            
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            reactome_out = {}
            reactome = pd.DataFrame(self.input_data['REACTOME'])
            
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            reactome_out['pathway'] = []
            reactome_out['pathway_genes'] = []
            reactome_out['complex'] = []
            reactome_out['pathway_pval_FISH'] = []
            reactome_out['pathway_pval_BIN'] = []
            reactome_out['pathway_n'] = []
            reactome_out['pathway_pct'] = []
    
            reactome_out['top_level_pathway'] = []
            reactome_out['top_level_pathway_genes'] = []
            reactome_out['top_level_pathway_pval_FISH'] = []
            reactome_out['top_level_pathway_pval_BIN'] = []
            reactome_out['top_level_pathway_n'] = []
            reactome_out['top_level_pathway_pct'] = []
    
    
    
            for i in tqdm(set(reactome['top_level_pathway'])):
            
                query = f'SELECT * FROM REACTOME WHERE "top_level_pathway" IN ("{i}");'
        
                df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
    
                res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                            K = len(set(df['id'])), 
                                            n = len(set(self.input_data['gene_info']['sid'])), 
                                            k = len(set(reactome['id'][reactome['top_level_pathway'] == i])))
                
                for c in set(reactome['pathway'][reactome['top_level_pathway'].isin([i])]):
                    
                    res2 = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                K = len(set(df['id'][df['pathway'] == c])), 
                                                n = len(set(self.input_data['gene_info']['sid'])), 
                                                k = len(set(reactome['id'][reactome['pathway'] == c])))
            
                    
                             
                    reactome_out['top_level_pathway'].append(i)
                    reactome_out['top_level_pathway_genes'].append(list(set(reactome['found_names'][reactome['top_level_pathway'] == i])))
                    reactome_out['top_level_pathway_pval_FISH'].append(res['fisher_p_value'])
                    reactome_out['top_level_pathway_pval_BIN'].append(res['binomial_p_value'])
                    reactome_out['top_level_pathway_n'].append(len(set(reactome['id'][reactome['top_level_pathway'] == i])))
                    reactome_out['top_level_pathway_pct'].append(len(set(reactome['id'][reactome['top_level_pathway'] == i])) / len(set(self.input_data['gene_info']['sid'])))
    
                    reactome_out['pathway'].append(c)
                    reactome_out['pathway_genes'].append(list(set(reactome['found_names'][reactome['pathway'] == c])))
                    reactome_out['complex'].append(list(set(reactome['complex'][reactome['pathway'] == c])))
                    reactome_out['pathway_pval_FISH'].append(res2['fisher_p_value'])
                    reactome_out['pathway_pval_BIN'].append(res2['binomial_p_value'])
                    reactome_out['pathway_n'].append(len(set(reactome['id'][reactome['pathway'] == c])))
                    reactome_out['pathway_pct'].append(len(set(reactome['id'][reactome['pathway'] == c])) / len(set(self.input_data['gene_info']['sid'])))
    
            reactome_out = pd.DataFrame(reactome_out)
            
            # top_level_pathway adjustment
            reactome_out['top_level_pathway_adj_pval_BIN-BF'] = reactome_out['top_level_pathway_pval_BIN'] * len(reactome_out['top_level_pathway_pval_BIN'])
            reactome_out['top_level_pathway_adj_pval_BIN-BF'][reactome_out['top_level_pathway_adj_pval_BIN-BF'] >= 1] = 1
            reactome_out['top_level_pathway_adj_pval_FISH-BF'] = reactome_out['top_level_pathway_pval_FISH'] * len(reactome_out['top_level_pathway_pval_FISH'])
            reactome_out['top_level_pathway_adj_pval_FISH-BF'][reactome_out['top_level_pathway_adj_pval_FISH-BF'] >= 1] = 1
            
            reactome_out = reactome_out.sort_values(by='top_level_pathway_pval_BIN',  ascending=True)
        
            n = len(reactome_out['top_level_pathway_pval_BIN'])
        
            reactome_out['top_level_pathway_pval_BIN-BH'] = (reactome_out['top_level_pathway_pval_BIN'] * n) / np.arange(1, n+1)
            
            reactome_out = reactome_out.sort_values(by='top_level_pathway_pval_FISH',  ascending=True)
        
            reactome_out['top_level_pathway_adj_pval_FISH-BH'] = (reactome_out['top_level_pathway_pval_FISH'] * n) / np.arange(1, n+1)
            
            reactome_out['top_level_pathway_adj_pval_FISH-BH'][reactome_out['top_level_pathway_adj_pval_FISH-BH'] >= 1] = 1
            reactome_out['top_level_pathway_pval_BIN-BH'][reactome_out['top_level_pathway_pval_BIN-BH'] >= 1] = 1
            
            
            # pathway adjustment
            reactome_out['pathway_adj_pval_BIN-BF'] = reactome_out['pathway_pval_BIN'] * len(reactome_out['pathway_pval_BIN'])
            reactome_out['pathway_adj_pval_BIN-BF'][reactome_out['pathway_adj_pval_BIN-BF'] >= 1] = 1
            reactome_out['pathway_adj_pval_FISH-BF'] = reactome_out['pathway_pval_FISH'] * len(reactome_out['pathway_pval_FISH'])
            reactome_out['pathway_adj_pval_FISH-BF'][reactome_out['pathway_adj_pval_FISH-BF'] >= 1] = 1
            
            reactome_out = reactome_out.sort_values(by='pathway_pval_BIN',  ascending=True)
        
            n = len(reactome_out['pathway_pval_BIN'])
        
            reactome_out['pathway_pval_BIN-BH'] = (reactome_out['pathway_pval_BIN'] * n) / np.arange(1, n+1)
            
            reactome_out = reactome_out.sort_values(by='pathway_pval_FISH',  ascending=True)
        
            reactome_out['pathway_adj_pval_FISH-BH'] = (reactome_out['pathway_pval_FISH'] * n) / np.arange(1, n+1)
            
            reactome_out['pathway_adj_pval_FISH-BH'][reactome_out['pathway_adj_pval_FISH-BH'] >= 1] = 1
            reactome_out['pathway_pval_BIN-BH'][reactome_out['pathway_pval_BIN-BH'] >= 1] = 1
                
            
            conn.close()
            
            self.REACTOME_stat = reactome_out.to_dict(orient = 'list')
        
        else:
            print("\nREACTOME enrichment analysis could not be performed due to missing REACTOME information in the input data.")

        
        
    def DISEASES_overrepresentation(self):
        
        
        if "DISEASES" in self.input_data.keys():
        
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            diseases_out = {}
            diseases = pd.DataFrame(self.input_data['DISEASES'])
           
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            diseases_out['disease'] = []
            diseases_out['genes'] = []
            diseases_out['pval_FISH'] = []
            diseases_out['pval_BIN'] = []
            diseases_out['n'] = []
            diseases_out['pct'] = []
            
    
            for i in tqdm(set(diseases['disease'])):           
        
                
                query = f'SELECT * FROM disease WHERE disease IN ("{i}");'
    
                df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
     
                res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                            K = len(set(df['id'][df['disease'] == i])), 
                                            n = len(set(self.input_data['gene_info']['sid'])), 
                                            k = len(set(diseases['id'][diseases['disease'] == i])))
                
                diseases_out['disease'].append(i)
                diseases_out['genes'].append(list(set(diseases['found_names'][diseases['disease'].isin([i])])))
                diseases_out['pval_FISH'].append(res['fisher_p_value'])
                diseases_out['pval_BIN'].append(res['binomial_p_value'])
                diseases_out['n'].append(len(set(diseases['id'][diseases['disease'] == i])))
                diseases_out['pct'].append(len(set(diseases['id'][diseases['disease'] == i])) / len(set(self.input_data['gene_info']['sid'])))
          
            diseases_out = pd.DataFrame(diseases_out)
            
            diseases_out['adj_pval_BIN-BF'] = diseases_out['pval_BIN'] * len(diseases_out['pval_BIN'])
            diseases_out['adj_pval_BIN-BF'][diseases_out['adj_pval_BIN-BF'] >= 1] = 1
            diseases_out['adj_pval_FISH-BF'] = diseases_out['pval_FISH'] * len(diseases_out['pval_FISH'])
            diseases_out['adj_pval_FISH-BF'][diseases_out['adj_pval_FISH-BF'] >= 1] = 1
            
            diseases_out = diseases_out.sort_values(by='pval_BIN',  ascending=True)
        
            n = len(diseases_out['pval_BIN'])
        
            diseases_out['pval_BIN-BH'] = (diseases_out['pval_BIN'] * n) / np.arange(1, n+1)
            
            diseases_out = diseases_out.sort_values(by='pval_FISH',  ascending=True)
        
            diseases_out['adj_pval_FISH-BH'] = (diseases_out['pval_FISH'] * n) / np.arange(1, n+1)
            
            diseases_out['adj_pval_FISH-BH'][diseases_out['adj_pval_FISH-BH'] >= 1] = 1
            diseases_out['pval_BIN-BH'][diseases_out['pval_BIN-BH'] >= 1] = 1
                  
            conn.close()
    
            self.DISEASE_stat = diseases_out.to_dict(orient = 'list')
            
        else:
            print("\nDISEASES enrichment analysis could not be performed due to missing DISEASES information in the input data.")

            
        
        
    def ViMIC_overrepresentation(self):
        
        if "ViMIC" in self.input_data.keys():
        
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            vimic_out = {}
            vimic = pd.DataFrame(self.input_data['ViMIC'])
           
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            vimic_out['virus'] = []
            vimic_out['genes'] = []
            vimic_out['group'] = []
    
            vimic_out['pval_FISH'] = []
            vimic_out['pval_BIN'] = []
            vimic_out['n'] = []
            vimic_out['pct'] = []
            
    
            for i in tqdm(set(vimic['virus'])):           
        
                query = f'SELECT * FROM ViMIC WHERE virus IN ("{i}");'
    
                df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
     
                res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                            K = len(set(df['id'][df['virus'] == i])), 
                                            n = len(set(self.input_data['gene_info']['sid'])), 
                                            k = len(set(vimic['id'][vimic['virus'] == i])))
                
                vimic_out['virus'].append(i)
                vimic_out['genes'].append(list(set(vimic['found_names'][vimic['virus'].isin([i])])))
                vimic_out['group'].append(list(set(vimic['group'][vimic['virus'].isin([i])])))
                vimic_out['pval_FISH'].append(res['fisher_p_value'])
                vimic_out['pval_BIN'].append(res['binomial_p_value'])
                vimic_out['n'].append(len(set(vimic['id'][vimic['virus'] == i])))
                vimic_out['pct'].append(len(set(vimic['id'][vimic['virus'] == i])) / len(set(self.input_data['gene_info']['sid'])))
          
            vimic_out = pd.DataFrame(vimic_out)
            
            vimic_out['adj_pval_BIN-BF'] = vimic_out['pval_BIN'] * len(vimic_out['pval_BIN'])
            vimic_out['adj_pval_BIN-BF'][vimic_out['adj_pval_BIN-BF'] >= 1] = 1
            vimic_out['adj_pval_FISH-BF'] = vimic_out['pval_FISH'] * len(vimic_out['pval_FISH'])
            vimic_out['adj_pval_FISH-BF'][vimic_out['adj_pval_FISH-BF'] >= 1] = 1
            
            vimic_out = vimic_out.sort_values(by='pval_BIN',  ascending=True)
        
            n = len(vimic_out['pval_BIN'])
        
            vimic_out['pval_BIN-BH'] = (vimic_out['pval_BIN'] * n) / np.arange(1, n+1)
            
            vimic_out = vimic_out.sort_values(by='pval_FISH',  ascending=True)
        
            vimic_out['adj_pval_FISH-BH'] = (vimic_out['pval_FISH'] * n) / np.arange(1, n+1)
            
            vimic_out['adj_pval_FISH-BH'][vimic_out['adj_pval_FISH-BH'] >= 1] = 1
            vimic_out['pval_BIN-BH'][vimic_out['pval_BIN-BH'] >= 1] = 1
                  
            conn.close()
    
            self.ViMIC_stat = vimic_out.to_dict(orient = 'list')
        
        else:
            print("\nViMIC enrichment analysis could not be performed due to missing ViMIC information in the input data.")

            
        
        
    def features_specificity(self):
        
        if "HPA" in self.input_data.keys():
            
            if self.occ == None:
                with open(os.path.join(self.path_in_inside, 'occ_dict.json'), 'r') as json_file:
                    self.occ = (json.load(json_file))
            
            HPA_out = {}
            
            HPA = self.input_data['HPA']
            
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
    
    
            for k in HPA.keys():
                if k != 'HPA_blood_markers':
                    tmp = pd.DataFrame(HPA[k])
                    
                    tmp_out = {}
                    
                    tmp_out['specificity'] = []
                    tmp_out['genes'] = []
                    tmp_out['pval_FISH'] = []
                    tmp_out['pval_BIN'] = []
                    tmp_out['n'] = []
                    tmp_out['pct'] = []
                    
                    col_name = None
                    
                    if 'location' in list(tmp.columns):
                        col_name = 'location'
                    else:
                        col_name = 'name'
    
                    for i in tqdm(set(tmp[col_name])):      
                        
                        query = f'SELECT * FROM {k} WHERE {col_name} IN ("{i}");'
    
                        df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
                        
                        
             
                        res = self.run_enrichment_tests(N = self.occ['Genes_Homo_sapiens'], 
                                                    K = len(set(df['id'][df[col_name] == i])), 
                                                    n = len(set(self.input_data['gene_info']['sid'])), 
                                                    k = len(set(tmp['id'][tmp[col_name] == i])))
                        
                        tmp_out['specificity'].append(i)
                        tmp_out['genes'].append(list(set(tmp['found_names'][tmp[col_name].isin([i])])))
                        tmp_out['pval_FISH'].append(res['fisher_p_value'])
                        tmp_out['pval_BIN'].append(res['binomial_p_value'])
                        tmp_out['n'].append(len(set(tmp['id'][tmp[col_name] == i])))
                        tmp_out['pct'].append(len(set(tmp['id'][tmp[col_name] == i])) / len(set(self.input_data['gene_info']['sid'])))
                  
                    tmp_out = pd.DataFrame(tmp_out)
                    
                    tmp_out['adj_pval_BIN-BF'] = tmp_out['pval_BIN'] * len(tmp_out['pval_BIN'])
                    tmp_out['adj_pval_BIN-BF'][tmp_out['adj_pval_BIN-BF'] >= 1] = 1
                    tmp_out['adj_pval_FISH-BF'] = tmp_out['pval_FISH'] * len(tmp_out['pval_FISH'])
                    tmp_out['adj_pval_FISH-BF'][tmp_out['adj_pval_FISH-BF'] >= 1] = 1
                    
                    tmp_out = tmp_out.sort_values(by='pval_BIN',  ascending=True)
                
                    n = len(tmp_out['pval_BIN'])
                
                    tmp_out['pval_BIN-BH'] = (tmp_out['pval_BIN'] * n) / np.arange(1, n+1)
                    
                    tmp_out = tmp_out.sort_values(by='pval_FISH',  ascending=True)
                
                    tmp_out['adj_pval_FISH-BH'] = (tmp_out['pval_FISH'] * n) / np.arange(1, n+1)
                    
                    tmp_out['adj_pval_FISH-BH'][tmp_out['adj_pval_FISH-BH'] >= 1] = 1
                    tmp_out['pval_BIN-BH'][tmp_out['pval_BIN-BH'] >= 1] = 1
                    
                    HPA_out[k] = tmp_out.to_dict(orient = 'list')
    
                          
            conn.close()
            
            self.specificity_stat = HPA_out
               
        else:
            print("\nHPA enrichment analysis could not be performed due to missing HPA information in the input data.")

                
               
    def gene_interaction(self):
        
 
        interaction_mapping = {
            "coexpression": ["gene -> gene"],
            "coexpression_transferred": ["gene -> gene"],
            "cooccurence": ["gene -> gene"],
            "database": ["protein -> protein"],
            "database_transferred": ["protein -> protein"],
            "experiments": ["protein -> protein"],
            "experiments_transferred": ["protein -> protein"],
            "fusion": ["gene -> gene"],
            "homology": ["protein -> protein"],
            "neighborhood": ["gene -> gene"],
            "neighborhood_transferred": ["gene -> gene"],
            "textmining": ["protein -> protein"],
            "textmining_transferred": ["protein -> protein"]
            }
        
        interactome = {}
    
        if "IntAct" in self.input_data.keys() and "STRING" in self.input_data.keys():
            
            ia = pd.DataFrame(self.input_data['IntAct']['gene_products'])
            ia = ia[ia['source'].isin(self.interaction_source)]
            ia = ia[ia['species'].isin(self.species_study)]

            
            interactome['A'] = list(ia['found_names_1'])
            interactome['B'] = list(ia['found_names_2'])
            interactome['interaction_type'] = list(ia['interaction_type'])
            interactome['connection_type'] = [f"{a} -> {b}" for a, b in zip(ia['interactor_type_1'], ia['interactor_type_2'])]
            interactome['source'] = list(ia['source'])
            
            
            ia = pd.DataFrame(self.input_data['STRING'])
            ia = ia[ia['species'].isin(self.species_study)]
            ia = ia[ia['combined_score'] >= self.interaction_strength]

            
            interactome['A'] = list(interactome['A']) + list(ia['found_names_1'])
            interactome['B'] = list(interactome['B']) + list(ia['found_names_2'])
            interactome['source'] = list(interactome['source']) + ['STRING'] * len(list(ia['source']))
            interactome['interaction_type'] = list(interactome['interaction_type']) + [None] * len(list(ia['source']))
            interactome['connection_type'] = list(interactome['connection_type']) + [self.map_interactions_flat(x, interaction_mapping) for x in ia['source']]

            interactome = pd.DataFrame(interactome)
            interactome = interactome[interactome['source'].isin(self.interaction_source)]
            interactome = interactome.explode('connection_type')
            interactome = interactome.to_dict(orient = 'list')
            
            self.features_interactions = interactome
            
        elif "IntAct" in self.input_data.keys():
            
            ia = pd.DataFrame(self.input_data['IntAct']['gene_products'])
            ia = ia[ia['source'].isin(self.interaction_source)]
            ia = ia[ia['species'].isin(self.species_study)]

            
            interactome['A'] = list(ia['found_names_1'])
            interactome['B'] = list(ia['found_names_2'])
            interactome['interaction_type'] = list(ia['interaction_type'])
            interactome['connection_type'] = [f"{a} -> {b}" for a, b in zip(ia['interactor_type_1'], ia['interactor_type_2'])]
            interactome['source'] = list(ia['source'])
            
            
            self.features_interactions = interactome
            
        
        elif "STRING" in self.input_data.keys():
            
            ia = pd.DataFrame(self.input_data['STRING'])
            ia = ia[ia['species'].isin(self.species_study)]
            ia = ia[ia['combined_score'] >= self.interaction_strength]

            
            interactome['A'] = list(ia['found_names_1'])
            interactome['B'] = list(ia['found_names_2'])
            interactome['source'] = ['STRING'] * len(list(ia['source']))
            interactome['interaction_type'] = [None] * len(list(ia['source']))
            interactome['connection_type'] = [self.map_interactions_flat(x, interaction_mapping) for x in ia['source']]
            
            interactome = pd.DataFrame(interactome)
            interactome = interactome[interactome['source'].isin(self.interaction_source)]
            interactome = interactome.explode('connection_type')
            interactome = interactome.to_dict(orient = 'list')
            
            self.features_interactions = interactome


            
        else:
            print("\nGene interaction analysis could not be performed due to missing STRING/IntAct information in the input data.")

                 

    
    
  
    def GO_network(self):
        
        if "GO-TERM" in self.input_data.keys():
        
            goh = pd.DataFrame(self.input_data['GO-TERM']['hierarchy'])
            
            relation_colors = {
                'is_a_ids': 'blue',
                'part_of_ids': 'orange',
                'has_part_ids': 'green',
                'regulates_ids': 'red',
                'negatively_regulates_ids': 'purple',
                'positively_regulates_ids': 'brown'
            }
            
            full_df = pd.DataFrame()
            for i in ['is_a_ids', 'part_of_ids', 'has_part_ids', 'regulates_ids',
                   'negatively_regulates_ids', 'positively_regulates_ids']:
                tmp = goh[[i, 'GO_id']]
                tmp = tmp[tmp[i].notna()]
                tmp.columns = ['parent', 'children']
                tmp['color'] = relation_colors[i]
                if len(tmp.index) > 0:
                    full_df = pd.concat([full_df, tmp])
           
            
            goh = pd.DataFrame(self.input_data['GO-TERM']['gene_info'])
    
            full_df = pd.merge(full_df, goh[['GO_id', 'found_names']], left_on = 'parent', right_on= 'GO_id',  how='left')
            full_df.pop('GO_id')
            
            full_df = pd.merge(full_df, goh[['GO_id', 'found_names']], left_on = 'children', right_on= 'GO_id',  how='left')
            full_df.pop('GO_id')
            
            full_df = full_df.reset_index(drop = True)
            
            full_df['features'] = full_df[['found_names_x', 'found_names_y']].apply(lambda row: list(set(row)), axis=1)
            
            full_df.pop('found_names_x')
            full_df.pop('found_names_y')
            
            del goh
            
            full_df = full_df.explode('features')
    
    
            go_out = pd.DataFrame(self.GO_stat)
            
            test_col = self.select_test(self.network_stat['test'], self.network_stat['adj'])
            
            go_out_parent = go_out[go_out[f'parent_{test_col}'] <= self.network_stat['p_val']]
            go_out_children = go_out[go_out[f'child_{test_col}'] <= self.network_stat['p_val']]

            full_names = list(set(go_out_parent['parent'])) + list(set(go_out_children['child']))
            
            full_df = full_df[full_df['parent'].isin(full_names)]
            full_df = full_df[full_df['children'].isin(full_names)]
            
            
            gn = pd.DataFrame(self.input_data['GO-TERM']['go_names'])

            name_mapping = dict(zip(gn['GO_id'], gn['name']))
            full_df['parent'] = full_df['parent'].map(name_mapping)
            full_df['children'] = full_df['children'].map(name_mapping)
            
            
            self.GO_net = full_df.to_dict(orient = 'list')

        else:
            print("\nMissing GO information in the input data.")

          
        

    def KEGG_network(self):
        
         
        if "KEGG" in self.input_data.keys():
        
            kegg = pd.DataFrame(self.input_data['KEGG'])
            
            full_df = pd.DataFrame()
            
            full_df['parent'] = kegg['2nd']
            full_df['children'] = kegg['3rd']
            full_df['features'] = kegg['found_names']
    
    
            kegg_out = pd.DataFrame(self.KEGG_stat)
            
            test_col = self.select_test(self.network_stat['test'], self.network_stat['adj'])
            
            kegg_out = kegg_out[kegg_out[f'3rd_{test_col}'] <= self.network_stat['p_val']]
            kegg_out = kegg_out[kegg_out[f'2nd_{test_col}'] <= self.network_stat['p_val']]
    
            
            full_df = full_df[full_df['parent'].isin(kegg_out['2nd'])]
            full_df = full_df[full_df['children'].isin(kegg_out['3rd'])]
            
            full_df['color'] = 'gold'
            
            self.KEGG_net = full_df.to_dict(orient = 'list')
        
        else:
            print("\nMissing KEGG information in the input data.")

    
    
    def REACTOME_network(self):
        
        if "REACTOME" in self.input_data.keys():

            reactome = pd.DataFrame(self.input_data['REACTOME'])
            
            
            full_df = pd.DataFrame()
    
            full_df['parent'] = reactome['top_level_pathway']
            full_df['children'] = reactome['pathway']
            full_df['features'] = reactome['found_names']
    
    
    
            reactome_out = pd.DataFrame(self.REACTOME_stat)
            
            test_col = self.select_test(self.network_stat['test'], self.network_stat['adj'])
            
            reactome_out = reactome_out[reactome_out[f'top_level_pathway_{test_col}'] <= self.network_stat['p_val']]
            reactome_out = reactome_out[reactome_out[f'pathway_{test_col}'] <= self.network_stat['p_val']]
    
            
            full_df = full_df[full_df['parent'].isin(reactome_out['top_level_pathway'])]
            full_df = full_df[full_df['children'].isin(reactome_out['pathway'])]
            
            full_df['color'] = 'silver'
            
            self.REACTOME_net = full_df.to_dict(orient = 'list')
        
        else:
            print("\nMissing REACTOME information in the input data.")



    # def DISEASES_network(self):
        
    #     if self.DISEASE_stat != None:
            
            
    #         full_df = pd.DataFrame()
    
    
    #         diseases_out = pd.DataFrame(self.DISEASE_stat)
            
    #         test_col = self.select_test(self.network_stat['test'], self.network_stat['adj'])
            
    #         diseases_out = diseases_out[diseases_out[test_col] <= self.network_stat['p_val']]
            
    
    
    #         full_df['parent'] = diseases_out['disease']
    #         full_df['features'] = diseases_out['genes']
    #         full_df = full_df.explode('features')
            
    #         full_df = pd.merge(full_df, full_df, left_on = 'features' , right_on = 'features' , how = 'left')

    #         full_df.columns = ['parent', 'features', 'children']
            
    #         full_df = full_df[full_df['parent'] != full_df['children']]
            
    #         full_df['color'] = 'yellow'
            
    #         self.DISEASES_net = full_df.to_dict(orient = 'list')
        
    #     else:
    #         print("\nMissing DISEAES overrepresentation  analysis!")
       
        
       
    
    # def ViMIC_network(self):
        
    #     if self.ViMIC_stat != None:
            
            
    #         full_df = pd.DataFrame()
    
    #         vimic_out = pd.DataFrame(self.ViMIC_stat)
            
    #         test_col = self.select_test(self.network_stat['test'], self.network_stat['adj'])
            
    #         vimic_out = vimic_out[vimic_out[test_col] <= self.network_stat['p_val']]
            
    
    
    #         full_df['parent'] = vimic_out['virus']
    #         full_df['features'] = vimic_out['genes']
    #         full_df = full_df.explode('features')
            
    #         full_df = pd.merge(full_df, full_df, left_on = 'features' , right_on = 'features' , how = 'left')

    #         full_df.columns = ['parent', 'features', 'children']
            
    #         full_df['color'] = 'lightyellow'
            
    #         full_df = full_df[full_df['parent'] != full_df['children']]

            
    #         self.ViMIC_net = full_df.to_dict(orient = 'list')
        
    #     else:
    #         print("\nMissing DISEAE overrepresentation  analysis!")
       
       
        
       
    def full_analysis(self):
        
        print('\nGO-TERM overrepresentation analysis...')
        self.GO_overrepresentation()
        
        print('\nKEGG overrepresentation analysis...')
        self.KEGG_overrepresentation()
        
        print('\nREACTOME overrepresentation analysis...')
        self.REACTOME_overrepresentation()
        
        print('\nViMIC overrepresentation analysis...')
        self.ViMIC_overrepresentation()
        
        print('\nDISEASES overrepresentation analysis...')
        self.DISEASES_overrepresentation()
        
        print('\nSpecificity overrepresentation analysis...')
        self.features_specificity()
        
        print('\nInteraction analysis...')
        self.gene_interaction()
        
        print('\nNetwork creating...')
        self.REACTOME_network()
        self.KEGG_network()
        self.GO_network()
    
        
        print('\nComplete!')
       
        
    @property
    def get_full_results(self):
        
        full_results = {}
        
        full_results['enrichment'] = self.input_data
        
        stats = {}
        networks = {}
        
        try:
            stats['KEGG'] = self.get_KEGG_statistics
        except:
            pass
        
        try:
            stats['REACTOME'] = self.get_REACTOME_statistics
        except:
            pass
        
        try:
            stats['GO-TERM'] = self.get_GO_statistics
        except:
            pass
        
        
        try:
            stats['ViMIC'] = self.get_ViMIC_statistics
        except:
            pass
        
        try:
            stats['DISEASES'] = self.get_DISEASE_statistics
        except:
            pass
        
        
        try:
            stats['specificity'] = self.get_specificity_statistics
        except:
            pass
        
        try:
            stats['interactions'] = self.get_features_interactions_statistics
        except:
            pass
        
        try:
            networks['KEGG'] = self.get_KEGG_network
        except:
            pass
        
        try:
            networks['REACTOME'] = self.get_REACTOME_network
        except:
            pass
        
        try:
            networks['GO-TERM'] = self.get_GO_network
        except:
            pass
        
      
        
    
        full_results['statistics'] = stats 
        
        full_results['statistics']['setup'] = {}
        
        full_results['statistics']['setup']['network_stat'] = self.network_stat
        
        full_results['statistics']['setup']['go_grade'] = self.go_grade
        
        full_results['statistics']['setup']['interaction_strength'] = self.interaction_strength

        full_results['statistics']['setup']['interaction_source'] = self.interaction_source
        
        
        full_results['networks'] = networks 
        
        
        
        return full_results
        
            
        
        
        

        

        

        

        
  # """
  # This function creates a bar plot of statistical / overrepresentation analysis terms from GO-TERM, PATHWAYS, DISEASES, and VIRAL (diseases related to viral infection)

  # Args:
  #     GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
  #     GOPa_metadata (dict) - metadata from load_GOPa_meta function
  #     p_val (float) - value of minimal p_val for statistical test
  #     test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
  #     adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
  #     n (int) - maximal number of bars on the graph
  #     side (str) - orientation of bars ['left' / 'right']
  #     color (str)- color of the bars
  #     width (float) - width of the graph
  #     bar_width (float) - width of the bars
  #     count_type (str) - type of amount of term representation on bars ['perc' - percent representation / 'p_val' - -log(p-val) of selected statistic test / 'num' - number representation]
  #     details (float) - degree detail for display GO-TERM and PATHWAYS where 1 is the maximal value and 0.1 is minimal [0.1-1]
  #     omit (list) - type of terms to omit in the graph eg. ['GO-TERM', 'KEGG', 'REACTOME', 'DISEASES', 'VIRAL'] or None
   
  # Returns:
  #     dict of graphs: Dictionary of bar plots of overrepresentation analysis for GO-TERMS, PATHWAYS, DISEASES, and VIRAL
  # """
  
       
      
        
        
class Visualization():

    def __init__(self, input_data:dict):
        
        self.input_data = input_data
               
        
        super().__init__()


    def select_test(self, test, adj):
        try:
            test_string = ''
    
            if adj != None and adj.upper() in ['BF','BH']:
                test_string = test_string + 'adj_pval_'
            else:
                test_string = test_string + 'pval_'
            
            
            if test != None and test.upper() == 'BIN':
                test_string = test_string + 'BIN'
            elif test != None and test.upper() == 'FISH':
                test_string = test_string + 'FISH'
            else:
                test_string = test_string + 'BIN'
                
            
            if adj != None and adj.upper() == 'BF':
                test_string = test_string + '-BF'
            elif adj != None and adj.upper() == 'BH':
                test_string = test_string + '-BH'
            else:
                test_string = test_string + ''
            
            return test_string
        except:
            print('\n')
            print('Provided wrong test input!')
            
      
    def bar_plot(self,
        data,
        n=25, 
        side='right', 
        color='blue', 
        width=10, 
        bar_width=0.5, 
        stat='p_val', 
        sets='GO-TERM',
        column='name',
        x_max=None,
        show_axis=True,
        title=None,
        ax=None  # Dodano argument ax
    ):
        
        tmp = pd.DataFrame(data)
        
        # Sortowanie i wybór danych
        if stat.upper() == 'perc'.upper():
            tmp = tmp.sort_values(by='n', ascending=False).reset_index(drop=True).iloc[0:n, :]
            x_label = 'Percent of genes [%]'
            values = tmp['pct']
        elif stat.upper() == 'p_val'.upper():
            tmp = tmp.sort_values(by='-log(p-val)', ascending=False).reset_index(drop=True).iloc[0:n, :]
            x_label = '-log(p-val)'
            values = tmp['-log(p-val)']
        else:
            tmp = tmp.sort_values(by='n', ascending=False).reset_index(drop=True).iloc[0:n, :]
            x_label = 'Number of genes'
            values = tmp['n']
        
        # Domyślnie twórz nową figurę, jeśli ax nie jest podany
        if ax is None:
            fig_1, ax = plt.subplots(figsize=(width, float(len(tmp[column]) / 2.5)))
        
        # Tworzenie wykresu
        ax.barh(tmp[column], values, color=color, height=bar_width)
        
        # Oś X
        if show_axis:
            ax.set_xlabel(x_label)
        else:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        
        # Oś Y
        ax.set_ylabel('')
        ax.invert_yaxis()
        
        # Tytuł
        if title:
            ax.set_title(title)
        
        # Ukrywanie zbędnych osi
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # Maksymalna wartość X
        if x_max is not None:
            ax.set_xlim(0, x_max)
        
        # Ustawienie położenia etykiet Y
        if side == 'right':
            ax.yaxis.tick_right()
        elif side == 'left':
            ax.invert_xaxis()
        
        try:
            return fig_1
        except:
            return ax

    

    def bar_plot_blood(self,
        data,
        n=25, 
        side='right', 
        color='red', 
        width=10, 
        bar_width=0.5, 
        stat=None, 
        sets=None,
        column=None,
        x_max=None,
        show_axis=True,
        title=None,
        ax=None  
    ):
        
        tmp = pd.DataFrame(data)
        

        tmp = tmp.sort_values(by=stat, ascending=False).reset_index(drop=True).iloc[0:n, :]
        x_label = stat
        values = tmp[stat]
      
        
        # Domyślnie twórz nową figurę, jeśli ax nie jest podany
        if ax is None:
            fig_1, ax = plt.subplots(figsize=(width, float(len(tmp[column]) / 2.5)))
        
        # Tworzenie wykresu
        ax.barh(tmp[column], values, color=color, height=bar_width)
        
        # Oś X
        if show_axis:
            ax.set_xlabel(x_label)
        else:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        
        # Oś Y
        ax.set_ylabel('')
        ax.invert_yaxis()
        
        # Tytuł
        if title:
            ax.set_title(title)
        
        # Ukrywanie zbędnych osi
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # Maksymalna wartość X
        if x_max is not None:
            ax.set_xlim(0, x_max)
        
        # Ustawienie położenia etykiet Y
        if side == 'right':
            ax.yaxis.tick_right()
        elif side == 'left':
            ax.invert_xaxis()
        
        try:
            return fig_1
        except:
            return ax

  


    def GO_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 25, 
                side = 'right', 
                color = 'blue', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
                ):
        
        
        sets = 'GO-TERM'
        column = 'child_name'
        
        
        test_string = self.select_test(test, adj)

        tmp_in = pd.DataFrame(self.input_data['statistics'][sets])
        tmp_in = tmp_in[tmp_in[f'parent_{test_string}'] <= p_val]
        tmp_in = tmp_in[tmp_in[f'child_{test_string}'] <= p_val]
        tmp_in[f'child_{test_string}'] = tmp_in[f'child_{test_string}'] + np.min(tmp_in[f'child_{test_string}'][tmp_in[f'child_{test_string}'] != 0])/2
        tmp_in['-log(p-val)'] = -np.log(tmp_in[f'child_{test_string}'])
        tmp_in = tmp_in.reset_index(drop = True)
        
        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(tmp_in['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(tmp_in['-log(p-val)'])
            
        else:
            x_max = np.max(tmp_in['n'])

            
       
        iter_set = list(set(tmp_in['parent_name']))
        
        hlist = []
        for it in iter_set:
            inx = [x for x in tmp_in.index if it in tmp_in['parent_name'][x]] 
            tmp = tmp_in.loc[inx]
            if float(len(tmp[column])) > n:
                tn = n
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))
                else:
                    hlist.append(tn/2.1)

            else:
                tn = float(len(tmp[column]))
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))

                else:
                    hlist.append(tn/2.1)
                
        
        # hlist = [x for x in hlist]
        fig = plt.figure(figsize=(width, sum(hlist))) 

        gs = GridSpec(len(hlist), 1, height_ratios=hlist)  
        
        
        gs.update(hspace=len(hlist)/50)

        for l, i in enumerate(iter_set):
            inx = [x for x in tmp_in.index if i in tmp_in['parent_name'][x]] 
            tmp = tmp_in.loc[inx]
            
            ax = fig.add_subplot(gs[l])
            
            show_axis = (l + 1 == len(iter_set))  
            
            self.bar_plot(
                data=tmp,
                n=n, 
                side=side, 
                color='blue', 
                width=width, 
                bar_width=bar_width, 
                stat=stat, 
                sets='GO-TERM',
                column=column,
                x_max=x_max,
                show_axis=show_axis,
                title=i,  
                ax=ax  
            )


        plt.show()
        
        return fig
    
    

    
   
    
    
    def KEGG_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 5, 
                side = 'right', 
                color = 'orange', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
                ):
        
        sets = 'KEGG'
        column = '3rd'
        
        
        test_string = self.select_test(test, adj)

          
        tmp_in = pd.DataFrame(self.input_data['statistics'][sets])
        tmp_in = tmp_in[tmp_in[f'2nd_{test_string}'] <= p_val]
        tmp_in = tmp_in[tmp_in[f'3rd_{test_string}'] <= p_val]
        tmp_in[f'3rd_{test_string}'] = tmp_in[f'3rd_{test_string}'] + np.min(tmp_in[f'3rd_{test_string}'][tmp_in[f'3rd_{test_string}'] != 0])/2
        tmp_in['-log(p-val)'] = -np.log(tmp_in[f'3rd_{test_string}'])
        tmp_in = tmp_in.reset_index(drop = True)
        
        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(tmp_in['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(tmp_in['-log(p-val)'])
            
        else:
            x_max = np.max(tmp_in['n'])

        # queue
        
        tmp_qq = tmp_in[['2nd', '-log(p-val)']]
        
        tmp_qq['amount'] = tmp_qq['2nd'].map(tmp_qq['2nd'].value_counts())

        tmp_qq = tmp_qq.groupby('2nd', as_index=False).agg(
            amount=('amount', 'first'),  
            avg_log_pval=('-log(p-val)', 'mean') 
        )
        

        tmp_qq = tmp_qq.sort_values(
            by=['amount', 'avg_log_pval'], 
            ascending=[False, False]     
        ).reset_index(drop=True) 
               
        
        #######################################################################
        iter_set = list(tmp_qq['2nd'])

        
        hlist = []
        for it in iter_set:
            inx = [x for x in tmp_in.index if it in tmp_in['2nd'][x]] 
            tmp = tmp_in.loc[inx]
            if float(len(tmp[column])) > n:
                tn = n
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))
                else:
                    hlist.append(tn/2.1)

            else:
                tn = float(len(tmp[column]))
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))

                else:
                    hlist.append(tn/2.1)
                
            
        fig = plt.figure(figsize=(width, sum(hlist))) 

        gs = GridSpec(len(hlist), 1, height_ratios=hlist)  
        gs.update(hspace=len(hlist)/50)

        for l, i in enumerate(iter_set):
            inx = [x for x in tmp_in.index if i in tmp_in['2nd'][x]] 
            tmp = tmp_in.loc[inx]
            
            ax = fig.add_subplot(gs[l])
            
            show_axis = (l + 1 == len(iter_set))  
            
            self.bar_plot(
                data=tmp,
                n=n, 
                side=side, 
                color=color, 
                width=width, 
                bar_width=bar_width, 
                stat=stat, 
                sets='KEGG',
                column=column,
                x_max=x_max,
                show_axis=show_axis,
                title=i,  
                ax=ax  
            )


        plt.show()
        
        
        return fig
    
    
    def REACTOME_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 5, 
                side = 'right', 
                color = 'silver', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
            
                ):
        
        sets = 'REACTOME'
        column = 'pathway'
        
        
        test_string = self.select_test(test, adj)

        tmp_in = pd.DataFrame(self.input_data['statistics'][sets])
        tmp_in = tmp_in[tmp_in[f'top_level_pathway_{test_string}'] <= p_val]
        tmp_in = tmp_in[tmp_in[f'pathway_{test_string}'] <= p_val]
        tmp_in[f'pathway_{test_string}'] = tmp_in[f'pathway_{test_string}'] + np.min(tmp_in[f'pathway_{test_string}'][tmp_in[f'pathway_{test_string}'] != 0])/2
        tmp_in['-log(p-val)'] = -np.log(tmp_in[f'pathway_{test_string}'])
        tmp_in = tmp_in.reset_index(drop = True)
        
        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(tmp_in['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(tmp_in['-log(p-val)'])
            
        else:
            x_max = np.max(tmp_in['n'])

        # queue
        
        tmp_qq = tmp_in[['top_level_pathway', '-log(p-val)']]
        
        tmp_qq['amount'] = tmp_qq['top_level_pathway'].map(tmp_qq['top_level_pathway'].value_counts())

        tmp_qq = tmp_qq.groupby('top_level_pathway', as_index=False).agg(
            amount=('amount', 'first'),  
            avg_log_pval=('-log(p-val)', 'mean') 
        )
        

        tmp_qq = tmp_qq.sort_values(
            by=['amount', 'avg_log_pval'], 
            ascending=[False, False]     
        ).reset_index(drop=True) 
               
        
        #######################################################################
        iter_set = list(tmp_qq['top_level_pathway'])
        
        hlist = []
        for it in iter_set:
            inx = [x for x in tmp_in.index if it in tmp_in['top_level_pathway'][x]] 
            tmp = tmp_in.loc[inx]
            if float(len(tmp[column])) > n:
                tn = n
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))
                else:
                    hlist.append(tn/2.1)

            else:
                tn = float(len(tmp[column]))
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))

                else:
                    hlist.append(tn/2.1)
                
            
        fig = plt.figure(figsize=(width, sum(hlist))) 

        gs = GridSpec(len(hlist), 1, height_ratios=hlist)
        gs.update(hspace=len(hlist)/50)


        for l, i in enumerate(iter_set):
            inx = [x for x in tmp_in.index if i in tmp_in['top_level_pathway'][x]] 
            tmp = tmp_in.loc[inx]
            
            ax = fig.add_subplot(gs[l])
            
            show_axis = (l + 1 == len(iter_set))  
            
            self.bar_plot(
                data=tmp,
                n=n, 
                side=side, 
                color=color, 
                width=width, 
                bar_width=bar_width, 
                stat=stat, 
                sets='REACTOME',
                column=column,
                x_max=x_max,
                show_axis=show_axis,
                title=i,  
                ax=ax  
            )


        plt.show()
        
        
        return fig
    
    
    
    def SPECIFICITY_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 5, 
                side = 'right', 
                color = 'bisque', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
            
                ):
        
        sets = 'specificity'
        column = 'specificity'
        
        test_string = self.select_test(test, adj)

                  
        full_df = pd.DataFrame()

        for si in self.input_data['statistics'][sets].keys():

            tmp_in = pd.DataFrame(self.input_data['statistics'][sets][si])
        
            tmp_in = tmp_in[tmp_in[test_string] <= p_val]
            tmp_in[test_string] = tmp_in[test_string] + np.min(tmp_in[test_string][tmp_in[test_string] != 0])/2
            tmp_in['-log(p-val)'] = -np.log(tmp_in[test_string])
            tmp_in = tmp_in.reset_index(drop = True)
            tmp_in['set'] = si
            
            full_df = pd.concat([full_df, tmp_in])
        
        full_df = full_df.reset_index(drop = True)
        full_df['specificity'] = [x[0].upper() + x[1:] if isinstance(x, str) and len(x) > 0 else x for x in full_df['specificity']]

        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(full_df['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(full_df['-log(p-val)'])
            
        else:
            x_max = np.max(full_df['n'])

        # queue
        
        tmp_qq = full_df[['set', '-log(p-val)']]
        
        tmp_qq['amount'] = tmp_qq['set'].map(tmp_qq['set'].value_counts())

        tmp_qq = tmp_qq.groupby('set', as_index=False).agg(
            amount=('amount', 'first'),  
            avg_log_pval=('-log(p-val)', 'mean') 
        )
        

        tmp_qq = tmp_qq.sort_values(
            by=['amount', 'avg_log_pval'], 
            ascending=[False, False]     
        ).reset_index(drop=True) 
               
        
        #######################################################################
        iter_set = list(tmp_qq['set'])
        # colors = ['darkblue', 'blue', 'lightblue']
        
        hlist = []
        for it in iter_set:
            print(it)
            inx = [x for x in full_df.index if it in full_df['set'][x]] 
            tmp = full_df.loc[inx]
            if float(len(tmp[column])) > n:
                tn = n
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))
                else:
                    hlist.append(tn/2.1)

            else:
                tn = float(len(tmp[column]))
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))

                else:
                    hlist.append(tn/2.1)
                
            
        fig = plt.figure(figsize=(width, sum(hlist))) 

        gs = GridSpec(len(hlist), 1, height_ratios=hlist) 
        gs.update(hspace=len(hlist)/50)


        for l, i in enumerate(iter_set):
            inx = [x for x in full_df.index if i in full_df['set'][x]] 
            tmp = full_df.loc[inx]
            
            ax = fig.add_subplot(gs[l])
            
            show_axis = (l + 1 == len(iter_set))  
            
            self.bar_plot(
                data=tmp,
                n=n, 
                side=side, 
                color=color, 
                width=width, 
                bar_width=bar_width, 
                stat=stat, 
                sets='REACTOME',
                column=column,
                x_max=x_max,
                show_axis=show_axis,
                title=i,  
                ax=ax  
            )


        plt.show()
        
        
        return fig
    
    
    def DISEASES_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 5, 
                side = 'right', 
                color = 'thistle', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
            
                ):
        
        sets = 'DISEASES'
        column = 'disease'
        
        test_string = self.select_test(test, adj)

        
        tmp_in = pd.DataFrame(self.input_data['statistics'][sets])
    
        tmp_in = tmp_in[tmp_in[test_string] <= p_val]
        tmp_in[test_string] = tmp_in[test_string] + np.min(tmp_in[test_string][tmp_in[test_string] != 0])/2
        tmp_in['-log(p-val)'] = -np.log(tmp_in[test_string])
        tmp_in = tmp_in.reset_index(drop = True)
        tmp_in['disease'] = [x[0].upper() + x[1:] if isinstance(x, str) and len(x) > 0 else x for x in tmp_in['disease']]

        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(tmp_in['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(tmp_in['-log(p-val)'])
            
        else:
            x_max = np.max(tmp_in['n'])

        # queue
            
        fig = self.bar_plot(
            data=tmp_in,
            n=n, 
            side=side, 
            color=color, 
            width=width, 
            bar_width=bar_width, 
            stat=stat, 
            sets='DISEASES',
            column=column,
            x_max=x_max,
            show_axis=True,
            title='DISEASES',  
            ax=None  
        )


        plt.show()
        
        
        return fig
    
    
    def ViMIC_plot(self,   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 5, 
                side = 'right', 
                color = 'aquamarine', 
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val'
            
                ):
        
        sets = 'ViMIC'
        column = 'virus'
        
        test_string = self.select_test(test, adj)

        
        tmp_in = pd.DataFrame(self.input_data['statistics'][sets])
    
        tmp_in = tmp_in[tmp_in[test_string] <= p_val]
        tmp_in[test_string] = tmp_in[test_string] + np.min(tmp_in[test_string][tmp_in[test_string] != 0])/2
        tmp_in['-log(p-val)'] = -np.log(tmp_in[test_string])
        tmp_in = tmp_in.reset_index(drop = True)

        
        if stat.upper() == 'perc'.upper():
            x_max = np.max(tmp_in['pct'])
            
        elif stat.upper() == 'p_val'.upper():
            x_max = np.max(tmp_in['-log(p-val)'])
            
        else:
            x_max = np.max(tmp_in['n'])

        # queue
            
        fig = self.bar_plot(
            data=tmp_in,
            n=n, 
            side=side, 
            color=color, 
            width=width, 
            bar_width=bar_width, 
            stat=stat, 
            sets='ViMIC',
            column=column,
            x_max=x_max,
            show_axis=True,
            title='ViMIC',  
            ax=None  
        )


        plt.show()
        
        
        return fig
    
    
    def blod_markers_plot(self,   
                n = 10, 
                side = 'right', 
                color = 'red', 
                width = 10, 
                bar_width = 0.5            
                ):
        
               
        
        tmp_in = pd.DataFrame(self.input_data['enrichment']['HPA']['HPA_blood_markers'])
        
    
        x_max = max(np.log10(np.max(tmp_in['blood_concentration_IM[pg/L]'])), 
                       np.log10(np.max(tmp_in['blood_concentration_MS[pg/L]'])))
            
        
      
        iter_set = ['blood_concentration_IM[pg/L]', 'blood_concentration_MS[pg/L]']
        
        hlist = []
        for it in iter_set:
            print(it)
            
            tmp_len = len(tmp_in[it][tmp_in[it] == tmp_in[it]])
            
            if float(tmp_len) > n:
                tn = n
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))
                else:
                    hlist.append(tn/2.1)

            else:
                tn = float(tmp_len)
                if tn < 6:
                    if tn < 2:
                        hlist.append(tn/(2.5 + ((6-tn)/10)))
                    else:
                        hlist.append(tn/(2.12 + ((6-tn)/10)))

                else:
                    hlist.append(tn/2.1)
                
        
        fig = plt.figure(figsize=(width, sum(hlist))) 

        gs = GridSpec(len(hlist), 1, height_ratios=hlist)  
        
        
        gs.update(hspace=len(hlist)/30)

        for l, i in enumerate(iter_set):
           
            tmp = tmp_in[tmp_in[i] == tmp_in[i]]
            
            tmp[i] = np.log10(tmp[i])
            
            ax = fig.add_subplot(gs[l])
            
            show_axis = (l + 1 == len(iter_set))  
            
            self.bar_plot_blood(
                data=tmp,
                n=n, 
                side=side, 
                color=color, 
                width=width, 
                bar_width=bar_width, 
                stat=i, 
                sets='Blood markers',
                column='found_names',
                x_max=x_max,
                show_axis=show_axis,
                title=f'log({i})',  
                ax=ax  
            )


        plt.show()
        
        return fig
        
    
    
  
        
    
    def GOPa_network_create(self, data_set:str = 'GO-TERM', 
                            genes_inc:int = 10, 
                            gene_int:bool = True, 
                            genes_only:bool = True, 
                            min_con:int = 2, 
                            children_con:bool = False,
                            include_childrend:bool = True,
                            selected_parents:list = [],
                            selected_genes:list = []):
        

                           
       
        GOPa = pd.DataFrame(self.input_data['networks'][data_set])
        genes_list = list(set(GOPa['features']))
        
        
        if len(selected_genes) > 0:
            to_select_genes = []
            for p in selected_genes:
                if p in list(GOPa['features']):
                    to_select_genes.append(p)     
                else:
                    print('\Could not find {p} gene!')    
                    
            if len(to_select_genes) != 0:
                GOPa = GOPa[GOPa['features'].isin(to_select_genes)]
                genes_inc = max(genes_inc, len(to_select_genes))
            else:
                print('\nCould not use provided set of genes!')
       

        if data_set in ['GO-TERM', 'KEGG', 'REACTOME']:
            
         
            GOPa_drop = GOPa[['parent', 'children']].drop_duplicates()
            
            GOPa_drop = Counter(list(GOPa_drop['parent']))

            GOPa_drop = pd.DataFrame(GOPa_drop.items(), columns=['GOPa', 'n'])
            
            GOPa_drop = list(GOPa_drop['GOPa'][GOPa_drop['n'] >= min_con])
        
            GOPa = GOPa[GOPa['parent'].isin(GOPa_drop)]
            
            del GOPa_drop
            
           
            
            if genes_inc > 0:

                genes_list = GOPa['features']
                
                inter = None
                tmp_genes_list = []
                
                if gene_int == True:
                    inter = pd.DataFrame(self.input_data['statistics']['interactions'])
                    inter = inter[inter['A'] != inter['B']]
                    inter = inter[inter['A'].isin(genes_list)]
                    inter = inter[inter['B'].isin(genes_list)]
                    tmp_genes_list = list(set(list(inter['B']) + list(inter['A'])))
                    
                    if len(tmp_genes_list) > 0:
                        genes_list = tmp_genes_list

               
                genes_list = Counter(genes_list)
                
                genes_list = pd.DataFrame(genes_list.items(), columns=['features', 'n'])
                
                genes_list = genes_list.sort_values('n', ascending=False)
                
                
                gene_GOPa_p = GOPa[['parent', 'features']][GOPa['features'].isin(list(genes_list['features'][:genes_inc]))]
                gene_GOP_c = GOPa[['features','children']][GOPa['features'].isin(list(genes_list['features'][:genes_inc]))]
                genes_list = list(genes_list['features'][:genes_inc])
                
                if genes_only == True:
                    GOPa = GOPa[GOPa['features'].isin(genes_list)]

                gene_GOPa_p.columns = ['parent','children']
                
                gene_GOPa_p['color'] = 'gray'
                
                GOPa = pd.concat([GOPa[['parent','children', 'color']], gene_GOPa_p])
                
                if len(tmp_genes_list) > 0:
                    if isinstance(inter, pd.DataFrame):
                        inter = inter[inter['A'].isin(genes_list)]
                        inter = inter[inter['B'].isin(genes_list)]
                        inter = inter[['A', 'B']]
                        inter.columns = ['parent', 'children']
                        inter['color'] = 'red'
                        
                        GOPa = pd.concat([GOPa, inter])
                        
                
                
                if children_con == True:
                    
                    gene_GOP_c.columns = ['parent','children']
                    
                    gene_GOP_c['color'] = 'gray'
                    
                    GOPa = pd.concat([GOPa[['parent','children', 'color']], gene_GOP_c])
                    
                
                del gene_GOP_c, gene_GOPa_p
            
            
            
            gopa_list = list(GOPa['parent']) + list(GOPa['children'])

           
            gopa_list = Counter(gopa_list)
            
            gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'weight'])
            
            
            
            if len(selected_parents) > 0:
                to_select = []
                to_select_genes = []

                for p in selected_parents:
                    if p in list(GOPa['parent']):
                        to_select.append(p)   
                        t = list(GOPa['children'][GOPa['parent'] == p])
                        for i in t:
                            tg = [x for x in genes_list if x in list(GOPa['children'])]
                            if i in tg:
                                to_select_genes.append(i)   

                        
                    else:
                        print('\Could not find {p} parent term!')    
                        
                if len(to_select) != 0:
                    GOPa = GOPa[GOPa['parent'].isin(to_select + to_select_genes) & GOPa['children'].isin(list(GOPa['children'][GOPa['parent'].isin(to_select)]))]
                    gopa_list = gopa_list[gopa_list['GOPa'].isin(list(GOPa['parent']) + list(GOPa['children']))]

                else:
                    print('\nCould not use provided set of parent terms!')
              
                
            if include_childrend == False:
                GOPa = GOPa[GOPa['children'].isin(list(GOPa['parent']) + genes_list)]
                gopa_list = gopa_list[gopa_list['GOPa'].isin(list(GOPa['parent']) + genes_list)]
            
            
              
     
                
            G = nx.Graph() 
        
         
            for _, row in gopa_list.iterrows():
                node = row['GOPa']
                    
                if node in genes_list:
                    color = 'orange'
                    weight = np.log2(row['weight']*1000)

                elif node in list(GOPa['parent']):
                    color = 'cyan'
                    weight = np.log2(row['weight']*1000)*2
                else:
                    color = 'silver'
                    weight = np.log2(row['weight']*1000)


                G.add_node(node, size = weight, color = color)
                
            for _, row in GOPa.iterrows():
                source = row['parent']
                target = row['children']
                color = row['color']
                G.add_edge(source, target, color = color)
    
            
            return G
            
        else:
            
            print('\nWrong data set selected!')
            print('\nAvaiable data sets are included in:')

            for i in ['GO-TERM', 'KEGG', 'DISEASES', 'ViMIC', 'REACTOME']:
                print(f'\n{i}')
                
            

    def GI_network_create(self, data_set:str = 'GO-TERM', min_con:int = 2):
       
        inter = pd.DataFrame(self.input_data['statistics']['interactions'])
        inter = inter[['A', 'B', 'connection_type']]
        

        dict_meta = pd.DataFrame({
            'interactions': [
                ['gene -> gene'],
                ['protein -> protein'],
                ['gene -> protein'],
                ['protein -> gene'],
                ['gene -> gene', 'protein -> protein'],
                ['gene -> gene', 'gene -> protein'],
                ['gene -> gene', 'protein -> gene'],
                ['protein -> protein', 'gene -> protein'],
                ['protein -> protein', 'protein -> gene'],
                ['gene -> protein', 'protein -> gene'],
                ['gene -> gene', 'protein -> protein', 'gene -> protein'],
                ['gene -> gene', 'protein -> protein', 'protein -> gene'],
                ['gene -> gene', 'gene -> protein', 'protein -> gene'],
                ['protein -> protein', 'gene -> protein', 'protein -> gene'],
                ['gene -> gene', 'protein -> protein', 'gene -> protein', 'protein -> gene']
            ],
            'color': [
                '#f67089',
                '#f47832',
                '#ca9213',
                '#ad9d31',
                '#8eb041',
                '#4fb14f',
                '#33b07a',
                '#35ae99',
                '#36acae',
                '#38a9c5',
                '#3aa3ec',
                '#957cf4',
                '#cd79f4',
                '#f35fb5',
                '#f669b7'
            ]
        })
        
       
        
        
        genes_list = list(inter['A']) + list(inter['B'])
        
        genes_list = Counter(genes_list)
        
        genes_list = pd.DataFrame(genes_list.items(), columns=['features', 'n'])
        
        genes_list = genes_list.sort_values('n', ascending=False)
        
        genes_list = genes_list[genes_list['n'] >= min_con]
        
        inter = inter[inter['A'].isin(list(genes_list['features']))]
        inter = inter[inter['B'].isin(list(genes_list['features']))]
        
        

        inter = inter.groupby(['A', 'B']).agg({'connection_type': list}).reset_index()
        
        inter['color'] = 'black'
        
        for inx in inter.index:
            for inx2 in dict_meta.index:
                if set(inter['connection_type'][inx]) == set(dict_meta['interactions'][inx2]):
                    inter['color'][inx] = dict_meta['color'][inx2]
                    break

        G = nx.Graph() 
    
     
        for _, row in genes_list.iterrows():
            node = row['features']
            color = 'khaki'
            weight = np.log2(row['n']*500)
            G.add_node(node, size = weight, color = color)
            
        for _, row in inter.iterrows():
            source = row['A']
            target = row['B']
            color = row['color']
            G.add_edge(source, target, color = color)

        
        return G
            
       
        
       
         
    def AUTO_ML_network(self, 
                        genes_inc:int = 10, 
                        gene_int:bool = True, 
                        genes_only:bool = True, 
                        min_con:int = 2, 
                        children_con:bool = False, 
                        include_childrend:bool = False,
                        selected_parents:list = [],
                        selected_genes:list = []):
       
        
        full_genes = []
        genes_sets = []
        GOPa = pd.DataFrame()
        for s in ['GO-TERM', 'KEGG', 'REACTOME']:
            if s in self.input_data['networks'].keys():                
                genes_sets.append(set(self.input_data['networks'][s]['features']))
                full_genes += list(set(self.input_data['networks'][s]['features']))
                tmp = pd.DataFrame(self.input_data['networks'][s])
                tmp['set'] = s
                tmp['color'] = 'gray'
                GOPa = pd.concat([GOPa, tmp])
                
        
        common_elements = set.intersection(*genes_sets)
        
        del genes_sets
        
        
        inter = pd.DataFrame(self.input_data['statistics']['interactions'])
        inter = inter[inter['A'].isin(full_genes)]
        inter = inter[inter['B'].isin(full_genes)]
        
        if len(common_elements) > 0:
            inter = inter[inter['A'].isin(common_elements) | inter['B'].isin(common_elements)]
            

        selection_list = list(set(list(set(inter['B'])) + list(set(inter['A']))))
        
        if len(selected_genes) > 0:
            to_select_genes = []
            for p in selected_genes:
                if p in list(GOPa['features']):
                    to_select_genes.append(p)     
                else:
                    print('\Could not find {p} gene!')    
                    
            if len(to_select_genes) != 0:
                GOPa = GOPa[GOPa['features'].isin(to_select_genes)]
                genes_inc = max(genes_inc, len(to_select_genes))

            else:
                print('\nCould not use provided set of genes!')
                
        else:
           GOPa = GOPa[GOPa['features'].isin(selection_list)]

     
        GOPa_drop = GOPa[['parent', 'children']].drop_duplicates()
        
        GOPa_drop = Counter(list(GOPa_drop['parent']))

        GOPa_drop = pd.DataFrame(GOPa_drop.items(), columns=['GOPa', 'n'])
        
        GOPa_drop = list(GOPa_drop['GOPa'][GOPa_drop['n'] >= min_con])
    
        GOPa = GOPa[GOPa['parent'].isin(GOPa_drop)]
        
        del GOPa_drop

        
        if genes_inc > 0:

            genes_list = GOPa['features']
            
            inter = None
            tmp_genes_list = []


            if gene_int == True:
                inter = pd.DataFrame(self.input_data['statistics']['interactions'])
                inter = inter[inter['A'] != inter['B']]
                inter = inter[inter['A'].isin(genes_list)]
                inter = inter[inter['B'].isin(genes_list)]
                tmp_genes_list = list(set(list(inter['B']) + list(inter['A'])))
                
                if len(tmp_genes_list) > 0:
                    genes_list = tmp_genes_list

           
            genes_list = Counter(genes_list)
            
            genes_list = pd.DataFrame(genes_list.items(), columns=['features', 'n'])
            
            genes_list = genes_list.sort_values('n', ascending=False)
            
            
            gene_GOPa_p = GOPa[['parent', 'features', 'set', 'color']][GOPa['features'].isin(list(genes_list['features'][:genes_inc]))]
            gene_GOP_c = GOPa[['features','children', 'set', 'color']][GOPa['features'].isin(list(genes_list['features'][:genes_inc]))]
            genes_list = list(genes_list['features'][:genes_inc])
            
            if genes_only == True:
                GOPa = GOPa[GOPa['features'].isin(genes_list)]

            gene_GOPa_p.columns = ['parent','children', 'set', 'color']
            
            
            GOPa = pd.concat([GOPa[['parent','children', 'set', 'color']], gene_GOPa_p])
            
            if len(tmp_genes_list) > 0:
                if isinstance(inter, pd.DataFrame):
                    inter = inter[inter['A'].isin(genes_list)]
                    inter = inter[inter['B'].isin(genes_list)]
                    inter = inter[['A', 'B']]
                    inter.columns = ['parent', 'children']
                    inter['set'] = 'gene'
                    inter['color'] = 'red'

                    
                    GOPa = pd.concat([GOPa, inter])
                    
                    
            
            
            if children_con == True:
                
                gene_GOP_c.columns = ['parent','children', 'set', 'color']
                
                
                GOPa = pd.concat([GOPa[['parent','children', 'set', 'color']], gene_GOP_c])
                
            
            del gene_GOP_c, gene_GOPa_p
        
        
        
        
        gopa_list = list(GOPa['parent']) + list(GOPa['children'])

       
        gopa_list = Counter(gopa_list)
        
        gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'weight']).reset_index(drop = True)
      
        gopa_list['set'] = None
        
        for inx in gopa_list.index:
            if gopa_list['GOPa'][inx] in list(GOPa['parent']):
                gopa_list['set'][inx] = list(GOPa['set'][GOPa['parent'] == gopa_list['GOPa'][inx]])[0]
            elif gopa_list['GOPa'][inx] in list(GOPa['children']):
                gopa_list['set'][inx] = list(GOPa['set'][GOPa['children'] == gopa_list['GOPa'][inx]])[0]
            
            
            
        if len(selected_parents) > 0:
            to_select = []
            to_select_genes = []

            for p in selected_parents:
                if p in list(GOPa['parent']):
                    to_select.append(p)   
                    t = list(GOPa['children'][GOPa['parent'] == p])
                    for i in t:
                        tg = [x for x in genes_list if x in list(GOPa['children'])]
                        if i in tg:
                            to_select_genes.append(i)   

                    
                else:
                    print('\Could not find {p} parent term!')    
                    
            if len(to_select) != 0:
                GOPa = GOPa[GOPa['parent'].isin(to_select + to_select_genes) & GOPa['children'].isin(list(GOPa['children'][GOPa['parent'].isin(to_select)]))]
                gopa_list = gopa_list[gopa_list['GOPa'].isin(list(GOPa['parent']) + list(GOPa['children']))]

            else:
                print('\nCould not use provided set of parent terms!')
          
          
            
          
        if include_childrend == False:
            GOPa = GOPa[GOPa['children'].isin(list(GOPa['parent']) + genes_list)]
            gopa_list = gopa_list[gopa_list['GOPa'].isin(list(GOPa['parent']) + genes_list)]
        
        
        
            
        G = nx.Graph() 
    
        for _, row in gopa_list.iterrows():
            node = row['GOPa']
            
            color = 'black'
            
            if node in genes_list:
                color = 'orange'
                weight = np.log2(row['weight']*1000)

            elif node in list(GOPa['parent']):
                color = 'cyan'
                weight = np.log2(row['weight']*1000)*2


            else:
                if row['set'] == 'GO-TERM':
                    color = 'bisque'
                    weight = np.log2(row['weight']*1000)

                elif row['set'] == 'KEGG':
                    color = 'mistyrose'
                    weight = np.log2(row['weight']*1000)

                elif row['set'] == 'REACTOME':
                    color = 'darkkhaki'
                    weight = np.log2(row['weight']*1000)

              

           
            G.add_node(node, size = weight, color = color)
            
        for _, row in GOPa.iterrows():
            source = row['parent']
            target = row['children']
            color = row['color']
            G.add_edge(source, target, color = color)

        
        return G
        
        
    
    
    def gene_scatter(self, 
                     colors = 'viridis', 
                     species = 'human', 
                     hclust = 'complete', 
                     img_width = None, 
                     img_high = None, 
                     label_size = None, 
                     x_lab = 'Genes', 
                     legend_lab = 'log(TPM + 1)'):
    
            
        """
        This function creates a graph in the format of a scatter plot for expression data prepared in data frame format.
        
        Args:
            data (data frame) - data frame of genes/protein expression where on row are the gene/protein names and on column grouping variable (tissue / cell / ect. names)
            color (str) - palette color available for matplotlib in python eg. viridis
            species (str) - species for upper() or lower() letter for gene/protein name depending on 
            hclust (str) - type of data clustering of input expression data eg. complete or None if  no clustering
            img_width (float) - width of the image or None for auto-adjusting
            img_high (float) - high of the image or None for auto-adjusting
            label_size (float) - labels size of the image or None for auto-adjusting
            x_lab (str) - tex for x axis label
            legend_lab (str) - description for legend label
           
           
        Returns:
            graph: Scatter plot of expression data
        """
        
        input_data = self.input_data['enrichment']['RNA-SEQ']
        
        return_dict = {}
        
        for i in input_data.keys():
            data = pd.DataFrame(input_data[i])
            data.index = data['tissue']
            data.pop('tissue')
        
        
      
            scatter_df = data        
    
            if img_width == None:
                img_width = len(scatter_df.columns)*1.2
            
            if img_high == None:
                img_high = len(scatter_df.index)*0.9
                
            
            if label_size == None:
                label_size = np.log(len(scatter_df.index)  *  len(scatter_df.index))*2.5
                
                if label_size < 7:
                    label_size = 7
            
            cm = 1/2.54
            
            if len(scatter_df) > 1:
                
               
                
            
                Z = linkage(scatter_df, method=hclust)
        
        
                # Get the order of features based on the dendrogram
                order_of_features = dendrogram(Z, no_plot=True)['leaves']
        
                indexes_sort = list(scatter_df.index)
                sorted_list_rows = []
                for n in order_of_features:
                    sorted_list_rows.append(indexes_sort[n])
                    
                
                
                scatter_df = scatter_df.transpose()
            
                Z = linkage(scatter_df, method=hclust)
        
                # Get the order of features based on the dendrogram
                order_of_features = dendrogram(Z, no_plot=True)['leaves']
        
                indexes_sort = list(scatter_df.index)
                sorted_list_columns = []
                for n in order_of_features:
                    sorted_list_columns.append(indexes_sort[n])
                            
                       
                scatter_df = scatter_df.transpose()
                
                scatter_df = scatter_df.loc[sorted_list_rows, sorted_list_columns]
                 
            scatter_df = np.log(scatter_df + 1)
            scatter_df[scatter_df <= np.mean(scatter_df.quantile(0.10))] = np.mean(np.mean(scatter_df, axis=1))/10
        
            if species.lower() == 'human':
                scatter_df.index = [x.upper() for x in scatter_df.index ]
            else:
                scatter_df.index  = [x.title() for x in scatter_df.index ]
                
            scatter_df.insert(0, '  ', 0)
    
            # Add a column of zeros at the end
            scatter_df[' '] = 0
                 
            fig, ax = plt.subplots(figsize=(img_width*cm,img_high*cm))
        
            plt.scatter(x = [*range(0, len(scatter_df.columns), 1)], y = [' '] * len(scatter_df.columns),s=0, cmap=colors,  edgecolors=None)
        
            
    
        
            for index, row in enumerate(scatter_df.index):
                x = [*range(0, len(np.array(scatter_df.loc[row,])), 1)]
                y = [row] * len(x)
                s = np.array(scatter_df.loc[row,])
                plt.scatter(x,y,s=np.log(s+1)*70, c=s, cmap=colors,  edgecolors='black', vmin = np.array(scatter_df).min() ,vmax = np.array(scatter_df).max(), linewidth=0.00001)
                sm = plt.cm.ScalarMappable(cmap=colors)
                sm.set_clim(vmin = np.array(scatter_df).min() ,vmax = np.array(scatter_df).max())
                plt.xticks(x, scatter_df.columns)
                plt.ylabel(str(x_lab), fontsize=label_size)
                
                
            plt.scatter(x = [*range(0, len(scatter_df.columns), 1)], y = [''] * len(scatter_df.columns),s=0, cmap=colors,  edgecolors=None)
    
    
            
        
            plt.xticks(rotation = 80) 
            plt.tight_layout()
            plt.margins(0.005)
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)
     
        
        
            len_bar = ax.get_position().height/5
            if len(scatter_df) < 15:
                len_bar = 0.65
                
                cbar = plt.colorbar(sm)
                cbar.ax.set_ylabel(str(legend_lab), fontsize=label_size*0.9)
                cbar.ax.yaxis.set_ticks_position('right')
                cbar.ax.set_position([ax.get_position().x1 + 0.05, (ax.get_position().y0 + ax.get_position().y1)/1.9 , ax.get_position().width/0.05, len_bar])
                cbar.ax.yaxis.set_label_position('right')
                cbar.ax.yaxis.set_tick_params(labelsize=label_size*0.8)
                cbar.outline.set_edgecolor('none')
            else:
                cbar = plt.colorbar(sm)
                cbar.ax.set_ylabel(str(legend_lab), fontsize=label_size*0.9)
                cbar.ax.yaxis.set_ticks_position('right')
                cbar.ax.set_position([ax.get_position().x1 + 0.05, (ax.get_position().y0 + ax.get_position().y1)/1.45 , ax.get_position().width/0.05, len_bar])
                cbar.ax.yaxis.set_label_position('right')
                cbar.ax.yaxis.set_tick_params(labelsize=label_size*0.8)
                cbar.outline.set_edgecolor('none')
        
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)     
            ax.xaxis.set_tick_params(length=0,labelbottom=True)
            ax.yaxis.set_tick_params(length=0,labelbottom=True)
            ax.grid(False)
        
        
            return_dict[i] = fig
            
        return return_dict
        
        





class DSA():
    
    """
    This function conducts full GOPa analysis on two gene / protein lists [gene_up, gene_down] including search_genes, gopa_analysis, gopa_interaction_analysis, and gopa_specificity_analysis.
    This function is useful when you have to compare two sides genes / proteins deregulations (upregulated vs. downregulated) and choose the most accurate results adjusted to the provided gene lists based on min_fc (minimal fold change) between GOPa projects.

    Args:
        gene_up (list)- list of genes eg. ['KIT', 'EDNRB', 'PAX3'] 
        gene_down (list)- list of genes eg. ['KIT', 'EDNRB', 'PAX3'] 
        GOPa_metadata (dict) - metadata from load_GOPa_meta function 
        species (str or None) - ['human' / 'mouse' / 'both' / None] 
        min_fc (float) - minimal value of fold change the normalized by number of genes in analysis results of GOPa between GOPa projects obtained from the provided lists of genes [gene_up / gene_down]
        p_val (float) - value of minimal p_val for statistical test
        test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
        adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
       
        If choose 'human' or 'mouse' you will obtain information about this species' genes. 
        If choose 'both' you will obtain information for genes that are available mutually for both species. 
        If choose None you will obtain information for all genes available in the metadata.       

    Returns:
        list of dict: A list of two dictionaries with analyzed GOPa projects obtained on the provided list of genes [gene_up, gene_down] corrected on specific occurrences in GOPa results with min_fc. The first dictionary is related to gene_up results and the second to gene_down results 
       
    """
    
    def __init__(self, set_1:dict, set_2:dict):
        
        
        self.set_1 = set_1
        self.set_2 = set_2
        
        if self.set_1['enrichment']['species']['species_genes'] != self.set_2['enrichment']['species']['species_genes']:
            raise ValueError("The 'self.species_genes' attribute used for enrichment analysis differed between set_1 and set_2.")


        if self.set_1['enrichment']['species']['species_study'] != self.set_2['enrichment']['species']['species_study']:
            raise ValueError("The 'self.species_study' attribute used for enrichment analysis differed between set_1 and set_2.")
           
            
        if self.set_1['statistics']['setup']['network_stat'] != self.set_2['statistics']['setup']['network_stat']:
            raise ValueError("The 'self.network_stat' attribute used for analysis differed between set_1 and set_2.")


        if self.set_1['statistics']['setup']['go_grade'] != self.set_2['statistics']['setup']['go_grade']:
            raise ValueError("The 'self.go_grade' attribute used for analysis differed between set_1 and set_2.")


        if self.set_1['statistics']['setup']['interaction_strength'] != self.set_2['statistics']['setup']['interaction_strength']:
            raise ValueError("The 'self.interaction_strength' attribute used for analysis differed between set_1 and set_2.")


        if self.set_1['statistics']['setup']['interaction_source'] != self.set_2['statistics']['setup']['interaction_source']:
            raise ValueError("The 'self.interaction_source' attribute used for analysis differed between set_1 and set_2.")


        self.min_fc = 0.75 
        self.s1_genes = len(self.set_1['enrichment']['gene_info']['sid'])
        self.s2_genes = len(self.set_2['enrichment']['gene_info']['sid'])
        self.GO = None
        self.KEGG = None
        self.REACTOME = None
        self.specificity = None
        self.GI = None
        self.networks = None


    @property
    def get_GO_diff(self):
        return self.GO.to_dict(orient = 'list')
    
    
    @property
    def get_KEGG_diff(self):
        
        return self.KEGG.to_dict(orient = 'list')
    
    
    @property
    def get_REACTOME_diff(self):
        
        return self.REACTOME.to_dict(orient = 'list')
    
    
    @property
    def get_specificity_diff(self):
        
        return self.specificity.to_dict(orient = 'list')
    
    
    @property
    def get_GI_diff(self):
        
        return self.GI.to_dict(orient = 'list')
    
    
    @property
    def get_networks_diff(self):
        
        return self.networks



            
        
    def GO_diff(self):
        
        # parent check
        
        parent_columns = 'parent_name'
        parent_n = 'parent_n'
        
        children_columns = 'child_name'
        children_n = 'child_n'

        sets = 'GO-TERM'

        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[parent_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[parent_n][s1_tmp[parent_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[parent_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[parent_n][s2_tmp[parent_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2:  
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else:  
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp['type'] = 'parent'
        
        
        
        # children check
        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[children_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[children_n][s1_tmp[children_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[children_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[children_n][s2_tmp[children_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2:  
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else:  
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp1 = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp1['type'] = 'children'
        
        full_values = pd.concat([tmp, tmp1])
        
        
        self.GO = full_values


 
    def KEGG_diff(self):
        
        # parent check
        
        parent_columns = '2nd'
        parent_n = '2nd_n'
        
        children_columns = '3rd'
        children_n = '3rd_n'

        sets = 'KEGG'

        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[parent_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[parent_n][s1_tmp[parent_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[parent_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[parent_n][s2_tmp[parent_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2:  
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else:  
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp['type'] = 'parent'
        
        
        
        # children check
        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[children_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[children_n][s1_tmp[children_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[children_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[children_n][s2_tmp[children_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2:  
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else:  
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp1 = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp1['type'] = 'children'
        
        full_values = pd.concat([tmp, tmp1])
        
        
        self.KEGG = full_values
        
        
    def REACTOME_diff(self):
        
        # parent check
        
        parent_columns = 'top_level_pathway'
        parent_n = 'top_level_pathway_n'
        
        children_columns = 'pathway'
        children_n = 'pathway_n'

        sets = 'REACTOME'

        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[parent_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[parent_n][s1_tmp[parent_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[parent_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[parent_n][s2_tmp[parent_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2: 
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else: 
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp['type'] = 'parent'
        
        
        
        # children check
        
        s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets])
        s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets])
        
        

        term = []
        norm_n = []
        
        for s1 in set(s1_tmp[children_columns]):
               
            term.append(s1)
            norm_n.append(float(np.mean(s1_tmp[children_n][s1_tmp[children_columns] == s1]))/self.s1_genes)
        
        s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        
        term = []
        norm_n = []
        
        for s2 in set(s2_tmp[children_columns]):
               
            term.append(s2)
            norm_n.append(float(np.mean(s2_tmp[children_n][s2_tmp[children_columns] == s2]))/self.s2_genes)
        
        s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})

        
        term = []
        FC = []
        dec = []
        
        terms_s1 = set(s1['term'])
        terms_s2 = set(s2['term'])
        
        norm_s1 = dict(zip(s1['term'], s1['norm_n']))
        norm_s2 = dict(zip(s2['term'], s2['norm_n']))
        
        all_terms = terms_s1.union(terms_s2)
        
        min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
        
        for g in all_terms:
            if g in terms_s1 and g in terms_s2:  
                fc_value = norm_s1[g] / norm_s2[g]
                if fc_value >= self.min_fc:
                    decision = 's1'
                elif fc_value <= 1 / self.min_fc:  
                    decision = 's2'
                else:
                    decision = 'equal'
            elif g in terms_s1:  
                fc_value = norm_s1[g] / min_value
                decision = 's1'
            else:  
                fc_value = min_value / norm_s2[g]
                decision = 's2'
        
            term.append(g)
            FC.append(fc_value)
            dec.append(decision)

        
        
        tmp1 = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
        tmp1['type'] = 'children'
        
        full_values = pd.concat([tmp, tmp1])
        
        
        self.REACTOME = full_values
        
        
        
    def spec_diff(self):
        
        
        parent_columns = 'specificity'
        parent_n = 'n'
      

        sets = 'specificity'
        
        key_list = self.set_1['statistics'][sets].keys()
        
        full_df = pd.DataFrame()
        
        for k in key_list:
            s1_tmp =  pd.DataFrame(self.set_1['statistics'][sets][k])
            s2_tmp =  pd.DataFrame(self.set_2['statistics'][sets][k])
            
    
            term = []
            norm_n = []
            
            for s1 in set(s1_tmp[parent_columns]):
                   
                term.append(s1)
                norm_n.append(float(np.mean(s1_tmp[parent_n][s1_tmp[parent_columns] == s1]))/self.s1_genes)
            
            s1 = pd.DataFrame({'term':term, 'norm_n': norm_n})
    
            
            
            term = []
            norm_n = []
            
            for s2 in set(s2_tmp[parent_columns]):
                   
                term.append(s2)
                norm_n.append(float(np.mean(s2_tmp[parent_n][s2_tmp[parent_columns] == s2]))/self.s2_genes)
            
            s2 = pd.DataFrame({'term':term, 'norm_n': norm_n})
    
            
            term = []
            FC = []
            dec = []
            
            terms_s1 = set(s1['term'])
            terms_s2 = set(s2['term'])
            
            norm_s1 = dict(zip(s1['term'], s1['norm_n']))
            norm_s2 = dict(zip(s2['term'], s2['norm_n']))
            
            all_terms = terms_s1.union(terms_s2)
            
            min_value = min(list(norm_s1.values()) + list(norm_s2.values())) / 10
            
            for g in all_terms:
                if g in terms_s1 and g in terms_s2: 
                    fc_value = norm_s1[g] / norm_s2[g]
                    if fc_value >= self.min_fc:
                        decision = 's1'
                    elif fc_value <= 1 / self.min_fc:  
                        decision = 's2'
                    else:
                        decision = 'equal'
                elif g in terms_s1:  
                    fc_value = norm_s1[g] / min_value
                    decision = 's1'
                else: 
                    fc_value = min_value / norm_s2[g]
                    decision = 's2'
            
                term.append(g)
                FC.append(fc_value)
                dec.append(decision)
    
            
            
            tmp = pd.DataFrame({'term':term, 'FC':FC, 'regulation':dec})
            tmp['set'] = k
            
            full_df = pd.concat([full_df, tmp])
        
        
        
        
        self.specificity = full_df
        
        
        
        
    def gi_diff(self):
        
        enr = Enrichment()
        
        enr.species_genes = self.set_1['enrichment']['species']['species_genes'] 

        enr.species_study = self.set_1['enrichment']['species']['species_study'] 

                
        enr.genome = pd.concat([pd.DataFrame(self.set_1['enrichment']['gene_info']), 
                                pd.DataFrame(self.set_2['enrichment']['gene_info'])])      
        
        enr.enriche_IntAct()
        
        enr.enriche_STRING()
        
        
        res = enr.get_results
        
        del enr
        
        ans = Analysis(res)
        
               
        ans.network_stat = self.set_1['statistics']['setup']['network_stat'] 

        ans.go_grade = self.set_1['statistics']['setup']['go_grade'] 
        
        ans.interaction_strength = self.set_1['statistics']['setup']['interaction_strength'] 
        
        ans.interaction_source = self.set_1['statistics']['setup']['interaction_source'] 

        ans.gene_interaction()
        
        full_ans = pd.DataFrame(ans.get_features_interactions_statistics)
        
        
        s1_ans = pd.DataFrame(self.set_1['statistics']['interactions'])
        s2_ans = pd.DataFrame(self.set_2['statistics']['interactions'])
        
        full_ans['set'] = 'inter'
        
        full_ans['set'][full_ans['A'].isin(list(s1_ans['A'])) & full_ans['B'].isin(list(s1_ans['B']))] = 's1'
        full_ans['set'][full_ans['A'].isin(list(s2_ans['A'])) & full_ans['B'].isin(list(s2_ans['B']))] = 's2'
        
        self.GI = full_ans



        
    def network_diff(self):
        
        networks = {}
        
        enr = Enrichment()
        
        enr.species_genes = self.set_1['enrichment']['species']['species_genes'] 

        enr.species_study = self.set_1['enrichment']['species']['species_study'] 

                
        enr.genome = pd.concat([pd.DataFrame(self.set_1['enrichment']['gene_info']), 
                                pd.DataFrame(self.set_2['enrichment']['gene_info'])])      
        
        enr.enriche_GOTERM()
        
        enr.enriche_KEGG()
        
        enr.enriche_REACTOME()

        res = enr.get_results
        
        del enr
        
        ans = Analysis(res)
        
        ans.network_stat = self.set_1['statistics']['setup']['network_stat'] 

        ans.go_grade = self.set_1['statistics']['setup']['go_grade'] 
        
        ans.interaction_strength = self.set_1['statistics']['setup']['interaction_strength'] 
        
        ans.interaction_source = self.set_1['statistics']['setup']['interaction_source'] 
        
        ans.REACTOME_overrepresentation()
        ans.REACTOME_network()

        ans.KEGG_overrepresentation()
        ans.KEGG_network()
        
        ans.GO_overrepresentation()
        ans.GO_network()
        
        
        kn = pd.DataFrame(ans.get_KEGG_network)
        
        s1_ans = pd.DataFrame(self.set_1['networks']['KEGG'])
        s2_ans = pd.DataFrame(self.set_2['networks']['KEGG'])
        
        kn['set'] = 'inter'
        
        kn['set'][kn['parent'].isin(list(s1_ans['parent'])) & kn['children'].isin(list(s1_ans['children']))] = 's1'
        kn['set'][kn['parent'].isin(list(s2_ans['parent'])) & kn['children'].isin(list(s2_ans['children']))] = 's2'
        
        
        networks['KEGG'] = kn.to_dict(orient = 'list')
        
      
            
        kn = pd.DataFrame(ans.get_REACTOME_network)
        
        s1_ans = pd.DataFrame(self.set_1['networks']['REACTOME'])
        s2_ans = pd.DataFrame(self.set_2['networks']['REACTOME'])
        
        kn['set'] = 'inter'
        
        kn['set'][kn['parent'].isin(list(s1_ans['parent'])) & kn['children'].isin(list(s1_ans['children']))] = 's1'
        kn['set'][kn['parent'].isin(list(s2_ans['parent'])) & kn['children'].isin(list(s2_ans['children']))] = 's2'
        
        
        networks['REACTOME'] = kn.to_dict(orient = 'list')
        
        
                
        kn = pd.DataFrame(ans.get_GO_network)
        
        s1_ans = pd.DataFrame(self.set_1['networks']['GO-TERM'])
        s2_ans = pd.DataFrame(self.set_2['networks']['GO-TERM'])
        
        kn['set'] = 'inter'
        
        kn['set'][kn['parent'].isin(list(s1_ans['parent'])) & kn['children'].isin(list(s1_ans['children']))] = 's1'
        kn['set'][kn['parent'].isin(list(s2_ans['parent'])) & kn['children'].isin(list(s2_ans['children']))] = 's2'
        
        
        networks['GO-TERM'] = kn.to_dict(orient = 'list')
        
       
        
        self.networks = networks
        
    
    @property
    def get_results(self):
        
        results = {}
        
        results['set_1'] = self.set_1
        results['set_2'] = self.set_2

        
        results['GO-TERM'] = self.get_GO_diff
        
        
        results['KEGG'] = self.get_KEGG_diff
        
        
        results['REACTOME'] = self.get_REACTOME_diff


        results['specificity'] = self.get_specificity_diff
            
        
        results['GI'] = self.get_GI_diff
            
           
        
        results['GI'] = self.get_networks_diff
            
            
        
        return results


        
        
        
        

        


