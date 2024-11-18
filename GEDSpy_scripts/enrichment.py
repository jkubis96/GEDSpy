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
    def get_RNA_SEQ(self):
        
        if self.RNA_SEQ != None:
            
            return self.RNA_SEQ
        
        else:
            raise ValueError('\nLack of enrichment of with RNA-SEQ data...')



    @property
    def get_HPA(self):
        
        if self.HPA != None:
            
            return self.HPA
        
        else:
            raise ValueError('\nLack of enrichment of with HPA data...')



    @property
    def get_STRING(self):
        
        if self.STRING != None:
            
            return self.STRING
        
        else:
            raise ValueError('\nLack of enrichment of with STRING data...')



    @property
    def get_GO_TERM(self):
        
        if self.GO_TERM != None:
            
            return self.GO_TERM
        
        else:
            raise ValueError('\nLack of enrichment of with GO-TERM data...')



    @property
    def get_REACTOME(self):
        
        if self.REACTOME != None:
            
            return self.REACTOME
        
        else:
            raise ValueError('\nLack of enrichment of with REACTOME data...')

     
  
    @property
    def get_CellCon(self):
        
        if self.CellCon != None:
            
            return self.CellCon
        
        else:
            raise ValueError('\nLack of enrichment of with CellCon data...')

     
    @property
    def get_IntAct(self):
        
        if self.IntAct != None:
            
            return self.IntAct
        
        else:
            raise ValueError('\nLack of enrichment of with IntAct data...')

     
        
    @property
    def get_KEGG(self):
        
        if self.KEGG != None:
            
            return self.KEGG
           
        else:
            raise ValueError('\nLack of enrichment of with KEGG data...')

        
    @property
    def get_DISEASES(self):
        
        if self.Diseases != None:
            
            return self.Diseases
        
        else:
            raise ValueError('\nLack of enrichment of with Diseases data...')

        
        
    @property
    def get_ViMIC(self):
        
        if self.ViMIC != None:
            
            return self.ViMIC
        
        else:
            raise ValueError('\nLack of enrichment of with ViMIC data...')

            
            
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
                
                
                try:
                    df['name'] = [x.capitalize() for x in  df['name']]
                except:
                    pass
                
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
            query = f"SELECT * FROM GO_go_names WHERE GO_id IN ({ids});"
            
            df = pd.read_sql_query(query, conn).applymap(self.deserialize_data)
            
            final_dict['go_names'] = df.to_dict(orient = 'list')
            
            query = f"SELECT * FROM GO_hierarchy WHERE GO_id IN ({ids});"
    
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
            
            df['disease'] = [x.capitalize() for x in  df['disease']]
    
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
        
            conn = sqlite3.connect(os.path.join(self.path_in_inside, 'GEDS_db.db'))
            
            ids = [int(x) for x in self.genome['id_cell_int'] if x == x]
            
            
             
            species = self.species_study 
            species = ', '.join(map(lambda x: f"'{x}'", species))
            
            
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
    
            
            conn.close()
    
            self.CellCon = df.to_dict(orient = 'list')
            
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
       
        
        
        
        
        
        

# class enrich_data(_tools):
    
 
    
#     #######################################ANALYSIS########################################
    
#     #Analysis fucntions
    
    
    
#     def load_GOPa_meta(path_in_use = _path_in_inside):
        
#         print('\n')
#         print('Metadata loading...')
    
#         with open(path_in_use + '/gene_dictionary_jbio_annotated.json', 'r') as json_file:
#             GOPa_metadata = (json.load(json_file))
        
            
#         return GOPa_metadata
        
    
#     # beginign variables
#     input_metadata = load_GOPa_meta(path_in_use = _path_in_inside)
    
    
    
#     epd_genes = pd.read_csv('EPD/epd_promoters_202405211235.csv', sep =';')
#     epd_genes = epd_genes.drop_duplicates()
#     epd_genes = epd_genes.reset_index(drop = True)
#     epd_genes['gene_name2'] = [x.upper() for x in epd_genes['gene_name']]
#     gene_list = list(epd_genes['gene_name'])

    
    
#     epd_genes = pd.merge(epd_genes, names_dict, how = 'left', left_on = 'gene_name2', right_on = 'fetures_names')

#     epd_genes = epd_genes.drop(['fetures_names', 'gene_name2'], axis = 1)
    
#     epd_genes = epd_genes.rename(columns={'dict_inx':'sc_id'})
    
#     epd_genes['is_nc'][epd_genes['is_nc'] == 1] = True
#     epd_genes['is_nc'][epd_genes['is_nc'] == 0] = False

    
#     epd_genes.to_csv('epd_with_id.csv', index=False)
    

 

    
    



    
    
    
#     def search_genes(gene_list:list, GOPa_metadata, species=None):
        
#         """
#         This function checks if the all genes provided in the list are in the data.
    
#         Args:
#            gene_list (list)- list of genes eg. ['KIT', 'EDNRB', 'PAX3'] 
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function 
#            species (str or None) - ['human' / 'mouse' / 'both' / None] 
           
#            If choose 'human' or 'mouse' you will obtain information about this species' genes. 
#            If choose 'both' you will obtain information for genes that are available mutually for both species. 
#            If choose None you will obtain information for all genes available in the metadata.       
    
#         Returns:
#            dict: A dictionary of all information available in metadata for genes provide in 'gene_list:'
#        """
       
#         print('\n')
#         print('Searching genes from list in the GEDSpy data...')
    
        
#         try:
            
#             gene_list = [re.sub('.chr.*','', x) for x in gene_list]
            
#             gene_dictionary = pd.DataFrame(copy.deepcopy(input_metadata))
            
#             gene_list = list(set([x.upper() for x in gene_list]))
            
#             gene_dictionary.reset_index(drop = True)     
            
#             ids = find_fetures_id(gene_dictionary, fetures_list = gene_list)
            
#             names_dict = pd.DataFrame({'fetures_names': ids['found_genes'], 'dict_inx':ids['found_ids']})
            
#             species_ids = spec_dic(gene_dictionary, species = 'both1')
            
#             names_dict = names_dict[names_dict['dict_inx'].isin(species_ids)]
#             names_dict = names_dict.groupby('dict_inx')[['fetures_names']].agg(list).reset_index()
#             names_dict['fetures_names'] = ['; '.join(x) for x in names_dict['fetures_names']]

#             gene_dictionary = gene_dictionary[gene_dictionary['sid'].isin(names_dict['dict_inx'])]
            
#             gene_dictionary = pd.merge(gene_dictionary, names_dict, how = 'left', left_on = 'sid', right_on = 'dict_inx')

#             gene_dictionary = gene_dictionary.drop(['dict_inx'], axis = 1)

            
            

            
            
           
            
            
    
            
            
#             GOPa_results = GOPa_metadata2
            
#             GOPa_results['gene_dictionary'] = genes_df.to_dict(orient = 'list')
            
#             GOPa = pd.DataFrame(GOPa_results['GOPa'])
            
#             GOPa = GOPa[GOPa['dictionary_id'].isin(list(genes_df['dictionary_id']))]
            
#             GOPa_results['GOPa'] = GOPa.to_dict(orient = 'list')
        
#             #
#             GOPa_interactions = pd.DataFrame(GOPa_results['GOPa_interactions'])
             
#             GOPa_interactions = GOPa_interactions[GOPa_interactions['GO_id'].isin(list(GOPa['relation_id'][GOPa['source'] == 'GO-TERM']))]
             
#             GOPa_results['GOPa_interactions'] = GOPa_interactions.to_dict(orient = 'list')
            
#             del GOPa_interactions
        
#             #
            
#             GOPa_genes = pd.DataFrame(GOPa_results['GOPa_gene_interaction'])
            
#             GOPa_genes = GOPa_genes[GOPa_genes['id_1'].isin(list(genes_df['dictionary_id']))]
            
#             GOPa_genes = GOPa_genes[GOPa_genes['id_2'].isin(list(genes_df['dictionary_id']))]
        
#             GOPa_results['GOPa_gene_interaction'] = GOPa_genes.to_dict(orient = 'list')
            
#             del GOPa_genes
            
            
#             #
            
#             GOPa_specificity = GOPa_results['GOPa_specificity']
            
#             SEQ = pd.DataFrame(GOPa_specificity['SEQ'])
            
#             SEQ = SEQ[SEQ['dictionary_id'].isin(list(genes_df['dictionary_id']))]
            
#             GOPa_specificity['SEQ'] = SEQ.to_dict(orient = 'list')
            
#             del SEQ
            
#             location = pd.DataFrame(GOPa_specificity['location'])
            
#             location = location[location['dictionary_id'].isin(list(genes_df['dictionary_id']))]
            
#             GOPa_specificity['location'] = location.to_dict(orient = 'list')
            
#             del location
            
            
#             blood_levels = pd.DataFrame(GOPa_specificity['blood_levels'])
            
#             blood_levels = blood_levels[blood_levels['dictionary_id'].isin(list(genes_df['dictionary_id']))]
            
#             GOPa_specificity['blood_levels'] = blood_levels.to_dict(orient = 'list')
            
#             del blood_levels
            
            
            
#             GOPa_results['GOPa_specificity'] = GOPa_specificity
            
#             if len(not_found) > 0:
#                 print('\n')
#                 print(str(len(not_found)) + ' out of ' + str(len(gene_list)) + ' values were not found')
                
#             return GOPa_results, not_found
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
#     def gopa_analysis(GOPa_data, GOPa_metadata):
        
#         """
#         This function conducts statistical / overrepresentation analysis on raw GOPa_data.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from search_genes function
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function
         
    
#         Returns:
#            dict: A GOPa_data results with appropriate statistics
#        """
       
#         print('\n')
#         print('Overrepresentation analysis of terms in GOPa data...')
        
#         try:
           
#             GOPa_metadata2 = copy.deepcopy(GOPa_metadata)
#             GOPa_res = pd.DataFrame(GOPa_data['GOPa'])
            
            
#             input_genes = int(len(set(GOPa_res['gene_name'])))
        
#             total_genes = int(len(set(GOPa_metadata2['GOPa']['gene_name'])))
            
#             GOPa_metadata2 = pd.DataFrame(GOPa_metadata2['GOPa'])
            
#             GOPa_ocr = GOPa_metadata2.groupby('GOPa').agg({'dictionary_id': list}).reset_index()
#             GOPa_ocr['ocr_n'] = [int(len(x)) for x in GOPa_ocr['dictionary_id']]
            
#             GOPa_res = GOPa_res.groupby('GOPa').agg({'relation_id': list, 'dictionary_id': list, 'gene_name': list, 'source':max}).reset_index()
            
#             GOPa_res['n'] = [int(len(x)) for x in GOPa_res['gene_name']]
            
#             GOPa_res = pd.merge(GOPa_res, GOPa_ocr[['GOPa', 'ocr_n']], on = 'GOPa', how = 'left')
        
#             GOPa_res['pct'] = GOPa_res['n']/GOPa_res['ocr_n']
            
            
#             GOPa_res['p-val[BIN]'] = float('nan')
#             GOPa_res['p-val[FISH]'] = float('nan')
            
#             GOPa_res = GOPa_res.reset_index(drop = True)
        
            
        
#             for n, p in enumerate(tqdm(GOPa_res['GOPa'])):  
#                 GOPa_res['p-val[BIN]'][n] = stats.binomtest(int(GOPa_res['n'][n]), int(input_genes), float(GOPa_res['ocr_n'][n]/total_genes), alternative='greater').pvalue
#                 observed_genes = int(GOPa_res['n'][n])
#                 not_observed_genes = int(GOPa_res['ocr_n'][n]) - observed_genes
#                 ontingency_table = [[observed_genes, not_observed_genes], [input_genes*(int(GOPa_res['ocr_n'][n])/total_genes), total_genes - (input_genes*(int(GOPa_res['ocr_n'][n])/total_genes))]]
#                 odds_ratio, GOPa_res['p-val[FISH]'][n] = stats.fisher_exact(ontingency_table, alternative='greater')
        
#             GOPa_res['p-adj[BIN-BF]'] = GOPa_res['p-val[BIN]'] * len(GOPa_res['p-val[BIN]'])
#             GOPa_res['p-adj[BIN-BF]'][GOPa_res['p-adj[BIN-BF]'] >= 1] = 1
#             GOPa_res['p-adj[FISH-BF]'] = GOPa_res['p-val[FISH]'] * len(GOPa_res['p-val[FISH]'])
#             GOPa_res['p-adj[FISH-BF]'][GOPa_res['p-adj[FISH-BF]'] >= 1] = 1
            
#             GOPa_res = GOPa_res.sort_values(by='p-val[BIN]',  ascending=True)
        
#             n = len(GOPa_res['p-val[BIN]'])
        
#             GOPa_res['p-adj[BIN-FDR]'] = (GOPa_res['p-val[BIN]'] * n) / np.arange(1, n+1)
            
#             GOPa_res = GOPa_res.sort_values(by='p-val[FISH]',  ascending=True)
        
#             GOPa_res['p-adj[FISH-FDR]'] = (GOPa_res['p-val[FISH]'] * n) / np.arange(1, n+1)
            
#             GOPa_res['p-adj[FISH-FDR]'][GOPa_res['p-adj[FISH-FDR]'] >= 1] = 1
#             GOPa_res['p-adj[BIN-FDR]'][GOPa_res['p-adj[BIN-FDR]'] >= 1] = 1
        
                    
#             GOPa_data['GOPa'] = GOPa_res.to_dict(orient = 'list')
         
#             return GOPa_data
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
      
    
    
    
    
#     def gopa_interaction_analysis(GOPa_data):
        
#         """
#         This function conducts interaction analysis on GOPa_data.
    
#         Args:
#            GOPa_data (dict) - GOPa_data from gopa_analysis function
          
    
#         Returns:
#            dict: A GOPa_data results with data interactions
#        """
       
#         print('\n')
#         print('Analysis of interactions of GOPa data...')
    
       
#         try:
        
#             GOPa_res = pd.DataFrame(GOPa_data['GOPa_interactions'])
#             GOPa_enr = pd.DataFrame(GOPa_data['GOPa'])['relation_id']
            
#             GOPa_list = []
#             for x in GOPa_enr:
#                 GOPa_list = GOPa_list + x
            
#             del GOPa_enr
        
        
#             #grouping variables
            
#             GOPa_is_a_ids = GOPa_res[['GO_id','is_a_ids']].explode('is_a_ids')
            
#             GOPa_is_a_ids = GOPa_is_a_ids[GOPa_is_a_ids['is_a_ids'].isin(GOPa_list)]
            
#             GOPa_is_a_ids.columns = ['A', 'B']
            
#             GOPa_is_a_ids['color'] = 'gray'
            
            
        
#             GOPa_part_of_ids = GOPa_res[['GO_id','part_of_ids']].explode('part_of_ids')
            
#             GOPa_part_of_ids = GOPa_part_of_ids[GOPa_part_of_ids['part_of_ids'].isin(GOPa_list)]
            
#             GOPa_part_of_ids.columns = ['A', 'B']
            
#             GOPa_part_of_ids['color'] = 'gray'
        
        
            
#             GOPa_has_part_ids = GOPa_res[['GO_id','has_part_ids']].explode('has_part_ids')
            
#             GOPa_has_part_ids = GOPa_has_part_ids[GOPa_has_part_ids['has_part_ids'].isin(GOPa_list)]
        
#             GOPa_has_part_ids.columns = ['A', 'B']
            
#             GOPa_has_part_ids['color'] = 'gray'
        
            
          
#             #path and disease
            
#             #
            
#             GOPa_disease_ids = GOPa_res[['GO_id','disease_ids']].explode('disease_ids')
            
#             GOPa_disease_ids = GOPa_disease_ids[GOPa_disease_ids['disease_ids'].isin(GOPa_list)]
            
            
#             GOPa_disease_ids.columns = ['A', 'B']
            
#             GOPa_disease_ids['color'] = 'gray'
        
#             #
            
        
#             GOPa_path_ids = GOPa_res[['GO_id','path_ids']].explode('path_ids')
            
#             GOPa_path_ids = GOPa_path_ids[GOPa_path_ids['path_ids'].isin(GOPa_list)]
            
            
#             GOPa_path_ids.columns = ['A', 'B']
            
#             GOPa_path_ids['color'] = 'gray'
            
            
#             # GOPa_network = pd.concat([GOPa_is_a_ids, GOPa_part_of_ids, GOPa_has_part_ids, GOPa_disease_ids, GOPa_path_ids])
            
            
            
            
#             #color variables 
        
#             GOPa_regulates_ids = GOPa_res[['GO_id','regulates_ids']].explode('regulates_ids')
            
#             GOPa_regulates_ids = GOPa_regulates_ids[GOPa_regulates_ids['regulates_ids'].isin(GOPa_list)]
            
#             GOPa_regulates_ids.columns = ['A', 'B']
            
#             GOPa_regulates_ids['color'] = 'gold'
            
            
#             GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['A'].isin([None])]
#             GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['B'].isin([None])]
            
            
#             GOPa_regulates_ids['regulation'] = GOPa_regulates_ids['A'] + GOPa_regulates_ids['B']
        
        
#             #
            
#             GOPa_negatively_regulates_ids = GOPa_res[['GO_id','negatively_regulates_ids']].explode('negatively_regulates_ids')
        
#             GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[GOPa_negatively_regulates_ids['negatively_regulates_ids'].isin(GOPa_list)]
        
              
#             GOPa_negatively_regulates_ids.columns = ['A', 'B']
            
#             GOPa_negatively_regulates_ids['color'] = 'red'
    
            
#             GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['A'].isin([None])]
#             GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['B'].isin([None])]
            
#             GOPa_negatively_regulates_ids['regulation'] = GOPa_negatively_regulates_ids['A'] + GOPa_negatively_regulates_ids['B']
            
#             #
        
#             GOPa_positively_regulates_ids = GOPa_res[['GO_id','positively_regulates_ids']].explode('positively_regulates_ids')
        
#             GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[GOPa_positively_regulates_ids['positively_regulates_ids'].isin(GOPa_list)]
              
#             GOPa_positively_regulates_ids.columns = ['A', 'B']
            
#             GOPa_positively_regulates_ids['color'] = 'green'
            
            
#             GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['A'].isin([None])]
#             GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['B'].isin([None])]
            
#             GOPa_positively_regulates_ids['regulation'] = GOPa_positively_regulates_ids['A'] + GOPa_positively_regulates_ids['B']
        
        
#             GOPa_network = pd.concat([GOPa_is_a_ids, GOPa_part_of_ids, GOPa_has_part_ids, GOPa_disease_ids, GOPa_path_ids, GOPa_positively_regulates_ids, GOPa_negatively_regulates_ids, GOPa_regulates_ids])
    
#             GOPa_network = GOPa_network[~GOPa_network['A'].isin([None])]
#             GOPa_network = GOPa_network[~GOPa_network['B'].isin([None])]
            
            
#             GOPa_network['regulation'] = GOPa_network['A'] + GOPa_network['B']
    
            
#             GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_regulates_ids['regulation']))] = 'gold'
#             GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_negatively_regulates_ids['regulation']))] = 'red'
#             GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_positively_regulates_ids['regulation']))] = 'green'
        
#             GOPa_network = GOPa_network.drop('regulation', axis = 1)
            
#             GOPa_network = GOPa_network.drop_duplicates()
            
#             GOPa_network = GOPa_network.to_dict(orient = 'list')
        
#             GOPa_data['GOPa_interactions'] = GOPa_network
            
#             return GOPa_data
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def gopa_specificity_analysis(GOPa_data, GOPa_metadata):
        
#         """
#         This function conducts statistical / overrepresentation analysis for 
#         potential blood markers, tissue / cell specificity and cellular localization of genes / protein.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function
         
    
#         Returns:
#            dict: A GOPa_data results with appropriate statistics
#        """
       
#         print('\n')
#         print('Overrepresentation analysis of tissue, cellular and location specificity data...')
        
        
#         try:
            
#             GOPa_metadata2 = copy.deepcopy(GOPa_metadata)
#             GOPa_res = GOPa_data['GOPa_specificity']
#             GOPa_metadata2 = GOPa_metadata2['GOPa_specificity']
        
#             # SEQ
#             SEQ = pd.DataFrame(GOPa_res['SEQ'])
#             SEQ_meta = pd.DataFrame(GOPa_metadata2['SEQ'])
        
#             SEQ = SEQ[SEQ['distribution'] != 'Detected in all']
            
#             SEQ = SEQ.reset_index(drop = True)
            
#             input_genes = int(len(set(SEQ['gene_name'])))
            
#             total_genes = int(len(set(SEQ_meta['gene_name'])))
        
#             SEQ_meta = SEQ_meta.groupby('name').agg({'dictionary_id': list}).reset_index()
#             SEQ_meta['ocr_n'] = [int(len(x)) for x in SEQ_meta['dictionary_id']]
            
#             SEQ = SEQ.groupby('name').agg({'gene_name': list, 'dictionary_id': list, 'distribution':list, 'standarized_specificity_score': list,  'TPM': list, 'log_standarized_specificity_score': list,  'log_TPM': list, 'source':list}).reset_index()
            
#             SEQ['n'] = [int(len(x)) for x in SEQ['dictionary_id']]
            
#             SEQ = pd.merge(SEQ, SEQ_meta[['name', 'ocr_n']], on = 'name', how = 'left')
        
#             SEQ['pct'] = SEQ['n']/SEQ['ocr_n']
            
            
#             SEQ['p-val[BIN]'] = float('nan')
#             SEQ['p-val[FISH]'] = float('nan')
            
#             SEQ = SEQ.reset_index(drop = True)
            
#             for n, p in enumerate(tqdm(SEQ['name'])):  
#                 SEQ['p-val[BIN]'][n] = stats.binomtest(int(SEQ['n'][n]), int(input_genes), float(SEQ['ocr_n'][n]/total_genes), alternative='greater').pvalue
#                 observed_genes = int(SEQ['n'][n])
#                 not_observed_genes = int(SEQ['ocr_n'][n]) - observed_genes
#                 ontingency_table = [[observed_genes, not_observed_genes], [input_genes*(int(SEQ['ocr_n'][n])/total_genes), total_genes - (input_genes*(int(SEQ['ocr_n'][n])/total_genes))]]
#                 odds_ratio, SEQ['p-val[FISH]'][n] = stats.fisher_exact(ontingency_table, alternative='greater')
                
        
#             SEQ['p-adj[BIN-BF]'] = SEQ['p-val[BIN]'] * len(SEQ['p-val[BIN]'])
#             SEQ['p-adj[BIN-BF]'][SEQ['p-adj[BIN-BF]'] >= 1] = 1
#             SEQ['p-adj[FISH-BF]'] = SEQ['p-val[FISH]'] * len(SEQ['p-val[FISH]'])
#             SEQ['p-adj[FISH-BF]'][SEQ['p-adj[FISH-BF]'] >= 1] = 1
            
#             SEQ = SEQ.sort_values(by='p-val[BIN]',  ascending=True)
        
#             n = len(SEQ['p-val[BIN]'])
        
#             SEQ['p-adj[BIN-FDR]'] = (SEQ['p-val[BIN]'] * n) / np.arange(1, n+1)
            
#             SEQ = SEQ.sort_values(by='p-val[FISH]',  ascending=True)
        
#             SEQ['p-adj[FISH-FDR]'] = (SEQ['p-val[FISH]'] * n) / np.arange(1, n+1)
            
#             SEQ['p-adj[FISH-FDR]'][SEQ['p-adj[FISH-FDR]'] >= 1] = 1
#             SEQ['p-adj[BIN-FDR]'][SEQ['p-adj[BIN-FDR]'] >= 1] = 1
                    
#             GOPa_res['SEQ'] = SEQ.to_dict(orient = 'list')
            
#             #location
#             location = pd.DataFrame(GOPa_res['location'])
#             location_meta = pd.DataFrame(GOPa_metadata2['location'])
        
#             input_genes = int(len(set(location['gene_name'])))
            
#             total_genes = int(len(set(location_meta['gene_name'])))
        
#             location_meta = location_meta.groupby('location').agg({'dictionary_id': list}).reset_index()
#             location_meta['ocr_n'] = [int(len(x)) for x in location_meta['dictionary_id']]
            
#             location = location.groupby('location').agg({'gene_name': list, 'dictionary_id': list, 'primary_location':list}).reset_index()
            
#             location['n'] = [int(len(x)) for x in location['dictionary_id']]
            
#             location = pd.merge(location, location_meta[['location', 'ocr_n']], on = 'location', how = 'left')
        
#             location['pct'] = location['n']/location['ocr_n']
            
            
#             location['p-val[BIN]'] = float('nan')
#             location['p-val[FISH]'] = float('nan')
            
#             location = location.reset_index(drop = True)
            
#             for n, p in enumerate(tqdm(location['location'])):  
#                 location['p-val[BIN]'][n] = stats.binomtest(int(location['n'][n]), int(input_genes), float(location['ocr_n'][n]/total_genes), alternative='greater').pvalue
#                 observed_genes = int(location['n'][n])
#                 not_observed_genes = int(location['ocr_n'][n]) - observed_genes
#                 ontingency_table = [[observed_genes, not_observed_genes], [input_genes*(int(location['ocr_n'][n])/total_genes), total_genes - (input_genes*(int(location['ocr_n'][n])/total_genes))]]
#                 odds_ratio, location['p-val[FISH]'][n] = stats.fisher_exact(ontingency_table, alternative='greater')
                
        
#             location['p-adj[BIN-BF]'] = location['p-val[BIN]'] * len(location['p-val[BIN]'])
#             location['p-adj[BIN-BF]'][location['p-adj[BIN-BF]'] >= 1] = 1
#             location['p-adj[FISH-BF]'] = location['p-val[FISH]'] * len(location['p-val[FISH]'])
#             location['p-adj[FISH-BF]'][location['p-adj[FISH-BF]'] >= 1] = 1
            
#             location = location.sort_values(by='p-val[BIN]',  ascending=True)
        
#             n = len(location['p-val[BIN]'])
        
#             location['p-adj[BIN-FDR]'] = (location['p-val[BIN]'] * n) / np.arange(1, n+1)
            
#             location = location.sort_values(by='p-val[FISH]',  ascending=True)
        
#             location['p-adj[FISH-FDR]'] = (location['p-val[FISH]'] * n) / np.arange(1, n+1)
            
#             location['p-adj[FISH-FDR]'][location['p-adj[FISH-FDR]'] >= 1] = 1
#             location['p-adj[BIN-FDR]'][location['p-adj[BIN-FDR]'] >= 1] = 1
                    
#             GOPa_res['location'] = location.to_dict(orient = 'list')
            
            
#             GOPa_data['GOPa_specificity'] = GOPa_res
            
         
#             return GOPa_data
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
      
    
    
    
#     def select_test(test, adj):
#         try:
#             test_string = ''
    
#             if adj != None and adj.upper() in ['BF','FDR']:
#                 test_string = test_string + 'p-adj['
#             else:
#                 test_string = test_string + 'p-val['
            
            
#             if test != None and test.upper() == 'BIN':
#                 test_string = test_string + 'BIN'
#             elif test != None and test.upper() == 'FISH':
#                 test_string = test_string + 'FISH'
#             else:
#                 test_string = test_string + 'BIN'
                
            
#             if adj != None and adj.upper() == 'BF':
#                 test_string = test_string + '-BF]'
#             elif adj != None and adj.upper() == 'FDR':
#                 test_string = test_string + '-FDR]'
#             else:
#                 test_string = test_string + ']'
            
#             return test_string
#         except:
#             print('\n')
#             print('Provided wrong test input!')
            
        
               
    
        
    
#     def GOPa_bar_plot(GOPa_data, GOPa_metadata, p_val = 0.05, test = 'FISH', adj = 'FDR', n = 25, side = 'right', color = 'blue', width = 10, bar_width = 0.5, count_type = 'p_val', details = 0.75, omit = None):
        
#         """
#         This function creates a bar plot of statistical / overrepresentation analysis terms from GO-TERM, PATHWAYS, DISEASES, and VIRAL (diseases related to viral infection)
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function
#            p_val (float) - value of minimal p_val for statistical test
#            test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
#            adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
#            n (int) - maximal number of bars on the graph
#            side (str) - orientation of bars ['left' / 'right']
#            color (str)- color of the bars
#            width (float) - width of the graph
#            bar_width (float) - width of the bars
#            count_type (str) - type of amount of term representation on bars ['perc' - percent representation / 'p_val' - -log(p-val) of selected statistic test / 'num' - number representation]
#            details (float) - degree detail for display GO-TERM and PATHWAYS where 1 is the maximal value and 0.1 is minimal [0.1-1]
#            omit (list) - type of terms to omit in the graph eg. ['GO-TERM', 'KEGG', 'REACTOME', 'DISEASES', 'VIRAL'] or None
         
#         Returns:
#            dict of graphs: Dictionary of bar plots of overrepresentation analysis for GO-TERMS, PATHWAYS, DISEASES, and VIRAL
#         """
        
#         try:
            
#             test_string = select_test(test, adj)
#             figure_dict = {}
            
            
#             sets = {'GO-TERM':['GO-TERM'], 'PATHWAYS':['KEGG', 'REACTOME'], 'DISEASES':['DISEASES'], 'VIRAL-INFECTIONS':['VIRAL']}
            
#             GOPa = pd.DataFrame(GOPa_data['GOPa'])
            
#             if omit != None:
#                 try:
#                     GOPa = GOPa[~GOPa['source'].isin(omit)]
#                 except:
#                     GOPa = GOPa[GOPa['source'] != omit]
            
#             genes_number = int(len(set(GOPa.explode('gene_name')['gene_name'])))
            
#             k_set = []
#             for nk, k in enumerate(sets.keys()):
#                 if len(GOPa[GOPa['source'].isin(sets[k])].copy()) > 0:
#                     k_set.append(k)
                    
#             sets = {key: sets[key] for key in k_set}
                
            
#             GOPa = GOPa[GOPa[str(test_string)] <= p_val]
            
            
#             # Create subplots in a single row
            
#             for nk, k in enumerate(sets.keys()):
                
                    
#                     tmp = GOPa[GOPa['source'].isin(sets[k])].copy()
#                     tmp[test_string] = tmp[test_string] + np.min(tmp[test_string][tmp[test_string] != 0])/2
#                     tmp['-log(p-val)'] = -np.log(tmp[test_string])
#                     tmp['%'] = tmp['n'] / genes_number * 100
                    
                    
#                     detailed_tmp = pd.DataFrame(GOPa_metadata['GOPa'])
#                     detailed_tmp = detailed_tmp[detailed_tmp['source'].isin(sets[k])].copy().drop_duplicates()
#                     detailed_tmp = detailed_tmp.drop_duplicates()
                    
#                     if k in ['PATHWAYS', 'GO-TERM'] and details != None:
#                         detailed_tmp = Counter(list(detailed_tmp['GOPa']))
                        
#                         detailed_tmp = pd.DataFrame(detailed_tmp.items(), columns=['GOPa', 'n'])
                                            
#                         detailed_tmp = detailed_tmp.reset_index(drop = True)
#                         detailed_tmp = list(detailed_tmp['GOPa'][detailed_tmp['n'] < np.quantile(detailed_tmp['n'][detailed_tmp['n'] > 100], details)])
                        
#                         tmp = tmp[tmp['GOPa'].isin(detailed_tmp)]
    
                    
                
#                     # Create a horizontal bar plot
#                     if count_type.upper() == 'perc'.upper():
#                         tmp = tmp.sort_values(by='n', ascending=False)
#                         tmp = tmp.reset_index(drop=True)
#                         tmp = tmp.iloc[0:n,:]
                        
#                         height = float(len(tmp['GOPa'])/2.5)
                        
#                         fig, ax = plt.subplots(figsize=(width, height))
                        
#                         ax.barh(tmp['GOPa'], tmp['%'], color=color, height = bar_width)
#                         ax.set_xlabel('Percentr of genes [%]')
                        
#                     elif count_type.upper() == 'p_val'.upper():
#                         tmp = tmp.sort_values(by='-log(p-val)', ascending=False)
#                         tmp = tmp.reset_index(drop=True)
#                         tmp = tmp.iloc[0:n,:]
                        
#                         height = float(len(tmp['GOPa'])/2.5)
                        
#                         fig, ax = plt.subplots(figsize=(width, height))
                        
#                         ax.barh(tmp['GOPa'], tmp['-log(p-val)'], color=color, height = bar_width)
#                         ax.set_xlabel('-log(p-val)')
                        
#                     else:
#                         tmp = tmp.sort_values(by='n', ascending=False)
#                         tmp = tmp.reset_index(drop=True)
#                         tmp = tmp.iloc[0:n,:]
                        
#                         height = float(len(tmp['GOPa'])/2.5)
                        
#                         fig, ax = plt.subplots(figsize=(width, height))
                        
#                         ax.barh(tmp['GOPa'], tmp['n'], color=color, height = bar_width)
#                         ax.set_xlabel('Number of genes')
                        
                
#                     # Set labels and title
                    
#                     ax.set_ylabel('')
#                     ax.set_title(k)
                
#                     # Invert y-axis to have names on the right
#                     ax.invert_yaxis()
             
                
#                     if side == 'right':
#                         ax.yaxis.tick_right()
#                         ax.set_yticks(range(len(tmp)))
#                     elif side == 'left':
#                         ax.invert_xaxis()
                    
                     
                    
#                     figure_dict[k] = fig
                    
#             return figure_dict
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
        
        
    
    
#     def save_plots(images:dict, path = '', prefix = '', img_format = 'svg'):
#         """
#         This function saves the plots which are included in the dictionary.
        
        
#         Args:
#            images (dict) - dictionary of graphs where dictionary.keys() will part of saved graph name
#            path (str) - path to save the graphs
#            prefix (str) - prefix for the saved graph names
#            img_format (str) - format of saved graphs ['svg' / 'png' / 'jpg']
          
          
#         Returns:
#            file: Saved graphs in the indicated directory
#         """
        
#         try:
            
#             for i in images.keys():
#                 images[i].savefig(path + prefix + '_' + str(i) + '.' + img_format,  format = img_format, bbox_inches = 'tight')
                
                
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def GOPa_network_vis(GOPa_data, GOPa_metadata, p_val = 0.05, test = 'FISH', adj = 'FDR', n_max = 10, list_of_terms = None, omit = None, details = 0.75, path = _path_tmp):
        
#         """
#         This function creates a visualization of the GOPa terms connections in the network format.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function
#            p_val (float) - value of minimal p_val for statistical test
#            test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
#            adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
#            n_max (int) - maximum number of interactions for each term
#            list_of_terms (list) - list of terms to visualization their interactions ['term1', 'term2'] or None for highly occurrence interactions in GOPa_data
#            omit (list) - type of terms to omit in the graph eg. ['KEGG', 'REACTOME', 'DISEASES', 'VIRAL'] or None
#            details (float) - degree detail for display GO-TERM and PATHWAYS where 1 is the maximal value and 0.1 is minimal [0.1-1]. Work only if list_of_terms = None
#            path (str) - path to temporarily save the visualization
           
           
#         Returns:
#            graph: Network graph for GO-TERMS, PATHWAYS, DISEASES, and VIRAL interactions
#         """
        
#         try:
            
#             #screen parameter 
            
#             root = tk.Tk()
#             screen_width = root.winfo_screenwidth()
#             screen_height = root.winfo_screenheight()
#             root.destroy()
            
#             # Calculate desired height and width based on screen size
#             desired_height = int(screen_height * 0.85)  
#             desired_width = int(screen_width * 0.7)  
            
            
#             #data network prepare
        
        
#             test_string = select_test(test, adj)
            
            
#             GOPa_interactions = pd.DataFrame(GOPa_data['GOPa_interactions'])
            
#             GOPa = pd.DataFrame(GOPa_data['GOPa'])
#             GOPa = GOPa[GOPa[str(test_string)] <= p_val]
            
#             if omit != None:
#                 try:
#                     GOPa = GOPa[~GOPa['source'].isin(omit)]
#                 except:
#                     GOPa = GOPa[GOPa['source'] != omit]
            
            
            
#             GOPa = GOPa[['GOPa','relation_id', 'source', 'ocr_n']].explode('relation_id')
#             GOPa = GOPa.drop_duplicates()
            
#             GOPa['color'] = float('nan')
#             GOPa['color'][GOPa['source'] == 'GO-TERM'] = 'lightblue'
#             GOPa['color'][GOPa['source'].isin(['REACTOME', 'KEGG'])] = 'orange'
#             GOPa['color'][GOPa['source'] == 'DISEASES'] = 'purple'
#             GOPa['color'][GOPa['source'] == 'VIRAL'] = 'gray'
            
#             if list_of_terms != None:
#                 GOPa['color'][GOPa['GOPa'].isin(list_of_terms)] = 'aqua'
        
        
            
            
#             GOPa_interactions = GOPa_interactions[GOPa_interactions['A'].isin(list(GOPa['relation_id']))]
#             GOPa_interactions = GOPa_interactions[GOPa_interactions['B'].isin(list(GOPa['relation_id']))]
            
              
#             name_mapping = dict(zip(GOPa['relation_id'], GOPa['source']))
#             GOPa_interactions['source'] = GOPa_interactions['B'].map(name_mapping)
                
#             GOPa_interactions['source'][GOPa_interactions['source'].isin(['KEGG','REACTOME'])] = 'PATHWAYS'
#             GOPa['source'][GOPa['source'].isin(['KEGG','REACTOME'])] = 'PATHWAYS'
    
    
     
#             name_mapping = dict(zip(GOPa['relation_id'], GOPa['GOPa']))
#             GOPa_interactions['A'] = GOPa_interactions['A'].map(name_mapping)
#             GOPa_interactions['B'] = GOPa_interactions['B'].map(name_mapping)
            
            
    
            
#             interactions_df = pd.DataFrame()
#             if list_of_terms != None:
                
#                 for tr in list_of_terms:
#                     if True in list(GOPa_interactions['A'].isin([tr])):
                        
#                         tmp = GOPa_interactions[GOPa_interactions['A'].isin([tr])]
                        
                        
#                         B = Counter(list(tmp['B']))
                        
#                         B = pd.DataFrame(B.items(), columns=['GOPa', 'n'])
                        
#                         B = B.sort_values(by = 'n', ascending = False)
                        
#                         B = B.reset_index(drop = True)
    
#                         if len(B['n']) > n_max:
#                             tmp = tmp[tmp['B'].isin(B['GOPa'][B['n'] >= int(B['n'][math.ceil(n_max/len(set(tmp['source'])))])])]
#                         else:
#                             tmp = tmp
                        
    
#                         tmp2 = GOPa_interactions[GOPa_interactions['B'].isin(list(tmp['B']))]
#                         tmp3 = GOPa_interactions[GOPa_interactions['A'].isin(list(tmp2['A']))]
                        
                        
#                         B = Counter(list(tmp3['B']))
                        
#                         B = pd.DataFrame(B.items(), columns=['GOPa', 'n'])
                        
#                         B = B.sort_values(by = 'n', ascending = False)
                        
#                         B = B.reset_index(drop = True)
                        
                        
                        
#                         if len(B['n']) > n_max:
#                             tmp3 = tmp3[tmp3['B'].isin(list(B['GOPa'][B['n'] >= int(B['n'][math.ceil(n_max/len(set(tmp3['source'])))])]))]
                           
#                         else:
#                             tmp3 = tmp3
                            
                            
#                         tmp = pd.concat([tmp, tmp3])
    
#                         interactions_df = pd.concat([interactions_df, tmp])
                        
    
#                     elif True in list(GOPa_interactions['B'].isin([tr])):
                        
#                         tmp = GOPa_interactions[GOPa_interactions['B'].isin([tr])]
                        
                         
#                         A = Counter(list(tmp['A']))
                        
#                         A = pd.DataFrame(A.items(), columns=['GOPa', 'n'])
                        
#                         A = A.sort_values(by = 'n', ascending = False)
                        
#                         A = A.reset_index(drop = True)
                        
                        
#                         if len(A['n']) > n_max:
#                             tmp = tmp[tmp['A'].isin(A['GOPa'][A['n'] >= int(A['n'][math.ceil(n_max/len(set(tmp['source'])))])])]
#                         else:
#                             tmp = tmp
                        
    
#                         tmp2 = GOPa_interactions[GOPa_interactions['A'].isin(list(tmp['A']))]
                        
#                         srcs = list(set(tmp2['source']))
#                         list_of_B_term = []
                        
#                         for s in srcs:
                            
#                             gopa_list = Counter(list(tmp2['B'][tmp2['source'] == s]))
#                             # Create a DataFrame from the Counter dictionary
#                             gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'n'])
                            
#                             gopa_list = gopa_list.sort_values(by = 'n', ascending = False)
                            
#                             gopa_list = gopa_list.reset_index(drop = True)
                                   
#                             if int(n_max) > len(gopa_list['n']):
#                                 n_t = len(gopa_list['n']) - 1
#                             else:
#                                 n_t = math.ceil(len(set(tmp2['source'])))
                            
                            
#                             list_of_B_term = list_of_B_term + list(tmp2['B'][tmp2['B'].isin(list(gopa_list['GOPa'][gopa_list['n'] >= math.ceil(gopa_list['n'][n_t])]))])
                            
                        
#                         tmp2 = tmp2[(tmp2['A'].isin(list(tmp['A']))) & (tmp2['B'].isin(list(list_of_B_term)))]
    
#                         tmp = pd.concat([tmp,tmp2]) 
                        
                       
#                         interactions_df = pd.concat([interactions_df,tmp])
                
                
                
                
                    
#                 A = Counter(list(interactions_df['A']))
                
#                 A = pd.DataFrame(A.items(), columns=['GOPa', 'n'])
                
#                 A = A.sort_values(by = 'n', ascending = False)
                
#                 A = A.reset_index(drop = True)
                
#                 srcs = list(set(tmp2['source']))
#                 list_of_B_term = []
                
#                 for s in srcs:
                    
#                     B = Counter(list(interactions_df['B'][interactions_df['source'] == s]))                
#                     B = pd.DataFrame(B.items(), columns=['GOPa', 'n'])
                    
#                     B = B.sort_values(by = 'n', ascending = False)
                    
#                     B = B.reset_index(drop = True)
                    
#                     if len(B['n']) > n_max:
#                         list_of_B_term = list_of_B_term + list(B['GOPa'][B['n'] >= int(B['n'][math.ceil(n_max/len(set(interactions_df['source'])))])]) 
    
#                     else:
#                         list_of_B_term = list_of_B_term + list(B['GOPa']) 
    
                    
                
                
#                 interactions_df = interactions_df[(interactions_df['A'].isin(list(A['GOPa'][A['n'] >= int(A['n'][math.ceil(n_max/len(set(interactions_df['source'])))])]) + list_of_terms)) & (interactions_df['B'].isin(list_of_B_term+ list_of_terms))]
                  
                    
    
                
#             else:
                
#                 GOPa = GOPa[GOPa['GOPa'].isin(list(GOPa_interactions['A']) + list(GOPa_interactions['B']))]
#                 for nk, k in enumerate(set(GOPa['source'])):
#                     tmp = GOPa[GOPa['source'].isin([k])]
                          
#                     detailed_tmp = pd.DataFrame(GOPa_metadata['GOPa'])
#                     detailed_tmp['source'][detailed_tmp['source'].isin(['KEGG','REACTOME'])] = 'PATHWAYS'
#                     detailed_tmp = detailed_tmp[detailed_tmp['source'].isin([k])].copy().drop_duplicates()
#                     detailed_tmp = detailed_tmp.drop_duplicates()
                    
#                     if k in ['PATHWAYS', 'GO-TERM'] and details != None:
#                         detailed_tmp = Counter(list(detailed_tmp['GOPa']))
                        
#                         detailed_tmp = pd.DataFrame(detailed_tmp.items(), columns=['GOPa', 'n'])
                                            
#                         detailed_tmp = detailed_tmp.reset_index(drop = True)
#                         detailed_tmp = list(detailed_tmp['GOPa'][detailed_tmp['n'] < np.quantile(detailed_tmp['n'][detailed_tmp['n'] > 100], details)])
                        
#                         tmp = tmp[tmp['GOPa'].isin(detailed_tmp)]
    
    
#                 list_a = []
#                 list_b = []
#                 for src in list(set(GOPa_interactions['source'])):
                    
#                     A = Counter(list(GOPa_interactions['A'][GOPa_interactions['source'] == src]))
                    
#                     # Create a DataFrame from the Counter dictionary
#                     A = pd.DataFrame(A.items(), columns=['GOPa', 'n'])
                    
#                     A = A.sort_values(by = 'n', ascending = False)
                    
#                     A = A.reset_index(drop = True)
                    
#                     if len(A['n']) > n_max:
#                        list_a = list_a + list(GOPa_interactions['A'][GOPa_interactions['A'].isin(A['GOPa'][A['n'] >= int(A['n'][math.ceil(n_max/len(set(GOPa_interactions['source'])))])])])
#                     else:
#                         list_a = list_a + list(GOPa_interactions['A'][GOPa_interactions['source'] == src])
                        
                        
    
#                     B = Counter(list(GOPa_interactions['B'][GOPa_interactions['source'] == src]))
                    
#                     # Create a DataFrame from the Counter dictionary
#                     B = pd.DataFrame(B.items(), columns=['GOPa', 'n'])
                    
#                     B = B.sort_values(by = 'n', ascending = False)
                    
#                     B = B.reset_index(drop = True)
                    
#                     if len(B['n']) > n_max:
#                        list_b = list_b + list(GOPa_interactions['B'][GOPa_interactions['B'].isin(B['GOPa'][B['n'] >= int(B['n'][math.ceil(n_max/len(set(GOPa_interactions['source'])))])])])
#                     else:
#                         list_b =  list_b + list(GOPa_interactions['B'][GOPa_interactions['source'] == src])
    
                
#                 interactions_df = GOPa_interactions[(GOPa_interactions['A'].isin(list_a)) & (GOPa_interactions['B'].isin(list_b))]
                
                
#                 A = Counter(list(interactions_df['A']))
                
#                 # Create a DataFrame from the Counter dictionary
#                 A = pd.DataFrame(A.items(), columns=['GOPa', 'n'])
                
#                 A = A.sort_values(by = 'n', ascending = False)
                
#                 A = A.reset_index(drop = True)
                
#                 interactions_df = interactions_df[interactions_df['A'].isin(A['GOPa'][A['n'] >= int(A['n'][math.ceil(n_max/len(set(interactions_df['source'])))])])]
                
                
#             if len(interactions_df) > 0:
                
                
    
#                 GOPa_interactions = interactions_df
#                 GOPa_interactions = GOPa_interactions.reset_index(drop = True)
                
#                 gopa_list = list(GOPa_interactions['B']) + list(GOPa_interactions['A'])
                
#                 GOPa = GOPa[['GOPa', 'color']].drop_duplicates()
                
#                 GOPa = GOPa[GOPa['GOPa'].isin(gopa_list)]
                
#                 # Count the occurrences of each element in the list
#                 gopa_list = Counter(gopa_list)
                
#                 # Create a DataFrame from the Counter dictionary
#                 gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'weight'])
            
               
#                 GOPa = pd.merge(GOPa, gopa_list, on = 'GOPa', how = 'left')
            
            
                    
#                 G = nx.Graph() 
            
             
#                 for _, row in GOPa.iterrows():
#                     node = row['GOPa']
#                     color = row['color']
#                     weight = np.log2(row['weight']*500)
#                     G.add_node(node, size = weight, color = color)
                    
#                 for index, row in GOPa_interactions.iterrows():
#                     source = row['A']
#                     target = row['B']
#                     color = row['color']
#                     G.add_edge(source, target, color = color)
            
                
                
#                 # Create a pyvis Network instance
#                 net = Network(notebook=True, height=f"{desired_height}px", width=f"{desired_width}px")
                
#                 net.from_nx(G)
                
            
#                 net.repulsion(node_distance=150, spring_length=200)
#                 net.show_buttons(filter_=['nodes', 'physics'])
            
                
#                 net.show(os.path.join(path, 'tmp.html'))
#                 webbrowser.open(os.path.abspath(os.path.join(path, 'tmp.html')))
            
                
#                 return G
            
#             else:
#                 print('\n \n')
        
#                 print('Lack of GO-TERM connections to all provided terms')
                
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
    
#     def create_netowrk_GOPa_static(network, width = 16, height = 10, font_size = 10):
        
#         """
#         This function creates a network graph in static format.
    
#         Args:
#            network (network) - network from GOPa_network_vis or show_go_term functions
#            width (float) - width of the graph
#            height (float) - height of the graph
#            font_size (float) - size of the fonts
          
         
#         Returns:
#            graph: Network in static format with the legend
#         """
        
        
#         try:
#             # Layout
#             pos = nx.spring_layout(network, seed = 123)  # You can choose a different layout algorithm
        
#             # Drawing nodes and edges with attributes
#             node_sizes = [data['size']*10 for node, data in network.nodes(data=True)]
#             node_colors = [data['color'] for node, data in network.nodes(data=True)]
#             edge_colors = [data['color'] for _, _, data in network.edges(data=True)]
        
#             # Create the plot
#             fig, ax = plt.subplots(figsize=(width, height))
#             nx.draw_networkx_nodes(network, pos, node_size=node_sizes, node_color=node_colors)
#             nx.draw_networkx_edges(network, pos, edge_color=edge_colors)
        
#             node_labels = {node: node for node in network.nodes()}  
#             node_labels = {node: node for node in network.nodes()}  
#             labels = nx.draw_networkx_labels(network, pos, labels=node_labels, font_size=font_size, font_color="black", verticalalignment='center')
           
         
#             texts = [labels[node] for node in labels]
#             adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=.5))
            
#             legend_elements = [
#                 Line2D([0], [0], marker='o', color='w', label='TARGET-TERMS',
#                        markerfacecolor='aqua', markersize=10),
#                 Line2D([0], [0], marker='o', color='w', label='GO-TERM',
#                        markerfacecolor='lightblue', markersize=10),
#                 Line2D([0], [0], marker='o', color='w', label='PATHWAYS',
#                        markerfacecolor='orange', markersize=10),
#                 Line2D([0], [0], marker='o', color='w', label='DISEASES',
#                        markerfacecolor='plum', markersize=10),
#                 Line2D([0], [0], marker='o', color='w', label='VIRAL-DISEASES',
#                        markerfacecolor='gray', markersize=10),
#                 Line2D([], [], linestyle='-', color='gold', label='regulate', linewidth=1),
#                 Line2D([], [], linestyle='-', color='red', label='negatively_regulates', linewidth=1),
#                 Line2D([], [], linestyle='-', color='green', label='positively_regulate', linewidth=1)
#                 ]
          
                      
                       
        
#             ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1))
        
        
               
        
#             plt.axis('off')  # Turn off axis
#             plt.show()
            
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
    
#     def create_netowrk_gene_static(network, width = 12, height = 10, font_size = 10):
        
#         """
#         This function creates a network graph in static format.
    
#         Args:
#            network (network) - network from GOPa_network_vis or show_go_term functions
#            width (float) - width of the graph
#            height (float) - height of the graph
#            font_size (float) - size of the fonts
          
         
#         Returns:
#            graph: Network in static format with the legend
#         """
        
        
#         try:
#             # Layout
#             pos = nx.spring_layout(network, seed = 123)  # You can choose a different layout algorithm
        
#             # Drawing nodes and edges with attributes
#             node_sizes = [data['size']*10 for node, data in network.nodes(data=True)]
#             node_colors = [data['color'] for node, data in network.nodes(data=True)]
#             edge_colors = [data['color'] for _, _, data in network.edges(data=True)]
        
#             # Create the plot
#             fig, ax = plt.subplots(figsize=(width, height))
#             nx.draw_networkx_nodes(network, pos, node_size=node_sizes, node_color=node_colors)
#             nx.draw_networkx_edges(network, pos, edge_color=edge_colors)
        
#             node_labels = {node: node for node in network.nodes()}  
#             node_labels = {node: node for node in network.nodes()}  
#             labels = nx.draw_networkx_labels(network, pos, labels=node_labels, font_size=font_size, font_color="black", verticalalignment='center')
           
         
#             texts = [labels[node] for node in labels]
#             adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=.5))
        
               
        
#             plt.axis('off')  # Turn off axis
#             plt.show()
            
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
    
#     def create_netowrk_html(network, options = None):
        
#         """
#         This function creates a network graph in interactive HTML format with adjustment.
    
#         Args:
#            network (network) - network from GOPa_network_vis or show_go_term functions
#            options (str) - string of options for interactive graph adjustment generated in GOPa_network_vis or show_go_term functions or None
          
         
#         Returns:
#            graph: Network in interactive HTML format
#         """
        
#         try:
#             #screen parameter 
            
#             root = tk.Tk()
#             screen_width = root.winfo_screenwidth()
#             screen_height = root.winfo_screenheight()
#             root.destroy()
            
#             # Calculate desired height and width based on screen size
#             desired_height = int(screen_height * 0.8)  
#             desired_width = int(screen_width * 0.95)  
            
#             #
            
#             net = Network(notebook=True, height=f"{desired_height}px", width=f"{desired_width}px")
        
#             net.from_nx(network)
            
#             if options != None:
#                 net.set_options(options)
            
#             return net
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
            
    
    
    
    
#     def show_go_term(GOPa_metadata, GO_id = None, GO_name = None, n_max = 10, omit = None, path = _path_tmp):
        
#         """
#         This function creates a visualization of the GOPa terms connections in the network format.
    
#         Args:
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function
#            GO_id (str) - id of GO-TERM or None
#            GO_name (str) - name of GO-TERM or None
#            n_max (int) - maximum number of interactions for the term
#            path (str) - path to temporarily save the visualization
#            omit (list) - type of terms to omit in the graph eg. ['KEGG', 'REACTOME', 'DISEASES', 'VIRAL'] or None
    
          
         
#         Returns:
#            graph: Network for GO-TERM with its interactions
#         """
        
#         try:
        
#             #screen parameter 
            
#             root = tk.Tk()
#             screen_width = root.winfo_screenwidth()
#             screen_height = root.winfo_screenheight()
#             root.destroy()
            
#             # Calculate desired height and width based on screen size
#             desired_height = int(screen_height * 0.85)  
#             desired_width = int(screen_width * 0.7)  
            
            
#             #go-term network prepare
        
        
            
#             GOPa_metadata2 = pd.DataFrame(copy.deepcopy(GOPa_metadata['GOPa_interactions']))
            
#             GOPa_metadata2['name_space'] = [x.upper() for x in GOPa_metadata2['name_space']]
            
#             name_mapping = dict(zip(['BIOLOGICAL_PROCESS', 'CELLULAR_COMPONENT', 'MOLECULAR_FUNCTION'], ['BP : ','CC : ','MF : ']))
            
            
#             GOPa_metadata2['role_sc'] = GOPa_metadata2['name_space'].map(name_mapping)
    
            
#             prefixes = ['rRNA','mRNA','snoRNA','lncRNA','piRNA','tRNA','miRNA','snRNA','siRNA', 'cGMP', 'mTOR', 'cAMP' , 'vRNA', 'snRNP']
#             GOPa_metadata2['name'] = [x[0].upper() + x[1:] if not any(x[:4] in prefix for prefix in prefixes) else x for x in GOPa_metadata2['name']]
#             GOPa_metadata2['name'] = GOPa_metadata2['role_sc'] + GOPa_metadata2['name']
            
#             GOPa = pd.DataFrame(GOPa_metadata['GOPa'])
        
            
#             if GO_id != None:
            
#                 GOPa_metadata2 = GOPa_metadata2[GOPa_metadata2['GO_id'] == GO_id]
            
#             elif GO_name != None:
            
#                 GOPa_metadata2 = GOPa_metadata2[GOPa_metadata2['name'] == GO_name]
#             else:
#                 GOPa_metadata2 = None
        
#             try:
#                 #grouping variables
                
#                 GOPa_is_a_ids = GOPa_metadata2[['GO_id','is_a_ids']].explode('is_a_ids')
                    
#                 GOPa_is_a_ids.columns = ['A', 'B']
                    
                
            
#                 GOPa_part_of_ids = GOPa_metadata2[['GO_id','part_of_ids']].explode('part_of_ids')
                
#                 GOPa_part_of_ids.columns = ['A', 'B']
            
            
            
                
#                 GOPa_has_part_ids = GOPa_metadata2[['GO_id','has_part_ids']].explode('has_part_ids')
                
#                 GOPa_has_part_ids.columns = ['A', 'B']
            
                
            
                
#                 GOPa_regulates_ids = GOPa_metadata2[['GO_id','regulates_ids']].explode('regulates_ids')
                    
#                 GOPa_regulates_ids.columns = ['A', 'B']
                
                
                
#                 GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['A'].isin([None])]
#                 GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_regulates_ids['regulation'] = GOPa_regulates_ids['A'] + GOPa_regulates_ids['B']
            
            
#                 #
                
#                 GOPa_negatively_regulates_ids = GOPa_metadata2[['GO_id','negatively_regulates_ids']].explode('negatively_regulates_ids')
            
#                 GOPa_negatively_regulates_ids.columns = ['A', 'B']
                
                
#                 GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['A'].isin([None])]
#                 GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_negatively_regulates_ids['regulation'] = GOPa_negatively_regulates_ids['A'] + GOPa_negatively_regulates_ids['B']
                
#                 #
            
#                 GOPa_positively_regulates_ids = GOPa_metadata2[['GO_id','positively_regulates_ids']].explode('positively_regulates_ids')
                  
#                 GOPa_positively_regulates_ids.columns = ['A', 'B']
                
                            
#                 GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['A'].isin([None])]
#                 GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_positively_regulates_ids['regulation'] = GOPa_positively_regulates_ids['A'] + GOPa_positively_regulates_ids['B']
                
                
#                 GOPa_network = pd.concat([GOPa_is_a_ids, GOPa_part_of_ids, GOPa_has_part_ids, GOPa_positively_regulates_ids, GOPa_negatively_regulates_ids, GOPa_regulates_ids])
                
#                 GOPa_network = GOPa_network[~GOPa_network['A'].isin([None])]
#                 GOPa_network = GOPa_network[~GOPa_network['B'].isin([None])]
            
            
                
            
#                 det_list = set(list(list(GOPa_network['A']) + list(GOPa_network['B'])))
                
#                 GOPa_metadata2 = pd.DataFrame(copy.deepcopy(GOPa_metadata['GOPa_interactions']))
                
#                 GOPa = pd.DataFrame(GOPa_metadata['GOPa'])
            
            
#                 GOPa_metadata2 = GOPa_metadata2[GOPa_metadata2['GO_id'].isin(det_list)]
                    
            
#                 #grouping variables
                
#                 GOPa_is_a_ids = GOPa_metadata2[['GO_id','is_a_ids']].explode('is_a_ids')
                    
#                 GOPa_is_a_ids.columns = ['A', 'B']
                
                
#                 GOPa_is_a_ids['color'] = 'gray'
                
                
            
#                 GOPa_part_of_ids = GOPa_metadata2[['GO_id','part_of_ids']].explode('part_of_ids')
                
#                 GOPa_part_of_ids.columns = ['A', 'B']
            
                
#                 GOPa_part_of_ids['color'] = 'gray'
            
            
                
#                 GOPa_has_part_ids = GOPa_metadata2[['GO_id','has_part_ids']].explode('has_part_ids')
                
#                 GOPa_has_part_ids.columns = ['A', 'B']
            
                
#                 GOPa_has_part_ids['color'] = 'gray'
            
                
              
#                 #path and disease
                
#                 #
                
#                 GOPa_disease_ids = GOPa_metadata2[['GO_id','disease_ids']].explode('disease_ids')    
                
#                 GOPa_disease_ids.columns = ['A', 'B']
            
                
#                 GOPa_disease_ids['color'] = 'gray'
            
#                 #
                
            
#                 GOPa_path_ids = GOPa_metadata2[['GO_id','path_ids']].explode('path_ids')    
                
#                 GOPa_path_ids.columns = ['A', 'B']
            
                
#                 GOPa_path_ids['color'] = 'gray'
                
                
                    
#                 #color variables 
            
#                 GOPa_regulates_ids = GOPa_metadata2[['GO_id','regulates_ids']].explode('regulates_ids')
                    
#                 GOPa_regulates_ids.columns = ['A', 'B']
                
#                 GOPa_regulates_ids['color'] = 'gold'
    
                
                
#                 GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['A'].isin([None])]
#                 GOPa_regulates_ids = GOPa_regulates_ids[~GOPa_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_regulates_ids['regulation'] = GOPa_regulates_ids['A'] + GOPa_regulates_ids['B']
            
            
#                 #
                
#                 GOPa_negatively_regulates_ids = GOPa_metadata2[['GO_id','negatively_regulates_ids']].explode('negatively_regulates_ids')
            
#                 GOPa_negatively_regulates_ids.columns = ['A', 'B']
                
#                 GOPa_negatively_regulates_ids['color'] = 'red'
    
                
#                 GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['A'].isin([None])]
#                 GOPa_negatively_regulates_ids = GOPa_negatively_regulates_ids[~GOPa_negatively_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_negatively_regulates_ids['regulation'] = GOPa_negatively_regulates_ids['A'] + GOPa_negatively_regulates_ids['B']
                
#                 #
            
#                 GOPa_positively_regulates_ids = GOPa_metadata2[['GO_id','positively_regulates_ids']].explode('positively_regulates_ids')
                  
#                 GOPa_positively_regulates_ids.columns = ['A', 'B']
                
                
#                 GOPa_positively_regulates_ids['color'] = 'green'
                
#                 GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['A'].isin([None])]
#                 GOPa_positively_regulates_ids = GOPa_positively_regulates_ids[~GOPa_positively_regulates_ids['B'].isin([None])]
                
                
#                 GOPa_positively_regulates_ids['regulation'] = GOPa_positively_regulates_ids['A'] + GOPa_positively_regulates_ids['B']
            
            
#                 GOPa_network = pd.concat([GOPa_is_a_ids, GOPa_part_of_ids, GOPa_has_part_ids, GOPa_disease_ids, GOPa_path_ids, GOPa_positively_regulates_ids, GOPa_negatively_regulates_ids, GOPa_regulates_ids])
    
#                 GOPa_network = GOPa_network[~GOPa_network['A'].isin([None])]
#                 GOPa_network = GOPa_network[~GOPa_network['B'].isin([None])]
                
            
#                 GOPa_network['regulation'] = GOPa_network['A'] + GOPa_network['B']
                
                
#                 GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_regulates_ids['regulation']))] = 'gold'
#                 GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_negatively_regulates_ids['regulation']))] = 'red'
#                 GOPa_network['color'][GOPa_network['regulation'].isin(list(GOPa_positively_regulates_ids['regulation']))] = 'green'
            
#                 GOPa_network = GOPa_network.drop('regulation', axis = 1)
                
#                 GOPa_network = GOPa_network.drop_duplicates()
                
    
#                 gopa_list = list(set(list(GOPa_network['A']) + list(GOPa_network['B'])))
                
#                 GOPa = GOPa[GOPa['relation_id'].isin(gopa_list)]    
#                 GOPa = GOPa[['GOPa','relation_id', 'source']].explode('relation_id')
#                 GOPa = GOPa.drop_duplicates()
                
#                 GOPa['color'] = float('nan')
#                 GOPa['color'][GOPa['source'] == 'GO-TERM'] = 'lightblue'
#                 GOPa['color'][GOPa['source'].isin(['REACTOME', 'KEGG'])] = 'orange'
#                 GOPa['color'][GOPa['source'] == 'DISEASES'] = 'purple'
#                 GOPa['color'][GOPa['source'] == 'VIRAL'] = 'gray'
                
#                 if omit != None:
#                     try:
#                         GOPa = GOPa[~GOPa['source'].isin(omit)]
#                     except:
#                         GOPa = GOPa[GOPa['source'] != omit]
                
#                 go_without_gene = pd.DataFrame(get_GO()['connections'])
#                 go_without_gene = go_without_gene[go_without_gene['obsolete'] == False]
                
#                 name_mapping = dict(zip(list(GOPa['relation_id']) + list(go_without_gene['GO_id']), list(GOPa['source']) + list(['GO-TERM']*len(list(go_without_gene['GO_id'])))))
#                 GOPa_network['source'] = GOPa_network['B'].map(name_mapping)
                
                
#                 go_without_gene['name_space'] = [x.upper() for x in go_without_gene['name_space']]
                
#                 name_mapping = dict(zip(['BIOLOGICAL_PROCESS', 'CELLULAR_COMPONENT', 'MOLECULAR_FUNCTION'], ['BP : ','CC : ','MF : ']))
                
                
#                 go_without_gene['role_sc'] = go_without_gene['name_space'].map(name_mapping)
    
                
#                 prefixes = ['rRNA','mRNA','snoRNA','lncRNA','piRNA','tRNA','miRNA','snRNA','siRNA', 'cGMP', 'mTOR', 'cAMP' , 'vRNA', 'snRNP']
#                 go_without_gene['name'] = [x[0].upper() + x[1:] if not any(x[:4] in prefix for prefix in prefixes) else x for x in go_without_gene['name']]
                
                
#                 name_mapping = dict(zip(list(go_without_gene['GO_id']) + list(GOPa['relation_id']),   list(go_without_gene['role_sc'] + go_without_gene['name']) + list(GOPa['GOPa'])))
                
#                 GOPa_network['A'] = GOPa_network['A'].map(name_mapping)
#                 GOPa_network['B'] = GOPa_network['B'].map(name_mapping)
                
#                 del go_without_gene, GOPa_metadata2
                
                
                
#                 GOPa_network = GOPa_network.dropna()
                
                
#                 if omit != None:
#                     try:
#                         GOPa_network = GOPa_network[~GOPa_network['source'].isin(omit)]
#                     except:
#                         GOPa_network = GOPa_network[GOPa_network['source'] != omit]
                        
                        
#                 GOPa_network['source'][GOPa_network['source'].isin(['KEGG','REACTOME'])] = 'PATHWAYS'
                
#                 if len(GOPa_network) > 0:  
                    
                    
#                     srcs = list(set(GOPa_network['source']))
#                     list_of_B_term = []
                    
#                     for s in srcs:
#                         if s != 'GO-TERM':
#                             for a_term in set(GOPa_network['A'][GOPa_network['source'].isin([s])]):
                               
#                                 gopa_list = Counter(list(GOPa_network['B'][GOPa_network['A'] == a_term]))
                                
#                                 # Create a DataFrame from the Counter dictionary
#                                 gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'n'])
                                
#                                 gopa_list = gopa_list.sort_values(by = 'n', ascending = False)
                                
#                                 gopa_list = gopa_list.reset_index(drop = True)
                                       
#                                 if int(n_max/len(set(GOPa_network['source']))) > len(gopa_list['n']):
#                                     n_t = len(gopa_list['n']) - 1
#                                 else:
#                                     n_t = math.ceil(n_max/len(set(GOPa_network['source'])))
                                    
                                    
#                                 list_of_B_term = list_of_B_term + list(gopa_list['GOPa'][gopa_list['n'] >= math.ceil(gopa_list['n'][n_t])])
                            
                    
            
#                     GOPa_network = GOPa_network[GOPa_network['B'].isin(list_of_B_term)]
                    
#                     GOPa_network1 = GOPa_network[GOPa_network['source'] == 'GO-TERM']
#                     GOPa_network2 = GOPa_network[GOPa_network['source'] != 'GO-TERM']
                    
                    
#                     srcs = list(set(GOPa_network2['source']))
#                     list_of_B_term = []
                    
#                     for s in srcs:
                        
#                         gopa_list = Counter(list(GOPa_network2['B'][GOPa_network2['source'] == s]))
#                         # Create a DataFrame from the Counter dictionary
#                         gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'n'])
                        
#                         gopa_list = gopa_list.sort_values(by = 'n', ascending = False)
                        
#                         gopa_list = gopa_list.reset_index(drop = True)
                               
#                         if int(n_max) > len(gopa_list['n']):
#                             n_t = len(gopa_list['n']) - 1
#                         else:
#                             n_t = math.ceil(len(set(GOPa_network2['source'])))
                        
                        
#                         list_of_B_term = list_of_B_term + list(gopa_list['GOPa'][gopa_list['n'] >= math.ceil(gopa_list['n'][n_t])])
                
#                     GOPa_network2 = GOPa_network2[GOPa_network2['B'].isin(list_of_B_term)]
                   
#                     GOPa_network = pd.concat([GOPa_network1, GOPa_network2])
#                     GOPa_network = GOPa_network.reset_index(drop = True)
                    
#                     gopa_list = list(GOPa_network['B']) + list(GOPa_network['A'])
                    
#                     GOPa = GOPa[['GOPa', 'color']].drop_duplicates()
                    
#                     GOPa = GOPa[GOPa['GOPa'].isin(gopa_list)]
                    
#                     # Count the occurrences of each element in the list
#                     gopa_list = Counter(gopa_list)
                    
#                     # Create a DataFrame from the Counter dictionary
#                     gopa_list = pd.DataFrame(gopa_list.items(), columns=['GOPa', 'weight'])
                
                    
#                     GOPa_network = GOPa_network.drop_duplicates()
            
#                     GOPa = pd.merge(GOPa, gopa_list, on = 'GOPa', how = 'right')
#                     GOPa['color'][GOPa['color'] != GOPa['color']] = 'lightblue'
                
                        
#                     G = nx.Graph() 
                
                 
#                     for _, row in tqdm(GOPa.iterrows()):
#                         node = row['GOPa']
#                         color = row['color']
#                         weight = np.log2(row['weight']*500)
#                         G.add_node(node, size = weight, color = color)
                        
#                     for index, row in tqdm(GOPa_network.iterrows()):
#                         source = row['A']
#                         target = row['B']
#                         color = row['color']
#                         G.add_edge(source, target, color = color)
                
                    
                    
#                     # Create a pyvis Network instance
#                     net = Network(notebook=True, height=f"{desired_height}px", width=f"{desired_width}px")
                    
#                     net.from_nx(G)
                    
                
#                     net.repulsion(node_distance=150, spring_length=200)
#                     net.show_buttons(filter_=['nodes', 'physics'])
                
                    
#                     net.show(os.path.join(path, 'tmp.html'))
#                     webbrowser.open(os.path.abspath(os.path.join(path, 'tmp.html')))
                
                    
#                     return G
                
#                 else:
#                     print('\n \n')
#                     print('Lack of GO-TERM or GO-TERM is obsoleted')
#             except:
#                 print('\n')
#                 print('GO-TERM not provided')
                
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
#     def create_netowrk_GO_static(network, width = 12, height = 10, font_size = 10):
        
#         """
#         This function creates a network graph in static format.
    
#         Args:
#            network (network) - network from GOPa_network_vis or show_go_term functions
#            width (float) - width of the graph
#            height (float) - height of the graph
#            font_size (float) - size of the fonts
          
         
#         Returns:
#            graph: Network in static format with the legend
#         """
        
        
#         try:
#             # Layout
#             pos = nx.spring_layout(network, seed = 123)  # You can choose a different layout algorithm
        
#             # Drawing nodes and edges with attributes
#             node_sizes = [data['size']*10 for node, data in network.nodes(data=True)]
#             node_colors = [data['color'] for node, data in network.nodes(data=True)]
#             edge_colors = [data['color'] for _, _, data in network.edges(data=True)]
        
#             # Create the plot
#             fig, ax = plt.subplots(figsize=(width, height))
#             nx.draw_networkx_nodes(network, pos, node_size=node_sizes, node_color=node_colors)
#             nx.draw_networkx_edges(network, pos, edge_color=edge_colors)
        
#             node_labels = {node: node for node in network.nodes()}  
#             node_labels = {node: node for node in network.nodes()}  
#             labels = nx.draw_networkx_labels(network, pos, labels=node_labels, font_size=font_size, font_color="black", verticalalignment='center')
           
         
#             texts = [labels[node] for node in labels]
#             adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=.5))
        
               
        
#             plt.axis('off')  # Turn off axis
#             plt.show()
            
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def gene_interactions_network_vis(GOPa_data, GOPa_metadata, target_interaction = None, n_min = 2, top_n = 25, species = None, display_name = 'human', color = 'darkviolet', path = _path_tmp):
        
#         """
#         This function creates a visualization of the genes/proteins connections in the network format.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function 
#            target_interaction (list) - list of genes / proteins as target for interactions network display
#            n_min (int) - minimal number of interactions for gene / protein
#            top_n (int) - maximal number of top abundant interactions for gene / protein
#            species (str) - species for gene - gene / protein - protein interaction study ['human' / 'mouse']. If None all interactions for both species.
#            display_name (str) - format for gene / protein name display ['human' , 'mouse']
#            color (str) - color of the nodes
#            path (str) - path to temporarily save the visualization
            
           
#         Returns:
#            graph: Network graph for genes/proteins interactions
#         """
        
#         try:
            
#             #screen parameter 
            
#             root = tk.Tk()
#             screen_width = root.winfo_screenwidth()
#             screen_height = root.winfo_screenheight()
#             root.destroy()
            
#             # Calculate desired height and width based on screen size
#             desired_height = int(screen_height * 0.85)  
#             desired_width = int(screen_width * 0.7)  
            
            
#             #data network prepare
            
            
#             GOPa_gene_interactions = pd.DataFrame(GOPa_data['GOPa_gene_interaction'])
            
#             GOPa_names = pd.DataFrame(GOPa_data['gene_dictionary'])
            
            
            
#             name_mapping = dict(zip(GOPa_names['dictionary_id'], GOPa_names['gene_name']))
            
#             del GOPa_names
            
            
#             if target_interaction != None:
#                 target_interaction, not_in = search_genes(target_interaction, GOPa_metadata, species=None)
#                 target_interaction = pd.DataFrame(target_interaction['gene_dictionary'])
#                 try:
#                     GOPa_gene_interactions = GOPa_gene_interactions[(GOPa_gene_interactions['id_1'].isin(list(target_interaction['dictionary_id']))) | (GOPa_gene_interactions['id_2'].isin(list(target_interaction['dictionary_id'])))]
#                     del target_interaction
#                 except:
#                     print('Genes / proteins not found in data!')
            
            
#             GOPa_gene_interactions['id_1'] = GOPa_gene_interactions['id_1'].map(name_mapping)
#             GOPa_gene_interactions['id_2'] = GOPa_gene_interactions['id_2'].map(name_mapping)
            
    
            
#             if species != None and species.upper() == 'human'.upper():
#                GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['species'] == 'Homo sapiens']
#             elif species != None and species.upper() == 'mouse'.upper():
#                GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['species'] == 'Mus musculus']
    
                
            
                    
                    
#             if display_name != None and display_name.upper() == 'human'.upper():
#                 GOPa_gene_interactions['id_1'] = [x.upper() for x in GOPa_gene_interactions['id_1']]
#                 GOPa_gene_interactions['id_2'] = [x.upper() for x in GOPa_gene_interactions['id_2']]
    
#             elif display_name != None and display_name.upper() == 'mouse'.upper():
#                 GOPa_gene_interactions['id_1'] = [x[0].upper() + x[1:].lower() for x in GOPa_gene_interactions['id_1']]
#                 GOPa_gene_interactions['id_2'] = [x[0].upper() + x[1:].lower() for x in GOPa_gene_interactions['id_2']]
    
                    
               
               
#             GOPa_gene_interactions = GOPa_gene_interactions[['id_1', 'id_2']].drop_duplicates()
    
    
#             mx_calc = Counter(list(GOPa_gene_interactions['id_1']) + list(GOPa_gene_interactions['id_2']))
            
#             mx_calc = pd.DataFrame(mx_calc.items(), columns=['gene', 'n'])
            
#             mx_calc = mx_calc.sort_values(by = 'n', ascending = False)
            
#             mx_calc = mx_calc.reset_index(drop = True)
            
            
          
#             if len(mx_calc['n']) > top_n:
#                 GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['id_1'].isin(list(mx_calc['gene'][mx_calc['n'] >= int(mx_calc['n'][top_n])]))]
#                 GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['id_2'].isin(list(mx_calc['gene'][mx_calc['n'] >= int(mx_calc['n'][top_n])]))]
                
#             else:
#                 GOPa_gene_interactions = GOPa_gene_interactions
        
            
#             if len(GOPa_gene_interactions) > 0:
                
          
                
#                 gopa_list = list(GOPa_gene_interactions['id_1']) + list(GOPa_gene_interactions['id_2'])
               
#                 # Count the occurrences of each element in the list
#                 gopa_list = Counter(gopa_list)
                
#                 # Create a DataFrame from the Counter dictionary
#                 gopa_list = pd.DataFrame(gopa_list.items(), columns=['gene_name', 'weight'])
                
#                 gopa_list = gopa_list[gopa_list['weight'] >= n_min]
    
#                 GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['id_1'].isin(list(gopa_list['gene_name']))]
#                 GOPa_gene_interactions = GOPa_gene_interactions[GOPa_gene_interactions['id_2'].isin(list(gopa_list['gene_name']))]
     
                    
#                 G = nx.Graph() 
            
             
#                 for _, row in tqdm(gopa_list.iterrows()):
#                     node = row['gene_name']
#                     color = color
#                     weight = np.log2(row['weight']*500)
#                     G.add_node(node, size = weight, color = color)
                    
#                 for index, row in tqdm(GOPa_gene_interactions.iterrows()):
#                     source = row['id_1']
#                     target = row['id_2']
#                     G.add_edge(source, target, color = 'gray')
            
                
                
#                 # Create a pyvis Network instance
#                 net = Network(notebook=True, height=f"{desired_height}px", width=f"{desired_width}px")
                
#                 net.from_nx(G)
                
            
#                 net.repulsion(node_distance=150, spring_length=200)
#                 net.show_buttons(filter_=['nodes', 'physics'])
            
                
#                 net.show(os.path.join(path, 'tmp.html'))
#                 webbrowser.open(os.path.abspath(os.path.join(path, 'tmp.html')))
            
                
#                 return G
            
#             else:
#                 print('\n \n')
#                 print('Lack of interactions!')
                
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
#     def gene_type_bar_plot(GOPa_data, color = 'gold', side = 'right', width = 10, bar_width = 0.5, count_type = 'p_val'):
        
#         """
#         This function creates a bar plot of distribution of gene types in the data
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            side (str) - orientation of bars ['left' / 'right']
#            color (str)- color of the bars
#            width (float) - width of the graph
#            bar_width (float) - width of the bars
#            count_type (str) - type of amount of term representation on bars ['perc' - percent representation / 'num' - number representation]
         
#         Returns:
#            graph: Distribution of gene types in the data
#         """
        
#         try:
    
    
#             gene_info = pd.DataFrame(GOPa_data['gene_dictionary'])
            
           
#             # Count the occurrences of each element in the list
#             gene_info = Counter(gene_info['gen_type'])
            
#             # Create a DataFrame from the Counter dictionary
#             gene_info = pd.DataFrame(gene_info.items(), columns=['gene_type', 'n'])
            
#             gene_info = gene_info.fillna('unknow')
            
#             gene_info['%'] = gene_info['n'] / sum(gene_info['n']) * 100
            
#             gene_info = gene_info.sort_values(by = 'n', ascending = False)
            
        
#             height = float(len(gene_info['gene_type'])/2.5)
            
#             fig, ax = plt.subplots(figsize=(width, height))
        
#             # Create a horizontal bar plot
#             if count_type == 'perc':
#                 ax.barh(gene_info['gene_type'], gene_info['%'], color=color, height = bar_width)
#                 ax.set_xlabel('Percent of genes [%]')
#             else:
#                 ax.barh(gene_info['gene_type'], gene_info['n'], color=color, height = bar_width)
#                 ax.set_xlabel('Number of genes')
        
#             # Set labels and title
           
#             ax.set_ylabel('')
        
#             # Invert y-axis to have names on the right
#             ax.invert_yaxis()
         
        
#             if side == 'right':
#                 ax.yaxis.tick_right()
#                 ax.set_yticks(range(len(gene_info)))
#             elif side == 'left':
#                 ax.invert_xaxis()
            
           
         
                    
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
        
        
        
#     def blood_level_markers_bar_plot(GOPa_data, side = 'right', color = 'red', experiment_type = 'MS', n_max = 10,  width = 10, standarize = False, bar_width = 0.5):
        
        
#         """
#         This function creates a bar plot of proteins level in the blood
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            side (str) - orientation of bars ['left' / 'right']
#            color (str)- color of the bars
#            experiment_type (str) - type of experiment used for protein level measurement in the blood ['MS' - mass spectometry / 'IM' - immuno]
#            n_max (int) - maximal number of bars on the graph
#            width (float) - width of the graph
#            standarize (bool) - if True; standardize Concentration via calculating log(Concentration + 1), which allows you to visualize a blood gene/protein marker level if there is a very large difference between them
#            bar_width (float) - width of the bars
         
#         Returns:
#            graph: Proteins level in the blood
#         """
    
#         try:        
            
#             blood_level = pd.DataFrame(GOPa_data['GOPa_specificity']['blood_levels'])
           
           
#             # Create a horizontal bar plot
#             if experiment_type.upper() == 'MS':
#                 blood_level = blood_level[blood_level['blood_concentration_MS[pg/L]'] == blood_level['blood_concentration_MS[pg/L]']]
#                 blood_level = blood_level[blood_level['blood_concentration_MS[pg/L]'] != 0]
    
#                 blood_level =  blood_level.sort_values(by='blood_concentration_MS[pg/L]',  ascending = False)
#                 blood_level = blood_level.reset_index(drop = True)
#                 blood_level = blood_level.iloc[0:n_max]
#                 height = float(len(blood_level['gene_name'])/2.5)
                
#                 fig, ax = plt.subplots(figsize=(width, height))
                
#                 if standarize == True:
#                     ax.barh(blood_level['gene_name'], np.log(blood_level['blood_concentration_MS[pg/L]']+1), color=color, height = bar_width)
#                     ax.set_xlabel('log(Con + 1)') 
#                 else:
#                     ax.barh(blood_level['gene_name'], blood_level['blood_concentration_MS[pg/L]'], color=color, height = bar_width)
#                     ax.set_xlabel('Concentration [pg/L]') 
#             elif experiment_type.upper() == 'IM':
#                 blood_level = blood_level[blood_level['blood_concentration_IM[pg/L]'] == blood_level['blood_concentration_IM[pg/L]']]
#                 blood_level = blood_level[blood_level['blood_concentration_IM[pg/L]'] != 0]
#                 blood_level =  blood_level.sort_values(by='blood_concentration_IM[pg/L]',  ascending = False)
#                 blood_level = blood_level.reset_index(drop = True)
#                 blood_level = blood_level.iloc[0:n_max]
                    
#                 height = float(len(blood_level['gene_name'])/2.5)
                
#                 fig, ax = plt.subplots(figsize=(width, height))
                
#                 if standarize == True:
#                     ax.barh(blood_level['gene_name'], np.log(blood_level['blood_concentration_IM[pg/L]'] + 1), color=color, height = bar_width)
#                     ax.set_xlabel('log(Con + 1)')
#                 else:
#                     ax.barh(blood_level['gene_name'], blood_level['blood_concentration_IM[pg/L]'], color=color, height = bar_width)
#                     ax.set_xlabel('Concentration [pg/L]')
    
                
          
#             # Set labels and title
           
#             ax.set_ylabel('')
        
#             # Invert y-axis to have names on the right
#             ax.invert_yaxis()
         
        
#             if side == 'right':
#                 ax.yaxis.tick_right()
#                 ax.set_yticks(range(len(blood_level)))
#             elif side == 'left':
#                 ax.invert_xaxis()
            
                    
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
        
    
#     def tissue_specificity_bar_plot(GOPa_data, p_val = 0.05, test = 'FISH', adj = 'FDR', n_max = 20, side = 'right', color = 'wheat', width = 10, bar_width = 0.5):
        
#         """
#         This function creates a bar plot of statistical / overrepresentation analysis of tissue specificity.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            p_val (float) - value of minimal p_val for statistical test
#            test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
#            adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
#            n_max (int) - maximal number of bars on the graph
#            side (str) - orientation of bars ['left' / 'right']
#            color (str)- color of the bars
#            width (float) - width of the graph
#            bar_width (float) - width of the bars
         
#         Returns:
#            graph: Bar plots of overrepresentation analysis of tissue specificity
#         """
    
#         try:
            
            
#             SEQ = pd.DataFrame(GOPa_data['GOPa_specificity']['SEQ'])
           
#             test_string = select_test(test, adj)
            
            
#             SEQ = SEQ[SEQ[test_string] <= p_val]
            
#             SEQ[test_string] = SEQ[test_string] + np.min(SEQ[test_string][SEQ[test_string] != 0])/2
    
#             SEQ['-log(p-val)'] = -np.log(SEQ[test_string])
            
            
#             SEQ = SEQ.sort_values(by = '-log(p-val)', ascending = False)
            
#             SEQ = SEQ.reset_index(drop = True)
    
           
#             SEQ = SEQ.iloc[0:n_max,:]
           
#             height = float(len(SEQ['-log(p-val)'])/2.5)
            
#             fig, ax = plt.subplots(figsize=(width, height))
            
#             ax.barh(SEQ['name'], SEQ['-log(p-val)'], color=color, height = bar_width)
#             ax.set_xlabel('-log(p-val)')
                
          
#             # Set labels and title
           
#             ax.set_ylabel('')
        
#             # Invert y-axis to have names on the right
#             ax.invert_yaxis()
         
        
#             if side == 'right':
#                 ax.yaxis.tick_right()
#                 ax.set_yticks(range(len(SEQ)))
#             elif side == 'left':
#                 ax.invert_xaxis()
            
           
         
                    
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
        
        
    
#     def cellular_role_specificity_bar_plot(GOPa_data, p_val = 0.05, test = 'FISH', adj = 'FDR', n_max = 20, side = 'right', color = 'sandybrown', width = 10, bar_width = 0.5):
        
#         """
#         This function creates a bar plot of statistical / overrepresentation analysis of tissue specificity.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
#            p_val (float) - value of minimal p_val for statistical test
#            test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
#            adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
#            n_max (int) - maximal number of bars on the graph
#            side (str) - orientation of bars ['left' / 'right']
#            color (str)- color of the bars
#            width (float) - width of the graph
#            bar_width (float) - width of the bars
         
#         Returns:
#            graph: Bar plots of overrepresentation analysis of tissue specificity
#         """
    
#         try:
            
            
#             location = pd.DataFrame(GOPa_data['GOPa_specificity']['location'])
           
#             test_string = select_test(test, adj)
            
            
#             location = location[location[test_string] <= p_val]
            
#             location[test_string] = location[test_string] + np.min(location[test_string][location[test_string] != 0])/2
            
#             location['-log(p-val)'] = -np.log(location[test_string])
            
#             location = location.sort_values(by = '-log(p-val)', ascending = False)
    
#             location = location.reset_index(drop = True)
           
#             location = location.iloc[0:n_max,:]
           
#             height = float(len(location['-log(p-val)'])/2.5)
            
#             fig, ax = plt.subplots(figsize=(width, height))
            
#             ax.barh(location['location'], location['-log(p-val)'], color=color, height = bar_width)
#             ax.set_xlabel('-log(p-val)')
                
          
#             # Set labels and title
           
#             ax.set_ylabel('')
        
#             # Invert y-axis to have names on the right
#             ax.invert_yaxis()
         
        
#             if side == 'right':
#                 ax.yaxis.tick_right()
#                 ax.set_yticks(range(len(location)))
#             elif side == 'left':
#                 ax.invert_xaxis()
            
    
                    
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
        
    
    
#     # get data functions
    
    
#     def get_gopa_results(GOPa_data):
        
#         """
#         This function gets the GOPa [GO-TERM, PATHWAYS, DISEASES and VIRAL-DISEASES] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: GOPa [GO-TERM, PATHWAYS, DISEASES and VIRAL-DISEASES] 
#         """
        
#         try:
            
#             GOPa = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa']))
#             return GOPa
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
#     def get_gene_info(GOPa_data):
        
#         """
#         This function gets the genes info from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: Genes info
#         """
        
#         try:
            
#             gene_info = pd.DataFrame(copy.deepcopy(GOPa_data['gene_dictionary']))
#             return gene_info
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
#     def get_gopa_interactions_results(GOPa_data):
        
#         """
#         This function gets the GOPa interactions [GO-TERM, PATHWAYS, DISEASES and VIRAL-DISEASES] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: GOPa interactions [GO-TERM, PATHWAYS, DISEASES and VIRAL-DISEASES] 
#         """
        
#         try:
            
#             GOPa = pd.DataFrame(GOPa_data['GOPa'])
#             GOPa_interactions = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa_interactions']))
#             GOPa = GOPa[['GOPa', 'relation_id']]
#             GOPa = GOPa.explode('relation_id')
#             GOPa = GOPa.drop_duplicates()
#             name_mapping = dict(zip(GOPa['relation_id'], GOPa['GOPa']))
#             GOPa_interactions['A_name'] = GOPa_interactions['A'].map(name_mapping)
#             GOPa_interactions['B_name'] = GOPa_interactions['B'].map(name_mapping)
#             GOPa_interactions = GOPa_interactions.drop('color', axis = 1)
#             del GOPa
        
#             return GOPa_interactions
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def get_gopa_gene_interaction(GOPa_data):
        
#         """
#         This function gets the gene interactions [IntAct, STRING] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: Gene interactions [IntAct, STRING]
#         """
        
#         try:
            
#             GOPa = pd.DataFrame(GOPa_data['gene_dictionary'])
#             GOPa_gene_interaction = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa_gene_interaction']))
#             GOPa = GOPa[['gene_name', 'dictionary_id']]
        
#             name_mapping = dict(zip(GOPa['dictionary_id'], GOPa['gene_name']))
#             GOPa_gene_interaction['id1_name'] = GOPa_gene_interaction['id_1'].map(name_mapping)
#             GOPa_gene_interaction['id2_name'] = GOPa_gene_interaction['id_2'].map(name_mapping)
#             del GOPa
        
#             return GOPa_gene_interaction
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def get_gopa_blood_markers(GOPa_data):
    
#         """
#         This function gets the blood markers [HPA] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: Blood markers [HPA]
#         """
        
#         try:
            
#             blood_levels = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa_specificity']['blood_levels']))
        
#             return blood_levels
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
#     def get_gopa_cellular_role_specificity(GOPa_data):
        
#         """
#         This function gets the cellular location specificity [HPA] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: Cellular location specificity [HPA]
#         """
        
#         try:
            
#             cellular_specificity = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa_specificity']['location']))
        
#             return cellular_specificity
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
#     def get_gopa_tissue_specificity(GOPa_data):
        
#         """
#         This function gets the tissue specificity [HPA] from GOPa_data dictionary and return in data frame.
    
#         Args:
#            GOPa_data (dict) - raw GOPa_data from gopa_interaction_analysis function
           
         
#         Returns:
#            data_frame: Tissue specificity [HPA]
#         """
        
#         try:
            
#             tissue_specificity = pd.DataFrame(copy.deepcopy(GOPa_data['GOPa_specificity']['SEQ']))
        
#             return tissue_specificity
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
    
#     #rnaseq functions
    
#     #visualization scatter expression
    
#     def gene_scatter(data, colors = 'viridis', species = 'human', hclust = 'complete', img_width = None, img_high = None, label_size = None, x_lab = 'Genes', legend_lab = 'log(TPM + 1)'):
    
            
#         """
#         This function creates a graph in the format of a scatter plot for expression data prepared in data frame format.
        
#         Args:
#            data (data frame) - data frame of genes/protein expression where on row are the gene/protein names and on column grouping variable (tissue / cell / ect. names)
#            color (str) - palette color available for matplotlib in python eg. viridis
#            species (str) - species for upper() or lower() letter for gene/protein name depending on 
#            hclust (str) - type of data clustering of input expression data eg. complete or None if  no clustering
#            img_width (float) - width of the image or None for auto-adjusting
#            img_high (float) - high of the image or None for auto-adjusting
#            label_size (float) - labels size of the image or None for auto-adjusting
#            x_lab (str) - tex for x axis label
#            legend_lab (str) - description for legend label
           
           
#         Returns:
#            graph: Scatter plot of expression data
#         """
        
        
#         try:
#             scatter_df = data        
    
#             if img_width == None:
#                 img_width = len(scatter_df.columns)*1.2
            
#             if img_high == None:
#                 img_high = len(scatter_df.index)*0.9
                
            
#             if label_size == None:
#                 label_size = np.log(len(scatter_df.index)  *  len(scatter_df.index))*2.5
                
#                 if label_size < 7:
#                     label_size = 7
            
#             cm = 1/2.54
            
#             if len(scatter_df) > 1:
                
               
                
            
#                 Z = linkage(scatter_df, method=hclust)
        
        
#                 # Get the order of features based on the dendrogram
#                 order_of_features = dendrogram(Z, no_plot=True)['leaves']
        
#                 indexes_sort = list(scatter_df.index)
#                 sorted_list_rows = []
#                 for n in order_of_features:
#                     sorted_list_rows.append(indexes_sort[n])
                    
                
                
#                 scatter_df = scatter_df.transpose()
            
#                 Z = linkage(scatter_df, method=hclust)
        
#                 # Get the order of features based on the dendrogram
#                 order_of_features = dendrogram(Z, no_plot=True)['leaves']
        
#                 indexes_sort = list(scatter_df.index)
#                 sorted_list_columns = []
#                 for n in order_of_features:
#                     sorted_list_columns.append(indexes_sort[n])
                            
                       
#                 scatter_df = scatter_df.transpose()
                
#                 scatter_df = scatter_df.loc[sorted_list_rows, sorted_list_columns]
                 
#             scatter_df = np.log(scatter_df + 1)
#             scatter_df[scatter_df <= np.mean(scatter_df.quantile(0.10))] = np.mean(np.mean(scatter_df, axis=1))/10
        
#             if species.lower() == 'human':
#                 scatter_df.index = [x.upper() for x in scatter_df.index ]
#             else:
#                 scatter_df.index  = [x.title() for x in scatter_df.index ]
                
#             scatter_df.insert(0, '  ', 0)
    
#             # Add a column of zeros at the end
#             scatter_df[' '] = 0
                 
#             fig, ax = plt.subplots(figsize=(img_width*cm,img_high*cm))
        
#             plt.scatter(x = [*range(0, len(scatter_df.columns), 1)], y = [' '] * len(scatter_df.columns),s=0, cmap=colors,  edgecolors=None)
        
            
    
        
#             for index, row in enumerate(scatter_df.index):
#                 x = [*range(0, len(np.array(scatter_df.loc[row,])), 1)]
#                 y = [row] * len(x)
#                 s = np.array(scatter_df.loc[row,])
#                 plt.scatter(x,y,s=np.log(s+1)*70, c=s, cmap=colors,  edgecolors='black', vmin = np.array(scatter_df).min() ,vmax = np.array(scatter_df).max(), linewidth=0.00001)
#                 sm = plt.cm.ScalarMappable(cmap=colors)
#                 sm.set_clim(vmin = np.array(scatter_df).min() ,vmax = np.array(scatter_df).max())
#                 plt.xticks(x, scatter_df.columns)
#                 plt.ylabel(str(x_lab), fontsize=label_size)
                
                
#             plt.scatter(x = [*range(0, len(scatter_df.columns), 1)], y = [''] * len(scatter_df.columns),s=0, cmap=colors,  edgecolors=None)
    
    
            
        
#             plt.xticks(rotation = 80) 
#             plt.tight_layout()
#             plt.margins(0.005)
#             plt.xticks(fontsize=label_size)
#             plt.yticks(fontsize=label_size)
     
        
        
#             len_bar = ax.get_position().height/5
#             if len(scatter_df) < 15:
#                 len_bar = 0.65
                
#                 cbar = plt.colorbar(sm)
#                 cbar.ax.set_ylabel(str(legend_lab), fontsize=label_size*0.9)
#                 cbar.ax.yaxis.set_ticks_position('right')
#                 cbar.ax.set_position([ax.get_position().x1 + 0.05, (ax.get_position().y0 + ax.get_position().y1)/1.9 , ax.get_position().width/0.05, len_bar])
#                 cbar.ax.yaxis.set_label_position('right')
#                 cbar.ax.yaxis.set_tick_params(labelsize=label_size*0.8)
#                 cbar.outline.set_edgecolor('none')
#             else:
#                 cbar = plt.colorbar(sm)
#                 cbar.ax.set_ylabel(str(legend_lab), fontsize=label_size*0.9)
#                 cbar.ax.yaxis.set_ticks_position('right')
#                 cbar.ax.set_position([ax.get_position().x1 + 0.05, (ax.get_position().y0 + ax.get_position().y1)/1.45 , ax.get_position().width/0.05, len_bar])
#                 cbar.ax.yaxis.set_label_position('right')
#                 cbar.ax.yaxis.set_tick_params(labelsize=label_size*0.8)
#                 cbar.outline.set_edgecolor('none')
        
#             ax.spines['top'].set_visible(False)
#             ax.spines['right'].set_visible(False)
#             ax.spines['bottom'].set_visible(False)
#             ax.spines['left'].set_visible(False)     
#             ax.xaxis.set_tick_params(length=0,labelbottom=True)
#             ax.yaxis.set_tick_params(length=0,labelbottom=True)
#             ax.grid(False)
        
        
#             return fig
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
#     #preapre rnaseq
    
#     def gene_specificity(gene_list, path_in_use = _path_in_inside):
        
#         """
#         This function creates a dictionary of RNA-SEQ data for studied genes.
        
#         Args:
#            gene_list (list) - list of genes to check in the context of rna-seq data 
    
           
#         Returns:
#            dictionary: Dictionary of RNA-SEQ data for studied genes
#         """
        
#         try:
            
#             print('Genes search start...')
            
#             with open(path_in_use + '/GOPa_metadata_dict.json', 'r') as json_file:
#                 GOPa_metadata = (json.load(json_file))
                
                
#             GOPa_data, not_found = search_genes(gene_list, GOPa_metadata, species='human')
            
#             GOPa_genes = list(GOPa_data['gene_dictionary']['gene_name'])
            
#             del GOPa_metadata, GOPa_data
        
#             print('RNA-SEQ data searching...')
                
#             with open(path_in_use + '/human_tissue_expression_fetal_development_circular.json', 'r') as json_file:
#                 human_tissue_expression_fetal_development_circular = (json.load(json_file))
                
#             human_tissue_expression_fetal_development_circular = pd.DataFrame.from_dict(dict({k:v for k,v in human_tissue_expression_fetal_development_circular.items() if k != 'tissue'}), orient='index',  columns = human_tissue_expression_fetal_development_circular['tissue'])
#             human_tissue_expression_fetal_development_circular = human_tissue_expression_fetal_development_circular.loc[[x for x in GOPa_genes if x in human_tissue_expression_fetal_development_circular.index], :]
        
#             with open(path_in_use + '/human_tissue_expression_HPA.json', 'r') as json_file:
#                 human_tissue_expression_HPA = (json.load(json_file))
                
#             human_tissue_expression_HPA = pd.DataFrame.from_dict(dict({k:v for k,v in human_tissue_expression_HPA.items() if k != 'tissue'}), orient='index',  columns = human_tissue_expression_HPA['tissue'])
#             human_tissue_expression_HPA = human_tissue_expression_HPA.loc[[x for x in GOPa_genes if x in human_tissue_expression_HPA.index], :]
            
            
#             with open(path_in_use + '/human_tissue_expression_illumina_bodyMap2.json', 'r') as json_file:
#                 human_tissue_expression_illumina_bodyMap2 = (json.load(json_file))
                
#             human_tissue_expression_illumina_bodyMap2 = pd.DataFrame.from_dict(dict({k:v for k,v in human_tissue_expression_illumina_bodyMap2.items() if k != 'tissue'}), orient='index',  columns = human_tissue_expression_illumina_bodyMap2['tissue'])
#             human_tissue_expression_illumina_bodyMap2 = human_tissue_expression_illumina_bodyMap2.loc[[x for x in GOPa_genes if x in human_tissue_expression_illumina_bodyMap2.index], :]
        
        
#             with open(path_in_use + '/human_tissue_expression_RNA_total_tissue.json', 'r') as json_file:
#                 human_tissue_expression_RNA_total_tissue = (json.load(json_file))
                
#             human_tissue_expression_RNA_total_tissue = pd.DataFrame.from_dict(dict({k:v for k,v in human_tissue_expression_RNA_total_tissue.items() if k != 'tissue'}), orient='index',  columns = human_tissue_expression_RNA_total_tissue['tissue'])
#             human_tissue_expression_RNA_total_tissue = human_tissue_expression_RNA_total_tissue.loc[[x for x in GOPa_genes if x in human_tissue_expression_RNA_total_tissue.index], :]
        
#             seq_dict = {'human_tissue_expression_fetal_development_circular':human_tissue_expression_fetal_development_circular,
#                         'human_tissue_expression_HPA':human_tissue_expression_HPA,
#                         'human_tissue_expression_illumina_bodyMap2':human_tissue_expression_illumina_bodyMap2,
#                         'human_tissue_expression_RNA_total_tissue':human_tissue_expression_RNA_total_tissue}
            
#             return seq_dict
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
    
#     def rna_seq_scatter(rna_seq):
        
          
#         """
#         This function creates a ditionary of graphs for RNA-SEQ data from gene_specificity function
        
#         Args:
#            rna_seq (dictionary) - data from gene_specificity function
           
#         Returns:
#            dictionary: Dictionary of graphs for different type of RNA-SEQ data for different tissue specificity
#         """
        
#         try:
            
#             image_dict = {}
#             for k in rna_seq.keys():
#                 print('Preparing graph for ' + str(k))
                
                
#                 image_dict[k] =gene_scatter(rna_seq[k], colors = 'viridis', species = 'human', hclust = 'complete', img_width = None, img_high = None, label_size = None)
        
#             return image_dict
        
#         except:
#             print('\n')
#             print("Something went wrong. Check the function input data and try again!")
    
    
    
#     #Deregulated GOPa function - compare two genes / proteins list 
    
    
#     def DGOPa(gene_up, gene_down, GOPa_metadata, species = None, min_fc = 1.5, p_val = 0.05, test = 'FISH', adj = 'FDR'):
        
#         """
#         This function conducts full GOPa analysis on two gene / protein lists [gene_up, gene_down] including search_genes, gopa_analysis, gopa_interaction_analysis, and gopa_specificity_analysis.
#         This function is useful when you have to compare two sides genes / proteins deregulations (upregulated vs. downregulated) and choose the most accurate results adjusted to the provided gene lists based on min_fc (minimal fold change) between GOPa projects.
    
#         Args:
#            gene_up (list)- list of genes eg. ['KIT', 'EDNRB', 'PAX3'] 
#            gene_down (list)- list of genes eg. ['KIT', 'EDNRB', 'PAX3'] 
#            GOPa_metadata (dict) - metadata from load_GOPa_meta function 
#            species (str or None) - ['human' / 'mouse' / 'both' / None] 
#            min_fc (float) - minimal value of fold change the normalized by number of genes in analysis results of GOPa between GOPa projects obtained from the provided lists of genes [gene_up / gene_down]
#            p_val (float) - value of minimal p_val for statistical test
#            test (str) - type of statistical test ['FISH' - Fisher's exact test / 'BIN' - Binomial test]
#            adj (str) - type of p_value correction ['BF' - Bonferroni correction / 'FDR' - False Discovery Rate (BH procedure)]
           
#            If choose 'human' or 'mouse' you will obtain information about this species' genes. 
#            If choose 'both' you will obtain information for genes that are available mutually for both species. 
#            If choose None you will obtain information for all genes available in the metadata.       
    
#         Returns:
#            list of dict: A list of two dictionaries with analyzed GOPa projects obtained on the provided list of genes [gene_up, gene_down] corrected on specific occurrences in GOPa results with min_fc. The first dictionary is related to gene_up results and the second to gene_down results 
           
#         """
        
#         test_string = select_test(test, adj)
    
    
#         GOPa_data_up, not_found = search_genes(gene_up, GOPa_metadata, species=species)
         
        
#         GOPa_data_up =  gopa_analysis(GOPa_data_up, GOPa_metadata)
        
          
#         GOPa_data_up =  gopa_interaction_analysis(GOPa_data_up)
        
#         GOPa_data_up =  gopa_specificity_analysis(GOPa_data_up, GOPa_metadata)
    
    
        
#         GOPa_data_down, not_found = search_genes(gene_down, GOPa_metadata, species=species)
         
        
#         GOPa_data_down =  gopa_analysis(GOPa_data_down, GOPa_metadata)
        
          
#         GOPa_data_down =  gopa_interaction_analysis(GOPa_data_down)
        
        
#         GOPa_data_down =  gopa_specificity_analysis(GOPa_data_down, GOPa_metadata)
    
        
#         #GOPa
    
#         GOPa_up = pd.DataFrame(GOPa_data_up['GOPa'])
#         GOPa_up = GOPa_up[GOPa_up[test_string] <= p_val]
#         GOPa_up['norm_n'] =  GOPa_up['n'] / len(gene_up)
        
        
#         GOPa_down = pd.DataFrame(GOPa_data_down['GOPa'])
#         GOPa_down = GOPa_down[GOPa_down[test_string] <= p_val]
#         GOPa_down['norm_n'] =  GOPa_down['n'] / len(gene_down)
        
        
#         gopa_down_list = []
#         for g in GOPa_down['GOPa']:
#             if g in list(GOPa_up['GOPa']):
#                 if(float(GOPa_down['norm_n'][GOPa_down['GOPa'] == g])/ float(GOPa_up['norm_n'][GOPa_up['GOPa'] == g]) >= min_fc):
#                     gopa_down_list.append(g)
#             else:
#                 gopa_down_list.append(g)
            
                
#         gopa_up_list = []
#         for g in GOPa_up['GOPa']:
#             if g in list(GOPa_down['GOPa']):
#                 if(float(GOPa_up['norm_n'][GOPa_up['GOPa'] == g])/ float(GOPa_down['norm_n'][GOPa_down['GOPa'] == g]) >= min_fc):
#                     gopa_up_list.append(g)
#             else:
#                 gopa_up_list.append(g)
                
                
#         GOPa_down = GOPa_down[GOPa_down['GOPa'].isin(gopa_down_list)]
#         GOPa_data_down['GOPa'] = GOPa_down.to_dict(orient = 'list')
        
#         GOPa_up = GOPa_up[GOPa_up['GOPa'].isin(gopa_up_list)]
#         GOPa_data_up['GOPa'] = GOPa_up.to_dict(orient = 'list')
        
#         #specificity
        
        
#         GOPa_up = pd.DataFrame(GOPa_data_up['GOPa_specificity']['SEQ'])
#         GOPa_up = GOPa_up[GOPa_up[test_string] <= p_val]
#         GOPa_up['norm_n'] =  GOPa_up['n'] / len(gene_up)
        
        
#         GOPa_down = pd.DataFrame(GOPa_data_down['GOPa_specificity']['SEQ'])
#         GOPa_down = GOPa_down[GOPa_down[test_string] <= p_val]
#         GOPa_down['norm_n'] =  GOPa_down['n'] / len(gene_down)
        
        
#         gopa_down_list = []
#         for g in GOPa_down['name']:
#             if g in list(GOPa_up['name']):
#                 if(float(GOPa_down['norm_n'][GOPa_down['name'] == g])/ float(GOPa_up['norm_n'][GOPa_up['name'] == g]) >= min_fc):
#                     gopa_down_list.append(g)
#             else:
#                 gopa_down_list.append(g)
            
#         gopa_up_list = []
#         for g in GOPa_up['name']:
#             if g in list(GOPa_down['name']):
#                 if(float(GOPa_up['norm_n'][GOPa_up['name'] == g])/ float(GOPa_down['norm_n'][GOPa_down['name'] == g]) >= min_fc):
#                     gopa_up_list.append(g)
#             else:
#                 gopa_up_list.append(g)
                
                
#         GOPa_down = GOPa_down[GOPa_down['name'].isin(gopa_down_list)]
#         GOPa_data_down['GOPa_specificity']['SEQ'] = GOPa_down.to_dict(orient = 'list')
        
#         GOPa_up = GOPa_up[GOPa_up['name'].isin(gopa_up_list)]
#         GOPa_data_up['GOPa_specificity']['SEQ'] = GOPa_up.to_dict(orient = 'list')
        
        
#         #location
        
        
#         GOPa_up = pd.DataFrame(GOPa_data_up['GOPa_specificity']['location'])
#         GOPa_up = GOPa_up[GOPa_up[test_string] <= p_val]
#         GOPa_up['norm_n'] =  GOPa_up['n'] / len(gene_up)
        
        
#         GOPa_down = pd.DataFrame(GOPa_data_down['GOPa_specificity']['location'])
#         GOPa_down = GOPa_down[GOPa_down[test_string] <= p_val]
#         GOPa_down['norm_n'] =  GOPa_down['n'] / len(gene_down)
        
        
        
#         gopa_down_list = []
#         for g in GOPa_down['location']:
#             if g in list(GOPa_up['location']):
#                 if(float(GOPa_down['norm_n'][GOPa_down['location'] == g])/ float(GOPa_up['norm_n'][GOPa_up['location'] == g]) >= min_fc):
#                     gopa_down_list.append(g)
#             else:
#                 gopa_down_list.append(g)
          
                
#         gopa_up_list = []
#         for g in GOPa_up['location']:
#             if g in list(GOPa_down['location']):
#                 if(float(GOPa_up['norm_n'][GOPa_up['location'] == g])/ float(GOPa_down['norm_n'][GOPa_down['location'] == g]) >= min_fc):
#                     gopa_up_list.append(g)
#             else:
#                 gopa_up_list.append(g)
                
                
#         GOPa_down = GOPa_down[GOPa_down['location'].isin(gopa_down_list)]
#         GOPa_data_down['GOPa_specificity']['location'] = GOPa_down.to_dict(orient = 'list')
        
#         GOPa_up = GOPa_up[GOPa_up['location'].isin(gopa_up_list)]
#         GOPa_data_up['GOPa_specificity']['location'] = GOPa_up.to_dict(orient = 'list')
        
        
#         return GOPa_data_up, GOPa_data_down
