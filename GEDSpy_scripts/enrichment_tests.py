from enrichment import GetDataRaw

gdr = GetDataRaw()


recteome = gdr.get_raw_REACTOME()

"""
This method gets the REACTOME data.

Returns:
   dict: REACTOME data
"""


ref_gen = gdr.get_raw_REF_GEN()


"""
This method gets the REF_GEN which is the combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes for scientific use.
 
Returns:
   dict: Combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes
"""


ref_gen_seq = gdr.get_raw_REF_GEN_RNA_SEQ()

"""
This method gets the tissue-specific RNA-SEQ data including:
    -human_tissue_expression_HPA
    -human_tissue_expression_RNA_total_tissue
    -human_tissue_expression_fetal_development_circular
    -human_tissue_expression_illumina_bodyMap2


Returns:
   dict: Tissue specific RNA-SEQ data
"""



HPA = gdr.get_raw_HPA()

"""
This method gets the HPA (Human Protein Atlas) data.
 
Returns:
   dict: HPA data
"""



diseases = gdr.get_raw_DISEASES()

"""
This method gets the DISEASES data.

Returns:
   dict: DISEASES data
"""



vimic = gdr.get_raw_ViMIC()

"""
This method gets the ViMIC data.
 
Returns:
   dict: ViMIC data
"""



kegg = gdr.get_raw_KEGG()

"""
This method gets the KEGG data.

Returns:
   dict: KEGG data
"""



go = gdr.get_raw_GO()

"""
This method gets the GO-TERM data.
 
Returns:
   data_frame: GO-TERM data
"""



intact = gdr.get_raw_IntAct()

"""
This method gets the IntAct data.

Returns:
   dict: IntAct data
"""
 


string = gdr.get_raw_STRING()

"""
This method gets the STRING data.
 
Returns:
   dict: STRING data
"""



cell_talk = gdr.get_raw_CellTalk()

"""
This method gets the CellTalk data.
 
Returns:
   dict: CellTalk data
"""




cell_phone = gdr.get_raw_CellPhone()

"""
This method gets the CellPhone data.
 
Returns:
   dict: CellPhone data
"""



from enrichment import GetData

gd = GetData()


reactome = gd.get_REACTOME()
     
"""
This method gets the REACTOME data including the id to connect with REF_GENE by id_reactome

Returns:
   data_frame: REACTOME data
"""
     
    


ref_gen = gd.get_REF_GEN()
     
"""
This method gets the REF_GEN which is the combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes for scientific use.
 
Returns:
   dict: Combination of Homo sapiens / Mus musculus / Rattus norvegicus genomes
"""
     
    

ref_gen_seq = gd.get_REF_GEN_RNA_SEQ()
     
"""
This method gets the tissue-specific RNA-SEQ data including:
    -human_tissue_expression_HPA
    -human_tissue_expression_RNA_total_tissue
    -human_tissue_expression_fetal_development_circular
    -human_tissue_expression_illumina_bodyMap2

  
Returns:
   dict: Tissue specific RNA-SEQ data
"""
     
   
         
         
HPA = gd.get_HPA()
     
"""
This method gets the HPA (Human Protein Atlas) data including the id to connect with REF_GENE by id_HPA
 
Returns:
   dict: HPA data
"""

    


disease = gd.get_DISEASES()
     
"""
This method gets the DISEASES data including the id to connect with REF_GENE by id_diseases

Returns:
   dict: DISEASES data
"""
     
    


vimic = gd.get_ViMIC()
     
"""
This method gets the ViMIC data including the id to connect with REF_GENE by id_viral_diseases
 
Returns:
   dict: ViMIC data
"""

    


kegg = gd.get_KEGG()
     
"""
This method gets the KEGG data including the id to connect with REF_GENE by id_kegg

Returns:
   dict: KEGG data
"""
     
    


go = gd.get_GO()
     
"""
This method gets the GO-TERM data including the id to connect with REF_GENE by id_go
 
Returns:
   data_frame: GO-TERM data
"""
   



intact = gd.get_IntAct()
     
"""
This method gets the IntAct data including the id to connect with REF_GENE by id_IntAct

Returns:
   data_frame: IntAct data
"""
     
     


string = gd.get_STRING()
     
"""
This method gets the STRING data including the id to connect with REF_GENE by id_string
 
Returns:
   data_frame: STRING data
"""




CellTalk = gd.get_CellTalk()

"""
This method gets the CellTalk data after adjustment.
 
Returns:
   dict: CellTalk data
"""




CellPhone = gd.get_CellPhone()

"""
This method gets the CellPhone data after adjustment.
 
Returns:
   dict: CellPhone data
"""




CellInteractions = gd.get_interactions()

"""
This method gets the CellPhone & CellTalk data including the id to connect with REF_GENE by id_go.
 
Returns:
   dict: CellPhone data
"""





# ENRICHMENT

from enrichment import Enrichment

enr = Enrichment()


enr.select_features(['KIT', 'MAPT', 'PAX3', 'CALB1', 'ATXN3', 'fd'])


enr.full_enrichment()


HPA = enr.get_HPA
      
    
String = enr.get_STRING
      
    
go = enr.get_GO_TERM
      
  
reactome = enr.get_REACTOME
      
    
   
cc = enr.get_CellCon
      
     
   
intact =  enr.get_IntAct
    

kegg = enr.get_KEGG

      
dis = enr.get_DISEASES
      
    
vimic = enr.get_ViMIC

       
rna_seq = enr.get_RNA_SEQ   


results = enr.get_results
