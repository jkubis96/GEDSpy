def enrichment():
    
    

from enrichment import GetDataRaw


gdr = GetDataRaw()


recteome = gdr.get_raw_REACTOME()


ref_gen = gdr.get_raw_REF_GEN()


ref_gen_seq = gdr.get_raw_REF_GEN_RNA_SEQ()


HPA = gdr.get_raw_HPA()


diseases = gdr.get_raw_DISEASES()


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


from enrichment import GetDataRaw

gdr = GetDataRaw()




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



gene_list = ['CACNA1I','CALD1','CAMK1G','CAMK2N1','CAMSAP1','CCL15','CCL16','CCNL2','CCT8P1','CD46','CDC14A','CDK18','CDK19','CES3','CHEK2',
			 'CHID1','COL6A3','CPVL','CYP3A43','CYP3A5','DBNL','DUSP10','DUSP9','ECHS1','EGFR','EGR2','ELL2','ERMP1','ESR1','F7','FAM171A1',
			 'FAM20C','FGFR2','FH','FLAD1','FUT3','GAA','GBA2','GGCX','GJB1','GLRX5','GNAI2','GNB2','GNB3','GPNMB','GRB10','GRHPR','HMGCS2',
			 'HSD17B4','HSP90AB1','IER3IP1','IGF2R','IL1R1','INF2','IRAK1','ITGA1','ITGA7','ITIH1','ITIH3','ITIH4','ITPR1','ITSN1','JAK1',
			 'KALRN','KCNQ2','KCNQ4','KDM3A','KMO','KRAS','KSR1','LAMA5','LAMB2','LCN2','MAP2K7','MAP4K2','MAP4K3',
			 'MAPK13','MARCO','MAST2','MAT1A','MATR3','MCM8','MFSD10','MGAT5','MTMR10','MUSK','MYO9B','NBAS','NCOA6','NCSTN','NDUFA4','NEK4',
			 'NPR2','NUDT2','NUP210','ORC3L','PAOX','PEMT','PEX14','PFKL','PHKA2','PIM1','PLXND1','PMM1','PON3','POR','PPARG','PPARGC1B',
			 'PPP2R1A','PRKCE','PTK2B','PTP4A1','PTPN23','PTPRF','PTPRK','RARA','RNF10','RNF14','RNF165','ROCK2','RRBP1','RREB1','SCN1A','SDC1',
			 'SEPHS1','SERPINA1','SERPINA10','SFXN5','SHROOM1','SIL1','SIRPA','SLC12A7','SLC13A3','SLC16A2','SLC17A7','SLC22A23','SLC22A9',
			 'SLC23A2','SLC25A11','SLC25A25','SLC38A3','SLC45A3','SLC4A5','SLC5A1','SLC7A2','SLC8A3','SLC9A6','SLCO1A2','SLCO1B3','SMARCA2',
			 'SNRK','SNX4','SORBS1','SPEN','SPR','SRF','STAB1','STAT1','SUCLG2','SULT1B1','SULT1E1','TBC1D2B','TCHP','TGFBI','TGOLN2','THPO',
			 'TIE1','TIMM13','TLK2','TMEM62','TNFSF14','TNK2','TNS1','TPI1','TRIB3','TRMT11','TTYH3']




# Split the list into two halves
half = len(gene_list) // 2
list1 = gene_list[:half]
list2 = gene_list[half:]



# ENRICHMENT

from enrichment import Enrichment

enr = Enrichment()


enr.select_features(list1)


enr.full_enrichment()


# HPA = enr.get_HPA
      
    
# String = enr.get_STRING
      
    
# go = enr.get_GO_TERM
      
  
# reactome = enr.get_REACTOME
      
    
   
# cc = enr.get_CellCon
      
     
   
# intact =  enr.get_IntAct
    

# kegg = enr.get_KEGG

      
# dis = enr.get_DISEASES
      
    
# vimic = enr.get_ViMIC

       
# rna_seq = enr.get_RNA_SEQ   


results = enr.get_results





from enrichment import Analysis


ans = Analysis(results)

ans.networks_metadata

ans.interactions_metadata


ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    
 
ans.networks_metadata


ans.full_analysis()


results2 = ans.get_full_results


import json

with open('tmp.json', 'w') as json_file:
    json.dump(results2, json_file)


   
import json

with open('tmp.json', 'r') as json_file:
    results2 = (json.load(json_file))



from JVG import JVG

from enrichment import Visualization


vis = Visualization(results2)


plot  =  vis.GO_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 20, 
                    side = 'left', 
                    color = 'blue', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()




plot  =  vis.SPECIFICITY_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 25, 
                    side = 'right', 
                    color = 'bisque', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()





plot  =  vis.KEGG_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 25, 
                    side = 'right', 
                    color = 'orange', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()




plot  =  vis.REACTOME_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 25, 
                    side = 'right', 
                    color = 'silver', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()




plot  =  vis.DISEASES_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 25, 
                    side = 'right', 
                    color = 'thistle', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()



plot  =  vis.ViMIC_plot(p_val = 0.05, 
                    test = 'FISH', 
                    adj = None, 
                    n = 25, 
                    side = 'right', 
                    color = 'aquamarine', 
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val')


graph = JVG.MplEditor(plot)

graph.edit()






plot  =  vis.blod_markers_plot()


graph = JVG.MplEditor(plot)

graph.edit()





go = vis.GOPa_network_create(
                        data_set = 'GO-TERM', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = True,
                        include_childrend = True,
                        selected_parents = [],
                        selected_genes = [])



nt = JVG.NxEditor(go)
nt.edit()




kegg = vis.GOPa_network_create(
                        data_set = 'KEGG', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = True,
                        include_childrend = True,
                        selected_parents = ['Nervous system', 'Cell motility'],
                        selected_genes = [])



nt = JVG.NxEditor(kegg)
nt.edit()




react = vis.GOPa_network_create(
                        data_set = 'REACTOME', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = False,
                        include_childrend = False,
                        selected_parents = [],
                        selected_genes = ['F7', 'LAMA5'])





nt = JVG.NxEditor(react)
nt.edit()





inter = vis.GI_network_create()


nt = JVG.NxEditor(inter)
nt.edit()



ml = vis.AUTO_ML_network( 
                    genes_inc = 10, 
                    gene_int = True, 
                    genes_only = True, 
                    min_con = 1, 
                    children_con = True, 
                    include_childrend = True,
                    selected_parents = ['Cell motility', 'Development and regeneration'],
                    selected_genes = ['JAK1', 'KRAS'])



nt = JVG.NxEditor(ml)
nt.edit()




seq = vis.gene_scatter( 
                 colors = 'viridis', 
                 species = 'human', 
                 hclust = 'complete', 
                 img_width = None, 
                 img_high = None, 
                 label_size = None, 
                 x_lab = 'Genes', 
                 legend_lab = 'log(TPM + 1)')



seq.keys()

graph = JVG.MplEditor(seq['tissue_expression_RNA_total_tissue'])

graph.edit()




###################################################################################



# DSA




gene_list = ['CACNA1I','CALD1','CAMK1G','CAMK2N1','CAMSAP1','CCL15','CCL16','CCNL2','CCT8P1','CD46','CDC14A','CDK18','CDK19','CES3','CHEK2',
			 'CHID1','COL6A3','CPVL','CYP3A43','CYP3A5','DBNL','DUSP10','DUSP9','ECHS1','EGFR','EGR2','ELL2','ERMP1','ESR1','F7','FAM171A1',
			 'FAM20C','FGFR2','FH','FLAD1','FUT3','GAA','GBA2','GGCX','GJB1','GLRX5','GNAI2','GNB2','GNB3','GPNMB','GRB10','GRHPR','HMGCS2',
			 'HSD17B4','HSP90AB1','IER3IP1','IGF2R','IL1R1','INF2','IRAK1','ITGA1','ITGA7','ITIH1','ITIH3','ITIH4','ITPR1','ITSN1','JAK1',
			 'KALRN','KCNQ2','KCNQ4','KDM3A','KMO','KRAS','KSR1','LAMA5','LAMB2','LCN2','MAP2K7','MAP4K2','MAP4K3',
			 'MAPK13','MARCO','MAST2','MAT1A','MATR3','MCM8','MFSD10','MGAT5','MTMR10','MUSK','MYO9B','NBAS','NCOA6','NCSTN','NDUFA4','NEK4',
			 'NPR2','NUDT2','NUP210','ORC3L','PAOX','PEMT','PEX14','PFKL','PHKA2','PIM1','PLXND1','PMM1','PON3','POR','PPARG','PPARGC1B',
			 'PPP2R1A','PRKCE','PTK2B','PTP4A1','PTPN23','PTPRF','PTPRK','RARA','RNF10','RNF14','RNF165','ROCK2','RRBP1','RREB1','SCN1A','SDC1',
			 'SEPHS1','SERPINA1','SERPINA10','SFXN5','SHROOM1','SIL1','SIRPA','SLC12A7','SLC13A3','SLC16A2','SLC17A7','SLC22A23','SLC22A9',
			 'SLC23A2','SLC25A11','SLC25A25','SLC38A3','SLC45A3','SLC4A5','SLC5A1','SLC7A2','SLC8A3','SLC9A6','SLCO1A2','SLCO1B3','SMARCA2',
			 'SNRK','SNX4','SORBS1','SPEN','SPR','SRF','STAB1','STAT1','SUCLG2','SULT1B1','SULT1E1','TBC1D2B','TCHP','TGFBI','TGOLN2','THPO',
			 'TIE1','TIMM13','TLK2','TMEM62','TNFSF14','TNK2','TNS1','TPI1','TRIB3','TRMT11','TTYH3']




# Split the list into two halves
half = len(gene_list) // 2
list1 = gene_list[:half]
list2 = gene_list[half:]






# ENRICHMENT

from enrichment import Enrichment
from enrichment import Analysis



enr = Enrichment()


enr.select_features(list1)


enr.full_enrichment()


results1 = enr.get_results


ans = Analysis(results1)

ans.networks_metadata

ans.interactions_metadata


ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    
 


ans.full_analysis()


results1 = ans.get_full_results





import json

with open('set1.json', 'w') as json_file:
    json.dump(results1, json_file)


   




enr = Enrichment()


enr.select_features(list2)


enr.full_enrichment()


results2 = enr.get_results



ans = Analysis(results2)

ans.networks_metadata

ans.interactions_metadata


ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    
 


ans.full_analysis()


results2 = ans.get_full_results



import json

with open('set2.json', 'w') as json_file:
    json.dump(results2, json_file)


   








import json

with open('set1.json', 'r') as json_file:
    set1 = (json.load(json_file))


import json

with open('set2.json', 'r') as json_file:
    set2 = (json.load(json_file))



from enrichment import DSA


d = DSA(set1,set2)
 

# d.GO_diff()

# t = d.get_GO_diff

# d.KEGG_diff()

# d.get_KEGG_diff

# d.REACTOME_diff()

# d.get_REACTOME_diff

# d.spec_diff()

# d.get_specificity_diff

# d.gi_diff()

# d.get_GI_diff

# d.network_diff()


# d.get_networks_diff

# d.inter_processes()

d.full_analysis()


results3 = d.get_results




with open('set1_set2.json', 'w') as json_file:
    json.dump(results3, json_file)





import json

with open('set1_set2.json', 'r') as json_file:
    results3 = (json.load(json_file))


input_data = results3


