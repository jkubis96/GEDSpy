# Example analysis gene lists

gene_list = ['CACNA1I','CALD1','CAMK1G','CAMK2N1','CAMSAP1','CCL15','CCL16','CCNL2','CCT8P1','CD46','CDC14A','CDK18','CDK19','CES3','CHEK2',
			 'CHID1','COL6A3','CPVL','CYP3A43','CYP3A5','DBNL','DUSP10','DUSP9','ECHS1','EGFR','EGR2','ELL2','ERMP1','ESR1','F7','FAM171A1',
			 'FAM20C','FGFR2','FH','FLAD1','FUT3','GAA','GBA2','GGCX','GJB1','GLRX5','GNAI2','GNB2','GNB3','GPNMB','GRB10','GRHPR','HMGCS2',
			 'HSD17B4','HSP90AB1','IER3IP1','IGF2R','IL1R1','INF2','IRAK1','ITGA1','ITGA7','ITIH1','ITIH3','ITIH4','ITPR1','ITSN1','JAK1',
			 'KALRN','KCNQ2','KCNQ4','KDM3A','KMO','KRAS','KSR1','LAMA5','LAMB2','LCN2','MAP2K7','MAP4K2','MAP4K3',
			 'MAPK13','MARCO','MAST2','MAT1A','MATR3','MCM8','MFSD10','MGAT5','MTMR10','MUSK','MYO9B','NBAS']

gene_list2 = ['NCOA6','NCSTN','NDUFA4','NEK4', 'UDT2','NUP210','ORC3L','PAOX','PEMT','PEX14','PFKL','PHKA2','PIM1','PLXND1','PMM1','PON3','POR','PPARG','PPARGC1B',
              'PRKCE','PTK2B','PTP4A1','PTPN23','PTPRF','PTPRK','RARA','RNF10','RNF14','RNF165','ROCK2','RRBP1','RREB1','SCN1A','SDC1','SERPINA1','SERPINA10','SFXN5',
              'SHROOM1','SIL1','SIRPA','SLC12A7','SLC13A3','SLC16A2','SLC17A7','SLC22A23','SLC22A9','SLC25A11','SLC25A25','SLC38A3','SLC45A3','SLC4A5','SLC5A1','SLC7A2',
              'SLC8A3','SLC9A6','SLCO1A2','SLCO1B3','SMARCA2','NX4','SORBS1','SPEN','SPR','SRF','STAB1','STAT1','SUCLG2','SULT1B1','SULT1E1','TBC1D2B','TCHP','TGFBI',
              'TGOLN2','THPO', 'IMM13','TLK2','TMEM62','TNFSF14','TNK2','TNS1','TPI1','TRIB3','TRMT11','TTYH3']
    
    



# load classes
from Enrichment import Enrichment
from Enrichment import Analysis


#SET1

# create class Enrichment
enr = Enrichment()

# select featrues from genes/proteins list for Homo sapiens / Mus musculus / Rattus norvegicus for first gene_list
enr.select_features(gene_list)

# conduct full enrichment
enr.full_enrichment()

# get results for first set
results1 = enr.get_results

    
# create class Analysis
ans = Analysis(results1)

# set parameters or leave default - see manual on GitHub
ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    

# ! parameters for both analyses (set1 and set2) must be the same ! 
 
# conduct full analyses
ans.full_analysis()

# get full results for first set
results1 = ans.get_full_results



#SET2

# create class Enrichment
enr = Enrichment()

# select featrues from genes/proteins list for Homo sapiens / Mus musculus / Rattus norvegicus for second gene_list
enr.select_features(gene_list2)

# conduct full enrichment
enr.full_enrichment()

# get results for second set
results2 = enr.get_results


# create class Analysis
ans = Analysis(results2)


# set parameters or leave default - see manual on GitHub
ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    
 
# ! parameters for both analyses (set1 and set2) must be the same ! 
 
# conduct full analyses
ans.full_analysis()

# get full results for second set
results2 = ans.get_full_results




# load class DSA
from Enrichment import DSA

# create class with results from get_full_results for set1 and set2
dsa_compare = DSA(results1, results2)


# DSA analysis
# There are two way of enrich data:
# * each data separately:
     
# GO-TERM DSA
dsa_compare.GO_diff()

go = dsa_compare.get_GO_diff

# KEGG DSA
dsa_compare.KEGG_diff()

kegg = dsa_compare.get_KEGG_diff

# Reactome DSA
dsa_compare.REACTOME_diff()

reactome = dsa_compare.get_REACTOME_diff

# Specificity (HPA) DSA
dsa_compare.spec_diff()

hpa = dsa_compare.get_specificity_diff

# Gene interactions DSA
dsa_compare.gi_diff()

GI = dsa_compare.get_GI_diff

# Networks DSA
dsa_compare.network_diff()

networks = dsa_compare.get_networks_diff

# Inter Terms DSA
dsa_compare.inter_processes()

inter_terms = dsa_compare.get_inter_terms

# Connections DSA
dsa_compare.connections_diff()

connections = dsa_compare.get_set_to_set_con



# * all data in parallel:

dsa_compare.full_analysis()


# get_results can be used for return all DSA data like above ^^^
results3 = dsa_compare.get_results



# load class VisualizationDES
from Enrichment import VisualizationDES

# load library JVG - display and adjustment Networks and Bar plots
from JVG import JVG

# create class DAA with results from get_results for set1 and set2 DSA results
vis_des = VisualizationDES(results3)



# Gene type
plot  =  vis_des.diff_gene_type_plot( 
                        set1_name = 'Set 1', 
                        set2_name = 'Set 2', 
                        image_width = 12, 
                        image_high = 6, 
                        font_size = 15)

graph = JVG.MplEditor(plot)
graph.edit()


# GO-TERM
plot  =  vis_des.diff_GO_plot(   
            p_val = 0.05, 
            test = 'FISH', 
            adj = 'BH', 
            n = 25, 
            min_terms = 5,
            selected_parent = [],
            width = 10, 
            bar_width = 0.5, 
            stat = 'p_val',
            sep_factor = 50)

graph = JVG.MplEditor(plot)
graph.edit()


# KEGG
plot  =  vis_des.diff_KEGG_plot(   
                p_val = 0.05, 
                test = 'FISH', 
                adj = 'BH', 
                n = 25, 
                min_terms = 2,
                selected_parent = [],
                width = 10, 
                bar_width = 0.5, 
                stat = 'p_val',
                sep_factor = 50)

graph = JVG.MplEditor(plot)
graph.edit()


# Reactome
plot  =  vis_des.diff_REACTOME_plot(   
                    p_val = 0.05, 
                    test = 'FISH', 
                    adj = 'BH', 
                    n = 25, 
                    min_terms = 5,
                    selected_parent = [],
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'n',
                    sep_factor = 50)

graph = JVG.MplEditor(plot)
graph.edit()


# Specificity (HPA)
plot  =  vis_des.diff_SPECIFICITY_plot(   
                    p_val = 0.05, 
                    test = 'FISH', 
                    adj = 'BH', 
                    n = 6, 
                    min_terms = 1,
                    selected_set = [],
                    width = 10, 
                    bar_width = 0.5, 
                    stat = 'p_val',
                    sep_factor = 15)

graph = JVG.MplEditor(plot)
graph.edit()


# GO-TERM network
go = vis_des.diff_GOPa_network_create(data_set = 'GO-TERM', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = False,
                        include_childrend = True,
                        selected_parents = [],
                        selected_genes = [])

nt = JVG.NxEditor(go)
nt.edit()


# KEGG network
kegg = vis_des.diff_GOPa_network_create(data_set = 'KEGG', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = False,
                        include_childrend = True,
                        selected_parents = [],
                        selected_genes = [])

nt = JVG.NxEditor(kegg)
nt.edit()


# Reactome network
reactome = vis_des.diff_GOPa_network_create(data_set = 'REACTOME', 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 2, 
                        children_con = False,
                        include_childrend = True,
                        selected_parents = [],
                        selected_genes = [])

nt = JVG.NxEditor(reactome)
nt.edit()


# GI network
gi = vis_des.diff_GI_network_create(min_con = 2)

nt = JVG.NxEditor(gi)
nt.edit()


# AMUTO-ML network
ml = vis_des.diff_AUTO_ML_network( 
                    genes_inc = 10, 
                    gene_int = True, 
                    genes_only = True, 
                    min_con = 2, 
                    children_con = False, 
                    include_childrend = False,
                    selected_parents = [],
                    selected_genes = [])
    
nt = JVG.NxEditor(ml)
nt.edit()


# RNAseq
seq = vis_des.diff_gene_scatter( 
                 set_num = 2,
                 colors = 'viridis', 
                 species = 'human', 
                 hclust = 'complete', 
                 img_width = None, 
                 img_high = None, 
                 label_size = None, 
                 x_lab = 'Genes', 
                 legend_lab = 'log(TPM + 1)',
                 selected_list = [])

   
graph1 = JVG.MplEditor(seq['tissue_expression_HPA'])
graph1.edit()

graph2 = JVG.MplEditor(seq['tissue_expression_RNA_total_tissue'])
graph2.edit()

graph3 = JVG.MplEditor(seq['tissue_expression_fetal_development_circular'])
graph3.edit()

  