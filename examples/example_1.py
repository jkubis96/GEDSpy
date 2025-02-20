# Example analysis gene list   

gene_list = ['CACNA1I','CALD1','CAMK1G','CAMK2N1','CAMSAP1','CCL15','CCL16','CCNL2','CCT8P1','CD46','CDC14A','CDK18','CDK19','CES3','CHEK2',
			 'CHID1','COL6A3','CPVL','CYP3A43','CYP3A5','DBNL','DUSP10','DUSP9','ECHS1','EGFR','EGR2','ELL2','ERMP1','ESR1','F7','FAM171A1',
			 'FAM20C','FGFR2','FH','FLAD1','FUT3','GAA','GBA2','GGCX','GJB1','GLRX5','GNAI2','GNB2','GNB3','GPNMB','GRB10','GRHPR','HMGCS2',
			 'HSD17B4','HSP90AB1','IER3IP1','IGF2R','IL1R1','INF2','IRAK1','ITGA1','ITGA7','ITIH1','ITIH3','ITIH4','ITPR1','ITSN1','JAK1',
			 'KALRN','KCNQ2','KCNQ4','KDM3A','KMO','KRAS','KSR1','LAMA5','LAMB2','LCN2','MAP2K7','MAP4K2','MAP4K3',
			 'MAPK13','MARCO','MAST2','MAT1A','MATR3','MCM8','MFSD10','MGAT5','MTMR10','MUSK','MYO9B','NBAS']

    
    
from Enrichment import Enrichment

# load class Enrichment
enr = Enrichment()


# select featrues from genes/proteins list for Homo sapiens / Mus musculus / Rattus norvegicus
enr.select_features(gene_list)


# get selected genes/proteins info
gene_info = enr.get_gene_info



# There are two way of enrich data:
# * each data separately:
    
# Human protein Atlas data
enr.enriche_specificiti()

HPA = enr.get_HPA

# KEGG data
enr.enriche_KEGG()

KEGG = enr.get_KEGG

# GO-TERM data
enr.enriche_GOTERM()

GOTERM = enr.get_GO_TERM

# Reactome data
enr.enriche_REACTOME()

REACTOME = enr.get_REACTOME

# Human Diseases data
enr.enriche_DISEASES()

DISEASES = enr.get_DISEASES

# ViMic data
enr.enriche_ViMIC()

ViMIC = enr.get_ViMIC

# IntAct data
enr.enriche_IntAct()

IntAct = enr.get_IntAct

# STRING data
enr.enriche_STRING()

STRING = enr.get_STRING

# CellCon data
enr.enriche_CellCon()

CellConnections = enr.get_CellCon

# RNAseq data
enr.enriche_RNA_SEQ()

RNASEQ = enr.get_RNA_SEQ   



# * all data in parallel:


enr.full_enrichment()


# get_results can be used for return all data enriched separately like above ^^^
results = enr.get_results
 


# load class Analysis
from Enrichment import Analysis

# create class with results from get_results
ans = Analysis(results)

# adjustment of analysis parameters
# for more details go to documentation on GitHub

# network parameters
ans.networks_metadata

# interactions parameters
ans.interactions_metadata

# set analysis parameter (default parameters check in documentation on GitHub)
ans.set_p_value(value = 0.05)
ans.set_test(test = 'FISH')       
ans.set_correction(correction = None)    
 
ans.networks_metadata


#########################################################################
# Overrepresentation analysis
# There are two way of enrich data:
# * each data separately:

# GO-TERM
ans.GO_overrepresentation()

go = ans.get_GO_statistics

# KEGG
ans.KEGG_overrepresentation()

kegg = ans.get_KEGG_statistics

# Reactome
ans.REACTOME_overrepresentation()

reactome = ans.get_REACTOME_statistics

# ViMic
ans.ViMIC_overrepresentation()

vimic = ans.get_ViMIC_statistics

# Human diseases
ans.DISEASES_overrepresentation()

diseases = ans.get_DISEASE_statistics

# Specificity (HPA)
ans.features_specificity()

spec = ans.get_specificity_statistics

# Interactions (STRING/IntAct)
ans.gene_interaction()

inter = ans.get_features_interactions_statistics

# Reactome paths network
ans.REACTOME_network()

reactome_net = ans.get_REACTOME_network

# KEGG paths network
ans.KEGG_network()

kegg_net = ans.get_KEGG_network

# GO-TERM terms network
ans.GO_network()

go_net = ans.get_GO_network



# * all data in parallel:


ans.full_analysis()


# get_results can be used for return all data enriched separately like above ^^^
results2 = ans.get_full_results






# load class Visualization
from Enrichment import Visualization

# load library JVG - display and adjustment Networks and Bar plots
from JVG import JVG

# create class with results from get_full_results
vis = Visualization(results2)


# gene type
plot  =  vis.gene_type_plot(cmap = 'summer', image_width = 6, image_high = 6, font_size = 15)

graph = JVG.MplEditor(plot)
graph.edit()


# GO-TERM
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


# Specificity (HPA)
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


# KEGG
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


# Reactome
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


# Human diseases
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


# ViMic
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


# blood markers
plot  =  vis.blod_markers_plot()

graph = JVG.MplEditor(plot)
graph.edit()


# GO-TERM network
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


# KEGG network
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


# Reactome network
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


# Gene interactions network
inter = vis.GI_network_create()

nt = JVG.NxEditor(inter)
nt.edit()


# AUTO-ML
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


# RNAseq
seq = vis.gene_scatter( 
                 colors = 'viridis', 
                 species = 'human', 
                 hclust = 'complete', 
                 img_width = None, 
                 img_high = None, 
                 label_size = None, 
                 x_lab = 'Genes', 
                 legend_lab = 'log(TPM + 1)')


graph1 = JVG.MplEditor(seq['tissue_expression_HPA'])
graph1.edit()

graph2 = JVG.MplEditor(seq['tissue_expression_RNA_total_tissue'])
graph2.edit()

graph3 = JVG.MplEditor(seq['tissue_expression_fetal_development_circular'])
graph3.edit()


