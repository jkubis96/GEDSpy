
def DataPrepare_tests():
    
    from DataPrepare import Donwload
    
    errored = []

    dw = Donwload()
    
    #
    
    ref = None
    
    ref = dw.download_ref()
    
    if not isinstance(ref, dict):
        errored.append('download_ref')
    
    del ref
    
    #
    
    rns_seq = None
    
    rns_seq = dw.download_rns_seq()
    
    if not isinstance(rns_seq, dict):
        errored.append('download_rns_seq')
    
    del rns_seq
    
    #
    
    IntAct = None

    IntAct = dw.download_IntAct_data()
    
    if not isinstance(IntAct, dict):
        errored.append('download_IntAct_data')
    
    del IntAct
    
    
    #
    
    Viral = None

    Viral = dw.download_viral_deiseases()
    
    if not isinstance(Viral, dict):
        errored.append('download_viral_deiseases')
    
    del Viral
    
    
    
    #
    
    HPA = None

    HPA = dw.download_HPA()
    
    if not isinstance(HPA, dict):
        errored.append('download_HPA')
    
    del HPA
    
    
    
    #
    
    string = None

    string = dw.download_string()
    
    if not isinstance(string, dict):
        errored.append('download_string')
    
    del string
    
    
    
    
    #
    
    kegg = None

    kegg = dw.download_kegg()
    
    if not isinstance(kegg, dict):
        errored.append('download_kegg')
    
    del kegg
    
    #
    
    dis = None

    dis = dw.download_diseases()
    
    if not isinstance(dis, dict):
        errored.append('download_diseases')
    
    del dis

    #
    
    reactome = None

    reactome = dw.download_reactome()
    
    if not isinstance(reactome, dict):
        errored.append('download_reactome')
    
    del reactome
    
    
    
    #
    
    go = None

    go = dw.download_go_term()
    
    if not isinstance(go, dict):
        errored.append('download_go_term')
    
    del go
    

    
    #
    
    cp = None

    cp = dw.download_cell_phone()
    
    if not isinstance(cp, dict):
        errored.append('download_cell_phone')
    
    del cp
    
    
    
    #
    
    ct = None

    ct = dw.download_cell_talk()
    
    if not isinstance(ct, dict):
        errored.append('download_cell_talk')
    
    del ct
    
    del dw
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored
    

DataPrepare_tests()




###############################################################################


def EnrichmentGetRaw_tests():


    from Enrichment import GetDataRaw
    
    errored = []
    
    
    gdr = GetDataRaw()
    
    
    recteome = None
    
    recteome = gdr.get_raw_REACTOME()
    
    if not isinstance(recteome, dict):
        errored.append('get_raw_REACTOME')
    
    del recteome
    
    
    ref_gen = None

    ref_gen = gdr.get_raw_REF_GEN()
    
    if not isinstance(ref_gen, dict):
        errored.append('get_raw_REF_GEN')
    
    del ref_gen
    
    
    ref_gen_seq = None

    ref_gen_seq = gdr.get_raw_REF_GEN_RNA_SEQ()
    
    if not isinstance(ref_gen_seq, dict):
        errored.append('get_raw_REF_GEN_RNA_SEQ')
    
    del ref_gen_seq
    
    
    HPA = None

    HPA = gdr.get_raw_HPA()
    
    if not isinstance(HPA, dict):
        errored.append('get_raw_HPA')
    
    del HPA
    
    
    diseases = None

    diseases = gdr.get_raw_DISEASES()
    
    if not isinstance(diseases, dict):
        errored.append('get_raw_DISEASES')
    
    del diseases
    
    
    vimic = None

    vimic = gdr.get_raw_ViMIC()
    
    if not isinstance(vimic, dict):
        errored.append('get_raw_ViMIC')
    
    del vimic
     
    
    kegg = None

    kegg = gdr.get_raw_KEGG()
    
    if not isinstance(kegg, dict):
        errored.append('get_raw_KEGG')
    
    del kegg
    
    
    go = None

    go = gdr.get_raw_GO()
    
    if not isinstance(go, dict):
        errored.append('get_raw_GO')
    
    del go
    
    
    intact = None

    intact = gdr.get_raw_IntAct()
    
    if not isinstance(intact, dict):
        errored.append('get_raw_IntAct')
    
    del intact
    
    
    string = None

    string = gdr.get_raw_STRING()
    
    if not isinstance(string, dict):
        errored.append('get_raw_STRING')
    
    del string
    
    
    cell_talk = None

    cell_talk = gdr.get_raw_CellTalk()
    
    if not isinstance(cell_talk, dict):
        errored.append('get_raw_CellTalk')
    
    del cell_talk
    
    
    cell_phone = None

    cell_phone = gdr.get_raw_CellPhone()
    
    if not isinstance(cell_phone, dict):
        errored.append('get_raw_CellPhone')
    
    del cell_phone
    
    
    del gdr
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored


EnrichmentGetRaw_tests()


###############################################################################







def EnrichmentGetData_tests():


    from Enrichment import GetData
    
    errored = []
    
    
    gd = GetData()
    
    
    reactome = None
    
    reactome = gd.get_REACTOME()
    
    if not isinstance(reactome, dict):
        errored.append('get_REACTOME')
         
    
    ref_gen = None
    
    ref_gen = gd.get_REF_GEN()
    
    if not isinstance(ref_gen, dict):
        errored.append('get_REF_GEN')
        
    
    ref_gen_seq = None
    
    ref_gen_seq = gd.get_REF_GEN_RNA_SEQ()
    
    if not isinstance(ref_gen_seq, dict):
        errored.append('get_REF_GEN_RNA_SEQ')
         
      
    HPA = None
    
    HPA = gd.get_HPA()
    
    if not isinstance(HPA, dict):
        errored.append('get_HPA')
         
    
    disease = None
    
    disease = gd.get_DISEASES()
    
    if not isinstance(disease, dict):
        errored.append('get_DISEASES')
    
    
    vimic = None
    
    vimic = gd.get_ViMIC()
    
    if not isinstance(vimic, dict):
        errored.append('get_ViMIC')
         
    
    kegg = None
    
    kegg = gd.get_KEGG()
    
    if not isinstance(kegg, dict):
        errored.append('get_KEGG')
         
    
    go = None
    
    go = gd.get_GO()
    
    if not isinstance(go, dict):
        errored.append('get_GO')
         
    
    intact = None
    
    intact = gd.get_IntAct()
    
    if not isinstance(intact, dict):
        errored.append('get_IntAct')
         
    
    string = None
    
    string = gd.get_STRING()
    
    if not isinstance(string, dict):
        errored.append('get_STRING')
         
    
    CellTalk = None
    
    CellTalk = gd.get_CellTalk()
    
    if not isinstance(CellTalk, dict):
        errored.append('get_CellTalk')
    
    
    CellPhone = None
    
    CellPhone = gd.get_CellPhone()
    
    if not isinstance(CellPhone, dict):
        errored.append('get_CellPhone')
    
    
    CellInteractions = None
    
    CellInteractions = gd.get_interactions()
    
    if not isinstance(CellInteractions, dict):
        errored.append('get_interactions')
    
      
    del gd
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored
    

EnrichmentGetData_tests()



def Enrichment_tests():
    
    errored = []

    
    
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
    
        
    from Enrichment import Enrichment
    
    enr = Enrichment()
    
    
    enr.select_features(list1)
    

    gene_info = enr.get_gene_info
    
    if not isinstance(gene_info, dict):
        errored.append('select_features')
 

    
    
    print('\nSpecificity enrichment...')

    enr.enriche_specificiti()
    
    HPA = enr.get_HPA
    
    if not isinstance(HPA, dict):
        errored.append('enriche_specificiti')

    
    print('\nEnrichment with KEGG information...')

    enr.enriche_KEGG()
    
    KEGG = enr.get_KEGG
    
    if not isinstance(KEGG, dict):
        errored.append('enriche_KEGG')

    
    print('\nEnrichment with GO-TERM information...')

    enr.enriche_GOTERM()
    
    GOTERM = enr.get_GO_TERM
    
    if not isinstance(GOTERM, dict):
        errored.append('enriche_GOTERM')

    
    print('\nEnrichment with REACTOME information...')

    enr.enriche_REACTOME()
    
    REACTOME = enr.get_REACTOME
    
    if not isinstance(REACTOME, dict):
        errored.append('enriche_REACTOME')

    
    print('\nEnrichment with DISEASES information...')

    enr.enriche_DISEASES()
    
    DISEASES = enr.get_DISEASES
    
    if not isinstance(DISEASES, dict):
        errored.append('enriche_DISEASES')

    
    print('\nEnrichment with ViMIC information...')

    enr.enriche_ViMIC()
    
    ViMIC = enr.get_ViMIC
    
    if not isinstance(ViMIC, dict):
        errored.append('enriche_ViMIC')

    
    print('\nEnrichment with IntAct information...') 

    enr.enriche_IntAct()
    
    IntAct = enr.get_IntAct
    
    if not isinstance(IntAct, dict):
        errored.append('enriche_IntAct')

    
    print('\nEnrichment with STRING information...')

    enr.enriche_STRING()
    
    STRING = enr.get_STRING
    
    if not isinstance(STRING, dict):
        errored.append('enriche_STRING')

    
    print('\nEnrichment with CellConnections information...')

    enr.enriche_CellCon()
    
    CellConnections = enr.get_CellCon
    
    if not isinstance(CellConnections, dict):
        errored.append('enriche_CellCon')

    
    print('\nEnrichment with tissue specific RNA-SEQ information...')

    enr.enriche_RNA_SEQ()
    
    RNASEQ = enr.get_RNA_SEQ   
    
    if not isinstance(RNASEQ, dict):
        errored.append('enriche_RNA_SEQ')



    del enr
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored




Enrichment_tests()




def AnVis_tests():
    
    errored = []

    
    
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
    
        
    from Enrichment import Enrichment
    
    enr = Enrichment()
    
    
    enr.select_features(list1)
    
    
    # full
    
    enr.full_enrichment()

    
    results = None
    
    results = enr.get_results
    
    if not isinstance(results, dict):
        errored.append('full_enrichment')

    
    

    from Enrichment import Analysis
    
    
    ans = Analysis(results)
    
    ans.networks_metadata
    
    ans.interactions_metadata
    
    ans.set_p_value(value = 0.05)
    ans.set_test(test = 'FISH')       
    ans.set_correction(correction = None)    
     
    ans.networks_metadata
    
    
    #########################################################################
    
    
    print('\nGO-TERM overrepresentation analysis...')
    ans.GO_overrepresentation()
    
    go = None
    
    go = ans.get_GO_statistics
    
    if not isinstance(go, dict):
        errored.append('GO_overrepresentation')
        
    del go
    
    print('\nKEGG overrepresentation analysis...')
    ans.KEGG_overrepresentation()
    
    kegg = None
    
    kegg = ans.get_KEGG_statistics
    
    if not isinstance(kegg, dict):
        errored.append('KEGG_overrepresentation')
        
    del kegg
    
    
    
    
    print('\nREACTOME overrepresentation analysis...')
    ans.REACTOME_overrepresentation()
    
    
    reactome = None
    
    reactome = ans.get_REACTOME_statistics
    
    if not isinstance(reactome, dict):
        errored.append('REACTOME_overrepresentation')
        
    del reactome
    
    
    print('\nViMIC overrepresentation analysis...')
    ans.ViMIC_overrepresentation()
    
    vimic = None
    
    vimic = ans.get_ViMIC_statistics
    
    if not isinstance(vimic, dict):
        errored.append('ViMIC_overrepresentation')
        
    del vimic
    

    
    print('\nDISEASES overrepresentation analysis...')
    ans.DISEASES_overrepresentation()
    
    diseases = None
    
    diseases = ans.get_DISEASE_statistics
    
    if not isinstance(diseases, dict):
        errored.append('DISEASES_overrepresentation')
        
    del diseases
    

    
    print('\nSpecificity overrepresentation analysis...')
    ans.features_specificity()
    
    
    spec = None
    
    spec = ans.get_specificity_statistics
    
    if not isinstance(spec, dict):
        errored.append('features_specificity')
        
    del spec
    

    
    print('\nInteraction analysis...')
    ans.gene_interaction()
    
    inter = None
    
    inter = ans.get_features_interactions_statistics
    
    if not isinstance(inter, dict):
        errored.append('gene_interaction')
        
    del inter
    

    
    print('\nNetwork creating...')
    ans.REACTOME_network()
    
    reactome_net = None
    
    reactome_net = ans.get_REACTOME_network
    
    if not isinstance(reactome_net, dict):
        errored.append('REACTOME_network')
        
    del reactome_net
    

    ans.KEGG_network()
    
    kegg_net = None
    
    kegg_net = ans.get_KEGG_network
    
    if not isinstance(kegg_net, dict):
        errored.append('KEGG_network')
        
    del kegg_net


    ans.GO_network()
    
    go_net = None
    
    go_net = ans.get_GO_network
    
    if not isinstance(go_net, dict):
        errored.append('GO_network')
        
    del go_net


    
    # full analysis
    
    ans.full_analysis()
    
    results2 = None
    
    results2 = ans.get_full_results
    
    if not isinstance(results2, dict):
        errored.append('full_analysis')
        
    del results

    from matplotlib import figure


    from Enrichment import Visualization


    vis = Visualization(results2)
    
    plot = None
    
    plot  =  vis.gene_type_plot(cmap = 'summer', image_width = 6, image_high = 6, font_size = 15)

    if not isinstance(plot, figure.Figure):
        errored.append('gene_type_plot')
        
    del plot



    plot = None
    
    plot  =  vis.GO_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 20, 
                        side = 'left', 
                        color = 'blue', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')


    if not isinstance(plot, figure.Figure):
        errored.append('GO_plot')
        
    del plot


    plot = None
    
    plot  =  vis.SPECIFICITY_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 25, 
                        side = 'right', 
                        color = 'bisque', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')

    if not isinstance(plot, figure.Figure):
        errored.append('SPECIFICITY_plot')
        
    del plot


    plot = None
    
    plot  =  vis.KEGG_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 25, 
                        side = 'right', 
                        color = 'orange', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')
    
    
    if not isinstance(plot, figure.Figure):
        errored.append('KEGG_plot')
        
    del plot



    plot = None

    plot  =  vis.REACTOME_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 25, 
                        side = 'right', 
                        color = 'silver', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')


    if not isinstance(plot, figure.Figure):
        errored.append('REACTOME_plot')
        
    del plot
 

    plot = None

    plot  =  vis.DISEASES_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 25, 
                        side = 'right', 
                        color = 'thistle', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')


    if not isinstance(plot, figure.Figure):
        errored.append('DISEASES_plot')
        
    del plot
 


    plot = None
    
    plot  =  vis.ViMIC_plot(p_val = 0.05, 
                        test = 'FISH', 
                        adj = None, 
                        n = 25, 
                        side = 'right', 
                        color = 'aquamarine', 
                        width = 10, 
                        bar_width = 0.5, 
                        stat = 'p_val')


    if not isinstance(plot, figure.Figure):
        errored.append('ViMIC_plot')
        
    del plot


    plot = None

    plot  =  vis.blod_markers_plot()
    
    
    if not isinstance(plot, figure.Figure):
        errored.append('blod_markers_plot')
        
    del plot

    import networkx as nx

    go = None

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


    if not isinstance(go, nx.Graph):
        errored.append('GOPa_network_create-GO-TERM')
        
    del go


    
    kegg = None
    
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


    if not isinstance(kegg, nx.Graph):
        errored.append('GOPa_network_create-KEGG')
        
    del kegg


    react = None

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


    if not isinstance(react, nx.Graph):
        errored.append('GOPa_network_create-Reactome')
        
    del react



    inter = None
    
    inter = vis.GI_network_create()
    
    if not isinstance(inter, nx.Graph):
        errored.append('GI_network_create')
        
    del inter


    ml = None
    
    ml = vis.AUTO_ML_network( 
                        genes_inc = 10, 
                        gene_int = True, 
                        genes_only = True, 
                        min_con = 1, 
                        children_con = True, 
                        include_childrend = True,
                        selected_parents = ['Cell motility', 'Development and regeneration'],
                        selected_genes = ['JAK1', 'KRAS'])

    
    if not isinstance(ml, nx.Graph):
        errored.append('AUTO_ML_network')
        
    del ml



    seq = None
    
    seq = vis.gene_scatter( 
                     colors = 'viridis', 
                     species = 'human', 
                     hclust = 'complete', 
                     img_width = None, 
                     img_high = None, 
                     label_size = None, 
                     x_lab = 'Genes', 
                     legend_lab = 'log(TPM + 1)')

    
    if not isinstance(seq, dict):
        errored.append('gene_scatter')
        
    del seq


      
    del enr
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored



AnVis_tests()



def AnDESVis_tests():


    # DSA
    errored = []

    
    
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
    
    from Enrichment import Enrichment
    from Enrichment import Analysis


    #SET1
    
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



    #SET2

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




    from Enrichment import DSA


    dsa_compare = DSA(results1, results2)
    
     
    print('\nGO-TERM differential analysis...')

    dsa_compare.GO_diff()
    
    tmp = dsa_compare.get_GO_diff

    print('\nKEGG differential analysis...')

    dsa_compare.KEGG_diff()
    
    tmp = dsa_compare.get_KEGG_diff
    
    print('\nREACTOME differential analysis...')

    dsa_compare.REACTOME_diff()
    
    tmp = dsa_compare.get_REACTOME_diff

    print('\nSpecificity differential analysis...')

    dsa_compare.spec_diff()
    
    tmp = dsa_compare.get_specificity_diff
    
    print('\nGene Interactions (GI) differential analysis...')

    dsa_compare.gi_diff()
    
    tmp = dsa_compare.get_GI_diff
    
    print('\nNetworks differential analysis...')

    dsa_compare.network_diff()
    
    tmp = dsa_compare.get_networks_diff
    
    print('\nInter Terms (IT) searching...')

    dsa_compare.inter_processes()
    
    tmp = dsa_compare.get_inter_terms
    
    
    print('\nInter CellConnections (ICC) searching...')

    dsa_compare.connections_diff()
    
    tmp = dsa_compare.get_set_to_set_con
    
    
    
    

    
    
    
    
    dsa_compare.full_analysis()


    results3 = dsa_compare.get_results
    
    
    
    del enr
    
    if len(errored) == 0:
        
        print('\nTest passed!')
        
    else:
        print('\nTest not passed!')
        print('\nError:')

        for r in errored:
            print(f'\n -{r}')

        
    return errored
    
    

 
    
 
    
 
    

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






with open('set1_set2.json', 'w') as json_file:
    json.dump(results3, json_file)





import json

with open('set1_set2.json', 'r') as json_file:
    results3 = (json.load(json_file))


input_data = results3
