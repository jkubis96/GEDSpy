import networkx as nx
import pytest
from matplotlib import figure

from gedspy import (
    DSA,
    Analysis,
    Enrichment,
    GetData,
    GetDataRaw,
    Visualization,
    VisualizationDES,
)


def test_GetDataRaw():

    gdr = GetDataRaw()

    recteome = None

    recteome = gdr.get_raw_REACTOME()

    assert isinstance(recteome, dict)


def test_GetData():

    gd = GetData()

    reactome = None

    reactome = gd.get_REACTOME()

    assert isinstance(reactome, dict)


def test_Enrichment():

    gene_list = [
        "CACNA1I",
        "CALD1",
        "CAMK1G",
        "CAMK2N1",
        "CAMSAP1",
        "CCL15",
        "CCL16",
        "CCNL2",
        "CCT8P1",
        "CD46",
        "CDC14A",
        "CDK18",
        "CDK19",
        "CES3",
        "CHEK2",
        "CHID1",
        "COL6A3",
        "CPVL",
        "CYP3A43",
        "CYP3A5",
        "DBNL",
        "DUSP10",
        "DUSP9",
        "ECHS1",
        "EGFR",
        "EGR2",
        "ELL2",
        "ERMP1",
        "ESR1",
        "F7",
        "FAM171A1",
        "FAM20C",
        "FGFR2",
        "FH",
        "FLAD1",
        "FUT3",
        "GAA",
        "GBA2",
        "GGCX",
        "GJB1",
        "GLRX5",
        "GNAI2",
        "GNB2",
        "GNB3",
        "GPNMB",
        "GRB10",
        "GRHPR",
        "HMGCS2",
        "HSD17B4",
        "HSP90AB1",
        "IER3IP1",
        "IGF2R",
        "IL1R1",
        "INF2",
        "IRAK1",
        "ITGA1",
        "ITGA7",
        "ITIH1",
        "ITIH3",
        "ITIH4",
        "ITPR1",
        "ITSN1",
        "JAK1",
        "KALRN",
        "KCNQ2",
        "KCNQ4",
        "KDM3A",
        "KMO",
        "KRAS",
        "KSR1",
        "LAMA5",
        "LAMB2",
        "LCN2",
        "MAP2K7",
        "MAP4K2",
        "MAP4K3",
        "MAPK13",
        "MARCO",
        "MAST2",
        "MAT1A",
        "MATR3",
        "MCM8",
        "MFSD10",
        "MGAT5",
        "MTMR10",
        "MUSK",
        "MYO9B",
        "NBAS",
        "NCOA6",
        "NCSTN",
        "NDUFA4",
        "NEK4",
        "NPR2",
        "NUDT2",
        "NUP210",
        "ORC3L",
        "PAOX",
        "PEMT",
        "PEX14",
        "PFKL",
        "PHKA2",
        "PIM1",
        "PLXND1",
        "PMM1",
        "PON3",
        "POR",
        "PPARG",
        "PPARGC1B",
        "PPP2R1A",
        "PRKCE",
        "PTK2B",
        "PTP4A1",
        "PTPN23",
        "PTPRF",
        "PTPRK",
        "RARA",
        "RNF10",
        "RNF14",
        "RNF165",
        "ROCK2",
        "RRBP1",
        "RREB1",
        "SCN1A",
        "SDC1",
        "SEPHS1",
        "SERPINA1",
        "SERPINA10",
        "SFXN5",
        "SHROOM1",
        "SIL1",
        "SIRPA",
        "SLC12A7",
        "SLC13A3",
        "SLC16A2",
        "SLC17A7",
        "SLC22A23",
        "SLC22A9",
        "SLC23A2",
        "SLC25A11",
        "SLC25A25",
        "SLC38A3",
        "SLC45A3",
        "SLC4A5",
        "SLC5A1",
        "SLC7A2",
        "SLC8A3",
        "SLC9A6",
        "SLCO1A2",
        "SLCO1B3",
        "SMARCA2",
        "SNRK",
        "SNX4",
        "SORBS1",
        "SPEN",
        "SPR",
        "SRF",
        "STAB1",
        "STAT1",
        "SUCLG2",
        "SULT1B1",
        "SULT1E1",
        "TBC1D2B",
        "TCHP",
        "TGFBI",
        "TGOLN2",
        "THPO",
        "TIE1",
        "TIMM13",
        "TLK2",
        "TMEM62",
        "TNFSF14",
        "TNK2",
        "TNS1",
        "TPI1",
        "TRIB3",
        "TRMT11",
        "TTYH3",
    ]

    # Split the list into two halves
    half = len(gene_list) // 2
    list1 = gene_list[:half]

    enr = Enrichment()

    enr.select_features(list1)

    gene_info = enr.get_gene_info()

    assert isinstance(gene_info, dict)

    enr.enriche_specificiti()

    HPA = enr.get_HPA()

    assert isinstance(HPA, dict)

    enr.enriche_KEGG()

    KEGG = enr.get_KEGG()

    assert isinstance(KEGG, dict)

    enr.enriche_GOTERM()

    GOTERM = enr.get_GO_TERM()

    assert isinstance(GOTERM, dict)

    enr.enriche_REACTOME()

    REACTOME = enr.get_REACTOME()

    assert isinstance(REACTOME, dict)

    enr.enriche_DISEASES()

    DISEASES = enr.get_DISEASES()

    assert isinstance(DISEASES, dict)

    enr.enriche_ViMIC()

    ViMIC = enr.get_ViMIC()

    assert isinstance(ViMIC, dict)

    enr.enriche_IntAct()

    IntAct = enr.get_IntAct()

    assert isinstance(IntAct, dict)

    enr.enriche_STRING()

    STRING = enr.get_STRING()

    assert isinstance(STRING, dict)

    enr.enriche_CellCon()

    CellConnections = enr.get_CellCon()

    assert isinstance(CellConnections, dict)

    enr.enriche_RNA_SEQ()

    RNASEQ = enr.get_RNA_SEQ()

    assert isinstance(RNASEQ, dict)

    del enr


def test_AnVis():

    gene_list = [
        "CACNA1I",
        "CALD1",
        "CAMK1G",
        "CAMK2N1",
        "CAMSAP1",
        "CCL15",
        "CCL16",
        "CCNL2",
        "CCT8P1",
        "CD46",
        "CDC14A",
        "CDK18",
        "CDK19",
        "CES3",
        "CHEK2",
        "CHID1",
        "COL6A3",
        "CPVL",
        "CYP3A43",
        "CYP3A5",
        "DBNL",
        "DUSP10",
        "DUSP9",
        "ECHS1",
        "EGFR",
        "EGR2",
        "ELL2",
        "ERMP1",
        "ESR1",
        "F7",
        "FAM171A1",
        "FAM20C",
        "FGFR2",
        "FH",
        "FLAD1",
        "FUT3",
        "GAA",
        "GBA2",
        "GGCX",
        "GJB1",
        "GLRX5",
        "GNAI2",
        "GNB2",
        "GNB3",
        "GPNMB",
        "GRB10",
        "GRHPR",
        "HMGCS2",
        "HSD17B4",
        "HSP90AB1",
        "IER3IP1",
        "IGF2R",
        "IL1R1",
        "INF2",
        "IRAK1",
        "ITGA1",
        "ITGA7",
        "ITIH1",
        "ITIH3",
        "ITIH4",
        "ITPR1",
        "ITSN1",
        "JAK1",
        "KALRN",
        "KCNQ2",
        "KCNQ4",
        "KDM3A",
        "KMO",
        "KRAS",
        "KSR1",
        "LAMA5",
        "LAMB2",
        "LCN2",
        "MAP2K7",
        "MAP4K2",
        "MAP4K3",
        "MAPK13",
        "MARCO",
        "MAST2",
        "MAT1A",
        "MATR3",
        "MCM8",
        "MFSD10",
        "MGAT5",
        "MTMR10",
        "MUSK",
        "MYO9B",
        "NBAS",
        "NCOA6",
        "NCSTN",
        "NDUFA4",
        "NEK4",
        "NPR2",
        "NUDT2",
        "NUP210",
        "ORC3L",
        "PAOX",
        "PEMT",
        "PEX14",
        "PFKL",
        "PHKA2",
        "PIM1",
        "PLXND1",
        "PMM1",
        "PON3",
        "POR",
        "PPARG",
        "PPARGC1B",
        "PPP2R1A",
        "PRKCE",
        "PTK2B",
        "PTP4A1",
        "PTPN23",
        "PTPRF",
        "PTPRK",
        "RARA",
        "RNF10",
        "RNF14",
        "RNF165",
        "ROCK2",
        "RRBP1",
        "RREB1",
        "SCN1A",
        "SDC1",
        "SEPHS1",
        "SERPINA1",
        "SERPINA10",
        "SFXN5",
        "SHROOM1",
        "SIL1",
        "SIRPA",
        "SLC12A7",
        "SLC13A3",
        "SLC16A2",
        "SLC17A7",
        "SLC22A23",
        "SLC22A9",
        "SLC23A2",
        "SLC25A11",
        "SLC25A25",
        "SLC38A3",
        "SLC45A3",
        "SLC4A5",
        "SLC5A1",
        "SLC7A2",
        "SLC8A3",
        "SLC9A6",
        "SLCO1A2",
        "SLCO1B3",
        "SMARCA2",
        "SNRK",
        "SNX4",
        "SORBS1",
        "SPEN",
        "SPR",
        "SRF",
        "STAB1",
        "STAT1",
        "SUCLG2",
        "SULT1B1",
        "SULT1E1",
        "TBC1D2B",
        "TCHP",
        "TGFBI",
        "TGOLN2",
        "THPO",
        "TIE1",
        "TIMM13",
        "TLK2",
        "TMEM62",
        "TNFSF14",
        "TNK2",
        "TNS1",
        "TPI1",
        "TRIB3",
        "TRMT11",
        "TTYH3",
    ]

    # Split the list into two halves
    half = len(gene_list) // 2
    list1 = gene_list[:half]

    enr = Enrichment()

    enr.select_features(list1)

    # full

    enr.full_enrichment()

    results = None

    results = enr.get_results()

    assert isinstance(results, dict)

    ans = Analysis(results)

    ans.networks_metadata()

    ans.interactions_metadata()

    ans.set_p_value(value=0.05)
    ans.set_test(test="FISH")
    ans.set_correction(correction=None)

    ans.networks_metadata()

    ans.GO_overrepresentation()

    go = None

    go = ans.get_GO_statistics()

    assert isinstance(go, dict)

    del go

    ans.KEGG_overrepresentation()

    kegg = None

    kegg = ans.get_KEGG_statistics()

    assert isinstance(kegg, dict)

    del kegg

    ans.REACTOME_overrepresentation()

    reactome = None

    reactome = ans.get_REACTOME_statistics()

    assert isinstance(reactome, dict)

    del reactome

    ans.ViMIC_overrepresentation()

    vimic = None

    vimic = ans.get_ViMIC_statistics()

    assert isinstance(vimic, dict)

    del vimic

    ans.DISEASES_overrepresentation()

    diseases = None

    diseases = ans.get_DISEASE_statistics()

    assert isinstance(diseases, dict)

    del diseases

    ans.features_specificity()

    spec = None

    spec = ans.get_specificity_statistics()

    assert isinstance(spec, dict)

    del spec

    ans.gene_interaction()

    inter = None

    inter = ans.get_features_interactions_statistics()

    assert isinstance(inter, dict)

    del inter

    ans.REACTOME_network()

    reactome_net = None

    reactome_net = ans.get_REACTOME_network()

    assert isinstance(reactome_net, dict)

    del reactome_net

    ans.KEGG_network()

    kegg_net = None

    kegg_net = ans.get_KEGG_network()

    assert isinstance(kegg_net, dict)

    del kegg_net

    ans.GO_network()

    go_net = None

    go_net = ans.get_GO_network()

    assert isinstance(go_net, dict)

    del go_net

    # full analysis

    ans.full_analysis()

    results2 = None

    results2 = ans.get_full_results()

    assert isinstance(results2, dict)

    del results

    vis = Visualization(results2)

    plot = None

    plot = vis.gene_type_plot(cmap="summer", image_width=6, image_high=6, font_size=15)

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.GO_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=20,
        side="left",
        color="blue",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.SPECIFICITY_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=25,
        side="right",
        color="bisque",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.KEGG_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=25,
        side="right",
        color="orange",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.REACTOME_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=25,
        side="right",
        color="silver",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.DISEASES_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=25,
        side="right",
        color="thistle",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.ViMIC_plot(
        p_val=0.05,
        test="FISH",
        adj=None,
        n=25,
        side="right",
        color="aquamarine",
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis.blod_markers_plot()

    assert isinstance(plot, figure.Figure)

    del plot

    go = None

    go = vis.GOPa_network_create(
        data_set="GO-TERM",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=True,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )

    assert isinstance(go, nx.Graph)

    del go

    kegg = None

    kegg = vis.GOPa_network_create(
        data_set="KEGG",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=True,
        include_childrend=True,
        selected_parents=["Nervous system", "Cell motility"],
        selected_genes=[],
    )

    assert isinstance(kegg, nx.Graph)

    del kegg

    react = None

    react = vis.GOPa_network_create(
        data_set="REACTOME",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=False,
        selected_parents=[],
        selected_genes=["F7", "LAMA5"],
    )

    assert isinstance(react, nx.Graph)

    del react

    inter = None

    inter = vis.GI_network_create()

    assert isinstance(inter, nx.Graph)

    del inter

    ml = None

    ml = vis.AUTO_ML_network(
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=1,
        children_con=True,
        include_childrend=True,
        selected_parents=["Cell motility", "Development and regeneration"],
        selected_genes=["JAK1", "KRAS"],
    )

    assert isinstance(ml, nx.Graph)

    del ml

    seq = None

    seq = vis.gene_scatter(
        colors="viridis",
        species="human",
        hclust="complete",
        img_width=None,
        img_high=None,
        label_size=None,
        x_lab="Genes",
        legend_lab="log(TPM + 1)",
    )

    assert isinstance(seq, dict)

    del seq

    del enr


def AnDESVis_tests():

    # DSA
    errored = []

    gene_list = [
        "CACNA1I",
        "CALD1",
        "CAMK1G",
        "CAMK2N1",
        "CAMSAP1",
        "CCL15",
        "CCL16",
        "CCNL2",
        "CCT8P1",
        "CD46",
        "CDC14A",
        "CDK18",
        "CDK19",
        "CES3",
        "CHEK2",
        "CHID1",
        "COL6A3",
        "CPVL",
        "CYP3A43",
        "CYP3A5",
        "DBNL",
        "DUSP10",
        "DUSP9",
        "ECHS1",
        "EGFR",
        "EGR2",
        "ELL2",
        "ERMP1",
        "ESR1",
        "F7",
        "FAM171A1",
        "FAM20C",
        "FGFR2",
        "FH",
        "FLAD1",
        "FUT3",
        "GAA",
        "GBA2",
        "GGCX",
        "GJB1",
        "GLRX5",
        "GNAI2",
        "GNB2",
        "GNB3",
        "GPNMB",
        "GRB10",
        "GRHPR",
        "HMGCS2",
        "HSD17B4",
        "HSP90AB1",
        "IER3IP1",
        "IGF2R",
        "IL1R1",
        "INF2",
        "IRAK1",
        "ITGA1",
        "ITGA7",
        "ITIH1",
        "ITIH3",
        "ITIH4",
        "ITPR1",
        "ITSN1",
        "JAK1",
        "KALRN",
        "KCNQ2",
        "KCNQ4",
        "KDM3A",
        "KMO",
        "KRAS",
        "KSR1",
        "LAMA5",
        "LAMB2",
        "LCN2",
        "MAP2K7",
        "MAP4K2",
        "MAP4K3",
        "MAPK13",
        "MARCO",
        "MAST2",
        "MAT1A",
        "MATR3",
        "MCM8",
        "MFSD10",
        "MGAT5",
        "MTMR10",
        "MUSK",
        "MYO9B",
        "NBAS",
        "NCOA6",
        "NCSTN",
        "NDUFA4",
        "NEK4",
        "NPR2",
        "NUDT2",
        "NUP210",
        "ORC3L",
        "PAOX",
        "PEMT",
        "PEX14",
        "PFKL",
        "PHKA2",
        "PIM1",
        "PLXND1",
        "PMM1",
        "PON3",
        "POR",
        "PPARG",
        "PPARGC1B",
        "PPP2R1A",
        "PRKCE",
        "PTK2B",
        "PTP4A1",
        "PTPN23",
        "PTPRF",
        "PTPRK",
        "RARA",
        "RNF10",
        "RNF14",
        "RNF165",
        "ROCK2",
        "RRBP1",
        "RREB1",
        "SCN1A",
        "SDC1",
        "SEPHS1",
        "SERPINA1",
        "SERPINA10",
        "SFXN5",
        "SHROOM1",
        "SIL1",
        "SIRPA",
        "SLC12A7",
        "SLC13A3",
        "SLC16A2",
        "SLC17A7",
        "SLC22A23",
        "SLC22A9",
        "SLC23A2",
        "SLC25A11",
        "SLC25A25",
        "SLC38A3",
        "SLC45A3",
        "SLC4A5",
        "SLC5A1",
        "SLC7A2",
        "SLC8A3",
        "SLC9A6",
        "SLCO1A2",
        "SLCO1B3",
        "SMARCA2",
        "SNRK",
        "SNX4",
        "SORBS1",
        "SPEN",
        "SPR",
        "SRF",
        "STAB1",
        "STAT1",
        "SUCLG2",
        "SULT1B1",
        "SULT1E1",
        "TBC1D2B",
        "TCHP",
        "TGFBI",
        "TGOLN2",
        "THPO",
        "TIE1",
        "TIMM13",
        "TLK2",
        "TMEM62",
        "TNFSF14",
        "TNK2",
        "TNS1",
        "TPI1",
        "TRIB3",
        "TRMT11",
        "TTYH3",
    ]

    # Split the list into two halves
    half = len(gene_list) // 2
    list1 = gene_list[:half]
    list2 = gene_list[half:]

    # ENRICHMENT

    # SET1

    enr = Enrichment()

    enr.select_features(list1)

    enr.full_enrichment()

    results1 = None

    results1 = enr.get_results()

    assert isinstance(results1, dict)

    ans = Analysis(results1)

    ans.networks_metadata()

    ans.interactions_metadata()

    ans.set_p_value(value=0.05)
    ans.set_test(test="FISH")
    ans.set_correction(correction=None)

    ans.full_analysis()

    results1 = None

    results1 = ans.get_full_results()

    assert isinstance(results1, dict)

    # SET2

    enr = Enrichment()

    enr.select_features(list2)

    enr.full_enrichment()

    results2 = None

    results2 = enr.get_results()

    assert isinstance(results1, dict)

    ans = Analysis(results2)

    ans.networks_metadata()

    ans.interactions_metadata()

    ans.set_p_value(value=0.05)
    ans.set_test(test="FISH")
    ans.set_correction(correction=None)

    ans.full_analysis()

    results2 = None

    results2 = ans.get_full_results()

    assert isinstance(results2, dict)

    dsa_compare = DSA(results1, results2)

    dsa_compare.GO_diff()

    tmp = None

    tmp = dsa_compare.get_GO_diff()

    assert isinstance(tmp, dict)

    dsa_compare.KEGG_diff()

    tmp = None

    tmp = dsa_compare.get_KEGG_diff()

    assert isinstance(tmp, dict)

    dsa_compare.REACTOME_diff()

    tmp = None

    tmp = dsa_compare.get_REACTOME_diff()

    assert isinstance(tmp, dict)

    dsa_compare.spec_diff()

    tmp = None

    tmp = dsa_compare.get_specificity_diff()

    assert isinstance(tmp, dict)

    dsa_compare.gi_diff()

    tmp = None

    tmp = dsa_compare.get_GI_diff()

    assert isinstance(tmp, dict)

    dsa_compare.network_diff()

    tmp = None

    tmp = dsa_compare.get_networks_diff()

    assert isinstance(tmp, dict)

    dsa_compare.inter_processes()

    tmp = None

    tmp = dsa_compare.get_inter_terms()

    assert isinstance(tmp, dict)

    dsa_compare.connections_diff()

    tmp = dsa_compare.get_set_to_set_con()

    assert isinstance(tmp, dict)

    del tmp

    dsa_compare.full_analysis()

    results3 = dsa_compare.get_results()

    vis_des = VisualizationDES(results3)

    plot = None

    plot = vis_des.diff_SPECIFICITY_plot(
        p_val=0.05,
        test="FISH",
        adj="BH",
        n=6,
        min_terms=1,
        selected_set=[],
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis_des.diff_gene_type_plot(
        set1_name="Set 1", set2_name="Set 2", image_width=12, image_high=6, font_size=15
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis_des.diff_GO_plot(
        p_val=0.05,
        test="FISH",
        adj="BH",
        n=25,
        min_terms=5,
        selected_parent=[],
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis_des.diff_KEGG_plot(
        p_val=0.05,
        test="FISH",
        adj="BH",
        n=25,
        min_terms=5,
        selected_parent=[],
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis_des.diff_REACTOME_plot(
        p_val=0.05,
        test="FISH",
        adj="BH",
        n=25,
        min_terms=5,
        selected_parent=[],
        width=10,
        bar_width=0.5,
        stat="n",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    plot = None

    plot = vis_des.diff_SPECIFICITY_plot(
        p_val=0.05,
        test="FISH",
        adj="BH",
        n=5,
        min_terms=1,
        selected_set=[],
        width=10,
        bar_width=0.5,
        stat="p_val",
    )

    assert isinstance(plot, figure.Figure)

    del plot

    go = None

    go = vis_des.diff_GOPa_network_create(
        data_set="GO-TERM",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )

    assert isinstance(go, nx.Graph)

    del go

    kegg = None

    kegg = vis_des.diff_GOPa_network_create(
        data_set="KEGG",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )

    assert isinstance(kegg, nx.Graph)

    del kegg

    reactome = None

    reactome = vis_des.diff_GOPa_network_create(
        data_set="REACTOME",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )

    assert isinstance(reactome, nx.Graph)

    del reactome

    gi = None

    gi = vis_des.diff_GI_network_create(min_con=2)

    assert isinstance(gi, nx.Graph)

    del gi

    aml = None

    aml = vis_des.diff_AUTO_ML_network(
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=False,
        selected_parents=[],
        selected_genes=[],
    )

    assert isinstance(aml, nx.Graph)

    del aml

    dgs = None

    dgs = vis_des.diff_gene_scatter(
        set_num=1,
        colors="viridis",
        species="human",
        hclust="complete",
        img_width=None,
        img_high=None,
        label_size=None,
        x_lab="Genes",
        legend_lab="log(TPM + 1)",
        selected_list=[],
    )

    assert isinstance(dgs, dict)

    del dgs
