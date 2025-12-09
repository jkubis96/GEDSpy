import matplotlib.figure as fig
import networkx as nx
import pytest

from gedspy import GetDataRaw
from gedspy import GetData
from gedspy import Enrichment


def test_GetDataRaw():

    gdr = GetDataRaw()

    recteome = None

    recteome = gdr.get_raw_REACTOME()

    assert isinstance(recteome, dict)




def test_GetDataRaw():

    gd = GetData()

    reactome = None

    reactome = gd.get_REACTOME()

    assert isinstance(reactome, dict)

    




def Enrichment_tests():

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


    enr = Enrichment()

    enr.select_features(list1)

    gene_info = enr.get_gene_info

    assert isinstance(gene_info, dict)

    enr.enriche_specificiti()

    HPA = enr.get_HPA

    assert isinstance(HPA, dict)

    enr.enriche_KEGG()

    KEGG = enr.get_KEGG

    assert isinstance(KEGG, dict)

    enr.enriche_GOTERM()

    GOTERM = enr.get_GO_TERM

    assert isinstance(GOTERM, dict)

    enr.enriche_REACTOME()

    REACTOME = enr.get_REACTOME

    assert isinstance(REACTOME, dict)

    enr.enriche_DISEASES()

    DISEASES = enr.get_DISEASES

    assert isinstance(DISEASES, dict)

    enr.enriche_ViMIC()

    ViMIC = enr.get_ViMIC

    assert isinstance(ViMIC, dict)

    enr.enriche_IntAct()

    IntAct = enr.get_IntAct

    assert isinstance(IntAct, dict)

    enr.enriche_STRING()

    STRING = enr.get_STRING

    assert isinstance(STRING, dict)

    enr.enriche_CellCon()

    CellConnections = enr.get_CellCon

    assert isinstance(CellConnections, dict)

    enr.enriche_RNA_SEQ()

    RNASEQ = enr.get_RNA_SEQ

    assert isinstance(RNASEQ, dict)

    del enr



