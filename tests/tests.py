import matplotlib.figure as fig
import networkx as nx
import pytest

from gedspy import GetDataRaw


def test_GetDataRaw():

    gdr = GetDataRaw()

    recteome = None

    recteome = gdr.get_raw_REACTOME()

    assert isinstance(recteome, dict)

    del recteome

    ref_gen = None

    ref_gen = gdr.get_raw_REF_GEN()

    assert isinstance(ref_gen, dict)

    del ref_gen

    ref_gen_seq = None

    ref_gen_seq = gdr.get_raw_REF_GEN_RNA_SEQ()

    assert isinstance(ref_gen_seq, dict)

    del ref_gen_seq

    HPA = None

    HPA = gdr.get_raw_HPA()

    assert isinstance(HPA, dict)

    del HPA

    diseases = None

    diseases = gdr.get_raw_DISEASES()

    assert isinstance(diseases, dict)

    del diseases

    vimic = None

    vimic = gdr.get_raw_ViMIC()

    assert isinstance(vimic, dict)

    del vimic

    kegg = None

    kegg = gdr.get_raw_KEGG()

    assert isinstance(kegg, dict)

    del kegg

    go = None

    go = gdr.get_raw_GO()

    assert isinstance(go, dict)

    del go

    intact = None

    intact = gdr.get_raw_IntAct()

    assert isinstance(intact, dict)

    del intact

    string = None

    string = gdr.get_raw_STRING()

    assert isinstance(string, dict)

    del string

    cell_talk = None

    cell_talk = gdr.get_raw_CellTalk()

    assert isinstance(cell_talk, dict)

    del cell_talk

    cell_phone = None

    cell_phone = gdr.get_raw_CellPhone()

    assert isinstance(cell_phone, dict)

    del cell_phone


###############################################################################

from gedspy import GetData


def test_GetDataRaw():

    gd = GetData()

    reactome = None

    reactome = gd.get_REACTOME()

    assert isinstance(reactome, dict)

    ref_gen = None

    ref_gen = gd.get_REF_GEN()

    assert isinstance(ref_gen, dict)

    ref_gen_seq = None

    ref_gen_seq = gd.get_REF_GEN_RNA_SEQ()

    assert isinstance(ref_gen_seq, dict)

    HPA = None

    HPA = gd.get_HPA()

    assert isinstance(HPA, dict)

    disease = None

    disease = gd.get_DISEASES()

    assert isinstance(disease, dict)

    vimic = None

    vimic = gd.get_ViMIC()

    assert isinstance(vimic, dict)

    kegg = None

    kegg = gd.get_KEGG()

    assert isinstance(kegg, dict)

    go = None

    go = gd.get_GO()

    assert isinstance(go, dict)

    intact = None

    intact = gd.get_IntAct()

    assert isinstance(intact, dict)

    string = None

    string = gd.get_STRING()

    assert isinstance(string, dict)

    CellTalk = None

    CellTalk = gd.get_CellTalk()

    assert isinstance(CellTalk, dict)

    CellPhone = None

    CellPhone = gd.get_CellPhone()

    assert isinstance(CellPhone, dict)

    CellInteractions = None

    CellInteractions = gd.get_interactions()

    assert isinstance(CellInteractions, dict)
