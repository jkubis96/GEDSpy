import matplotlib.figure as fig
import networkx as nx
import pytest

from gedspy import (
    DSA,
    Analysis,
    Enrichment,
    GetData,
    GetDataRaw,
    Visualization,
    VisualizationDES,
)


@pytest.fixture
def gdr():
    return GetDataRaw()


@pytest.mark.parametrize(
    "method_name",
    [
        "get_raw_REACTOME",
        "get_raw_REF_GEN",
        "get_raw_REF_GEN_RNA_SEQ",
        "get_raw_HPA",
        "get_raw_DISEASES",
        "get_raw_ViMIC",
        "get_raw_KEGG",
        "get_raw_GO",
        "get_raw_IntAct",
        "get_raw_STRING",
        "get_raw_CellTalk",
        "get_raw_CellPhone",
    ],
)
def test_GetDataRaw_methods_return_dict(gdr, method_name):
    """Test that all GetDataRaw.* methods return a dict."""
    method = getattr(gdr, method_name)
    result = method()

    assert isinstance(
        result, dict
    ), f"{method_name} should return dict, got {type(result)}"


###############################################################################


@pytest.fixture
def gd():
    return GetData()


@pytest.mark.parametrize(
    "method_name",
    [
        "get_REACTOME",
        "get_REF_GEN",
        "get_REF_GEN_RNA_SEQ",
        "get_HPA",
        "get_DISEASES",
        "get_ViMIC",
        "get_KEGG",
        "get_GO",
        "get_IntAct",
        "get_STRING",
        "get_CellTalk",
        "get_CellPhone",
        "get_interactions",
    ],
)
def test_GetData_methods_return_dict(gd, method_name):
    """Test that all GetData.* methods return a dict."""
    method = getattr(gd, method_name)
    result = method()

    assert isinstance(
        result, dict
    ), f"{method_name} should return dict, got {type(result)}"

