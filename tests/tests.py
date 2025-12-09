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

