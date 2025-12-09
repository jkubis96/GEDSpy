import matplotlib.figure as fig
import networkx as nx
import pytest

from gedspy import GetDataRaw


def test_GetDataRaw():

    gdr = GetDataRaw()

    recteome = None

    recteome = gdr.get_raw_REACTOME()

    assert isinstance(recteome, dict)




###############################################################################

from gedspy import GetData


def test_GetDataRaw():

    gd = GetData()

    reactome = None

    reactome = gd.get_REACTOME()

    assert isinstance(reactome, dict)

    


############################################################################

