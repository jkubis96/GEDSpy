import pytest

from gedspy import GetDataRaw



def test_GetDataRaw():

    gdr = GetDataRaw()

    recteome = None

    recteome = gdr.get_raw_REACTOME()

    assert  isinstance(recteome, dict)

    del recteome



