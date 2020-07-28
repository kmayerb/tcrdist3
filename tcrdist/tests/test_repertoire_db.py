import pytest
import numpy as np
from tcrdist import repertoire_db

def test_extract_cdrs_without_gaps():
    """ specific case """
    a = repertoire_db._extract_cdrs_without_gaps('GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
    '28-39;57-66;82-88;106-111')
    b = ['TSGFNG', 'VVLDGL', 'SRSNGY', 'CAVR']
    
    assert a==b

def test_extract_aligned_cdrs():
    """ specific case """
    a = repertoire_db._extract_aligned_cdrs('GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
    '28-39;57-66;82-88;106-111')
    b = ['TSG......FNG', 'VVL....DGL', 'SRSN.GY', 'CAVR..']
    
    assert a==b

