import pytest
"""
Unit Tests for translation.py
"""
def test_reverse_compliment():
    from tcrdist.translation import reverse_complement
    assert reverse_complement(seq = 'ATGC') == 'GCAT'
    assert reverse_complement(seq = 'AT.GC') == 'GC.AT'

