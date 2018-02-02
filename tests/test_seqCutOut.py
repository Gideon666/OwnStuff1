#! /usr/bin/env python

import pytest
from subScripts import seqCutOut as SCO


@pytest.fixture
def crsprSeq():
    """Creates a Crspr seq for cutting
    without errors
    crspr = >Acaca|T28_cd|834.32477646|-Pool:ALL
    'GCATGCGATCTATCCGTCGG'
    """
    seq = "TGTCCTCAGTATATTGCTGTTGACAGTGAGCGCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGCATGCGATCTATCCGTCGGGTTTTAGAGCTAGAAATAGCAAGTTAAAA"
    return seq

@pytest.fixture
def optionDir():
    """standardValue optionDir
    """
    Options = SCO.setDefaultValues()
    return Options

@pytest.fixture
def optionDirCrspr():
    """
    """
    Options = SCO.setDefaultValues()
    Options['patternCutLeft'] = -38
    Options['patternCutRight'] = 28
    return Options

class Test_regExLoop(object):

    def test_standard_matching(self, crsprSeq, optionDirCrspr):
        assert(SCO.regExLoop(crsprSeq,\
                "CTTGGCTTTATATATCTTGTGGAAAGGACG",\
                70, 90, optionDirCrspr) ==\
                'GCATGCGATCTATCCGTCGG')




class Test_pswmLoop(object):

    def test_standard_matching(self, crsprSeq, optionDirCrspr):
        optionDirCrspr['loopMM'] = 0
        assert(SCO.pswmLoop(crsprSeq,\
                "CTTGGCTTTATATATCTTGTGGAAAGGACG",\
                -38,\
                28,\
                optionDirCrspr\
                ) == 'GCATGCGATCTATCCGTCGG')

        
    def test_mm_matching(self, crsprSeq, optionDirCrspr):
        optionDirCrspr['loopMM'] = 2 # changed TT to AA
        assert(SCO.pswmLoop(crsprSeq,\
                "CTTGGCTTTATATATC"+\
                "AA"+\
                "GTGGAAAGGACG",\
                -38,\
                28,\
                optionDirCrspr\
                ) == 'GCATGCGATCTATCCGTCGG')





