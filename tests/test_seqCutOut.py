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
def sensorSeq():
    """
    >shDnmt3a.1171Spl|Dnmt3a|Dnmt3aPlus_Sensor-Pool:ALL
    """
    seq = "TGCTGTTGACAGTGAGCGAACTCCAGATGTTCTTTGCCAATAGTGAAGCCACAGAT"+\
          "GTATTGGCAAAGAACATCTGGAGTCTGCCGAAGAGCGCTCTTCTATCCTTCTCGAC"+\
          "TCCAGATGTTCTTTGCCAATAACCATGGTCCAGGCAGGTGAATTCAGCATGAGTAT"+\
          "TTCCGCGA"
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

@pytest.fixture
def optionDirSensor():
    """
    """
    Options = SCO.setDefaultValues()
    Options['patternCutLeft'] = -22
    Options['patternCutRight'] = 34
    return Options

class Test_regExLoop(object):

    def test_standard_matching(self, crsprSeq, optionDirCrspr):
        assert(SCO.regExLoop(crsprSeq,\
                "CTTGGCTTTATATATCTTGTGGAAAGGACG",\
                70, 90, optionDirCrspr) ==\
                'GCATGCGATCTATCCGTCGG')

    def test_hairpin_matching(self, sensorSeq, optionDirSensor):
        test_sensor = "CTTCTCGACTCCAGATGTTCTTTGCCAATAACCA"
        anchor = "TGCCGAAGAGCGCTCTTCTATC"
        returned_seq = SCO.regExLoop(sensorSeq, anchor, 116, 150, optionDirSensor)
        print returned_seq
        assert(returned_seq == test_sensor)
        #assert 0

class Test_pswmLoop(object):

    def test_standard_matching_crspr(self, crsprSeq, optionDirCrspr):
        optionDirCrspr['loopMM'] = 0
        assert(SCO.pswmLoop(crsprSeq,\
                "CTTGGCTTTATATATCTTGTGGAAAGGACG",\
                -38,\
                28,\
                optionDirCrspr\
                ) == 'GCATGCGATCTATCCGTCGG')

        
    def test_mm_matching_crspr(self, crsprSeq, optionDirCrspr):
        optionDirCrspr['loopMM'] = 2 # changed TT to AA
        assert(SCO.pswmLoop(crsprSeq,\
                "CTTGGCTTTATATATC"+\
                "AA"+\
                "GTGGAAAGGACG",\
                -38,\
                28,\
                optionDirCrspr\
                ) == 'GCATGCGATCTATCCGTCGG')




