#!/usr/bin/env python
"""
no longer turns the sequences arround (reverse complement)!!!!
necessary to have the same startpoint (barcode?)

does not cut the quality string!
now produce fasta instead of fastq
Needs the new regex module instead of the simple re for
fuzzy pattern matching!

usage: seqCutout.py [options] inputFiles
options:
    -m --mode: possible modes   n= normal
                                i= iTags
                                c= crspr
"""

import argparse
import sys, string, os, getopt
import datetime
import regex
import aux1.pswm as ps
import re
import pdb


######Level One##########

def processArgs(args, optionDir):
    """
    process options and arguments.
    """
    if 'configFile' in optionDir.keys():
        barcodeInfo, commandInfo = readConfig(optionDir)
        optionDir['barcodeInfo'] = barcodeInfo
        optionDir['commandInfo'] = commandInfo

    return args, optionDir

def manage(args, optionDir):
    """
    """
    #pdb.set_trace()
    resDict = {} 
    for a in args:
        resDict[string.split(a, sep=".")[0]] = cutOutFastq(a, optionDir)
    for k, v in resDict.items():
        #saveCutFastq(k, v, optionDir)
        saveCutFasta(k, v, optionDir)


######Level Two###########

#manage
def cutOutFastq(fileFastq, optionDir):
    """
    new cutOut Function for Fasta file format
    manager cutoutFunction
    """
    #### Initialize
    cuttedList = []
    fastaElements = [None]*4
    retElement = [None]*2
    lineCounter = 0
    startPos = optionDir['startCutout']
    endPos = optionDir['endCutout']
    with open(fileFastq, 'r') as ihandle:
        for line in ihandle:
            entryLine = string.strip(line)
            if lineCounter == 4:
                #dosomething
                retElement[0] = ">"+fastaElements[0]
                retElement[1] = smartCutting(fastaElements[1], startPos, endPos, optionDir)
                if(retElement[1]):
                    cuttedList.extend(retElement)

                lineCounter = 0
                fastaElements = [None]*4
                retElement = [None]*2
                
            fastaElements[lineCounter] = entryLine
            lineCounter += 1
    retElement[0] = ">"+fastaElements[0]
    retElement[1] = smartCutting(fastaElements[1], startPos, endPos, optionDir)
    if(retElement[1]):
        cuttedList.extend(retElement)
    return cuttedList

#manage
def saveCutFasta(filename, resList, optionDir):
    """
    """
    oname = filename + optionDir['oFile'][0] +"."+ optionDir['oFile'][1]
    with open(oname, 'w') as ohandle:
        for line in resList:
            ohandle.write(line+"\n")

#manage
def saveCutFastq(filename, reslist, optionDir):
    """
    """
    oname = filename + optionDir['oFile'][0] +"."+ optionDir['oFile'][1]
    with open(oname, 'w') as ohandle:
        for line in reslist:
            ohandle.write(line+"\n")

####Level three#######################

def smartCutting(seq, startPos, endPos, optionDir):
    """
    smartCutting matches the ANCHOR(loop for shRNAs) to the seq and 
    cuts left and right.
    The borders are hard (+22bp left and right)
    If allowed mismatches: pswms are used for the
    matching. Otherwise RegEx!

    For itags the itag is cutted from the end and
    glued back on, after cutting.

    If maskLoop is set the found loop is exchanged
    with the sequence for a perfect loop.
    """
    if optionDir['mode'] in ['N', 'I']:
        #startPos = 22
        #endPos = 22
        return smartCuttingShRNA(seq, startPos, endPos, optionDir)
    if optionDir['mode'] in ['C']:
        optionDir['patternCutLeft'] = -38
        optionDir['patternCutRight'] = 28
        startPos = 70
        endPos = 90
        return smartCuttingCrspr(seq, startPos, endPos, optionDir)

####Level four########################

def smartCuttingCrspr(seq, startPos, endPos, optionDir):
    """
    crspr Mode smartCutting
    """
    ## Init for Crspr
    anchor = optionDir['ANCHOR']
    retSeq = pswmLoop(seq, anchor, startPos, endPos, optionDir)
    return retSeq    

#cutOutFastq
def smartCuttingShRNA(seq, startPos, endPos, optionDir):
    """
    shRNA Mode smartCutting
    care regex cant find gabIns!!!
    """
    lS=22
    rS=22
    retSeq = ""
    mLen = optionDir['minLength']
    loop = "TAGTGAAGCCACAGATGTA"
    if (optionDir['gabIns'] != 0):
        retSeq = regExLoop(seq, loop, startPos, endPos, optionDir)
    else:
        retSeq = pswmLoop(seq, loop, startPos, endPos, optionDir)
    if optionDir['maskLoop']:
        retSeq = retSeq[:lS]+loop+retSeq[lS+loop:]
    ####
    if(len(retSeq) < mLen):
        return None
    ####itags
    if(optionDir['mode'] == 'I'):
        iTagId = re.compile("GAATT")
        reO = iTagId.search(seq[-20:])
        if(reO):
            patStart = reO.end()
            iTag = seq[-20:][patStart:patStart+7]
        elif(seq[-8]=="T"):
            iTag = seq[-7:]
        elif(seq[-11:-8] == "AAT"):
            iTag = seq[-7:]
        else:
            iTag = "NNNNNNN"
        while len(iTag) < 7:
            iTag+="N"
        retSeq+=iTag
    return retSeq

####Aux##############

#smartCutting
def pswmLoop(seq, loop, startPos, endPos, optionDir = {}):
    """
    """
    minScore = len(loop) - optionDir['loopMM']
    lS = optionDir['patternCutLeft']
    rS = optionDir['patternCutRight']
    #minScore = 18.0
    #print minScore
    #print loop
    loopmask = ps.PSWM([loop])
    sliceList = slwApproach(seq, len(loop))
    for pos in xrange(len(sliceList)):
        if loopmask.getScore(sliceList[pos]) >= minScore:
            #if len(seq[pos-lS:pos+len(loop)+rS]) == 63:
                #print loopmask.getScore(sliceList[pos])
                #print seq 
                #print len(seq)
                #print seq[pos-lS:pos+len(loop)+rS]
                #print len(seq[pos-lS:pos+len(loop)+rS])
            #print "found"
            return seq[max(0, pos-lS):pos+len(loop)+rS]
    return seq[startPos:endPos]

#pswmLoop
def slwApproach(seq, size, hard=True):
    """
    sliding window approach for strings
    """
    retList = []
    if hard:
        for i in xrange(len(seq)-size+1):
            retList.append(seq[i:i+size])
    else:
        pass
    return retList

#smartCutting
def regExLoop(seq, loop, startPos, endPos, optionDir = {}):
    """
    regEx Error, cant find gabIns
    """
    lS=optionDir['patternCutLeft']
    rS=optionDir['patternCutRight']
    inDels = optionDir['gabIns']
    looP = re.compile('{0}'.format(loop))
    #looP = re.compile(r'(?:{0}){{e<={1}}}'.format(loop, inDels))
    reL = looP.search(seq)
    if(reL):
        retSeq = seq[max(0, reL.start()-lS):reL.end()+rS]
    else:
        if optionDir['maskLoop']:
            nonMasked = seq[startPos:endPos]
            retSeq = retSeq[:lS]+loop+retSeq[lS+len(loop):]
            print "-Mask--"
            print nonMasked
            print retSeq
        else:    
            retSeq = seq[startPos:endPos]
        #sys.stdout.write("Warning in regExLoop!\nPattern not found!\n")
    return retSeq

def reverseComplement(seq):
        """ 
        """
        retSeq = ""
        swDict = {"A" : "T",\
                  "T" : "A",\
                  "G" : "C",\
                  "C" : "G",\
                  "N" : "N",\
                  "a" : "t",\
                  "t" : "a",\
                  "g" : "c",\
                  "c" : "g",\
                  "n" : "n"}
        revSeq = seq[::-1]
        for i in revSeq:
            retSeq += swDict[i]
        return retSeq

def readConfig(optionDir):
    """
    # comments
    @ special information
    @Header SampleNum   Barcode PoolId  Description     Groups(x to ignore)
    """
    commands = {}
    barcodes = {}
    with open(optionDir['configFile'], 'r') as ihandle:
        for line in ihandle:
            if line.startswith("#"):
                continue
            elif line.startswith("@"):
                continue
            else:
                #sys.stdout.write("ConfigFile read : \n\t\t{0}\n".format(string.strip(line)))
                atoms = string.split(string.strip(line))
                barcodes[atoms[1]] = atoms
    return barcodes, commands 

######Defaults###########

def translate_args_options(args_opts):
    opts = vars(args_opts)
    args = opts['args']
    del(opts['args'])

    return args, opts


def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['_cutout','fasta'],\
                'mode': "N",\
                'maskLoop': False,\
                'loopMM': 2,\
                'minLength':63,\
                'startCutout' : 34,\
                'endCutout' : 97,\
                'patternCutLeft' : 22,\
                'patternCutRight' : 22,\
                'gabIns' : 0,\
                'ANCHOR': 'TAGTGAAGCCACAGATGTA'\
                }
    return ValueDir

#######Main#######


def main2():
    ###var init
    global optionDir
    optionDir = setDefaultValues()
    ###
    parser = argparse.ArgumentParser(description=__doc__,\
            prog=string.split(sys.argv[0], sep=".").upper())

    parser.add_argument('-o', '--output', default=sys.stdout,\
            metavar="OFILE",\
            type=argparse.FileType('w'),\
            help="Path to output File")

    args_opts = parser.parse_args()
    args, opts = translate_args_options(args_opts)
    manage(args, optionDir)


def main():
    global optionDir
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, m:, c:, l:, a, g:, r:, L:, R:",\
                ["help", "output=", "mode=", "configFile=",\
                 "loopMissMatches=", "maskLoop", "gabIns=",\
                 "anchor=", "patternCutLeft=", "patternCutRight="])
    except getopt.Error, msg:
        sys.stdout.write(msg + "\n")
        sys.stdout.write("For help, use -h!\n")
        sys.exit(-1)
    for o, a in opts:
        if o in ['-h', '--help']:
            print __doc__
            sys.exit(1)
        if o in ['-o', '--output']:
            optionDir['oFile'][0] = string.split(string.strip(a), sep=".")[0]
            optionDir['oFile'][0] +=\
                    datetime.datetime.now().strftime('_%Hh_%d_%m_%y')
        if o in ['-L', '--patternCutLeft']:
            optionDir['patternCutLeft'] = int(string.strip(a))
        if o in ['-R', '--patternCutRight']:
            optionDir['patternCutRight'] = int(string.strip(a))
        if o in ['-m', '--mode']:
            mode = string.upper(string.strip(a))
            if mode in ['N', 'I', 'C']: #modes
                optionDir['mode'] = mode
        if o in ['-c', '--configFile']:
            optionDir['configFile'] = string.strip(a)
        if o in ['-l', '--loopMissMatches']:
            optionDir['loopMM'] = string.atoi(string.strip(a))
        if o in ['-a', '--maskLoop']:
            optionDir['maskLoop'] = True
        if o in ['-r', '--anchor']:
            optionDir['ANCHOR'] = string.strip(a).upper()
            sys.stdout.write("Anchor for cutout, set to: {0}\n".format(optionDir['ANCHOR']))
        if o in ['-g', '--gabIns']:
            print "gabs:"+str(a)
            optionDir['gabIns'] = string.atoi(string.strip(a))
    


    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()


