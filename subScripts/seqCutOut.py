#!/usr/bin/env python

import sys, string, os, getopt
import datetime
import re

"""
no longer turns the sequences arround (reverse complement)!!!!
necessary to have the same startpoint (barcode?)

does not cut the quality string!
now produce fasta instead of fastq


usage: seqCutout.py [options] inputFiles
options:
    -m --mode: possible modes   n= normal
                                i= iTags
"""
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
    resDict = {} 
    for a in args:
        resDict[string.split(a, sep=".")[0]] = cutOutFastq(a, optionDir)

    for k, v in resDict.items():
        #saveCutFastq(k, v, optionDir)
        saveCutFasta(k, v, optionDir)



######Level Two###########

def cutOutFastq(fileFastq, optionDir):
    """
    new cutOut Function for Fasta file format
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

def smartCutting(seq, startPos, endPos, optionDir):
    """
    """
    lS=22
    rS=22
    retSeq = ""
    mLen = optionDir['minLength']
    loop = "TAGTGAAGCCACAGATGTA"
    looP = re.compile(loop)
    reL = looP.search(seq)
    if(reL):
        retSeq = seq[reL.start()-lS:reL.end()+rS]
    else:
        retSeq = seq[startPos:endPos]
    if(len(retSeq) < mLen):
        return None
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

def saveCutFasta(filename, resList, optionDir):
    """
    """
    oname = filename + optionDir['oFile'][0] +"."+ optionDir['oFile'][1]
    with open(oname, 'w') as ohandle:
        for line in resList:
            ohandle.write(line+"\n")
    ohandle.close()

def saveCutFastq(filename, reslist, optionDir):
    """
    """
    oname = filename + optionDir['oFile'][0] +"."+ optionDir['oFile'][1]

    ohandle = open(oname, 'w')
    for line in reslist:
        ohandle.write(line+"\n")
    ohandle.close()

####Level three#######################

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

def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['_cutout','fasta'],\
                'mode': "N",\
                'minLength' : 63,\
                'startCutout' : 34,\
                'endCutout' : 97\
               }
    return ValueDir

#######Main#######

def main():
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, m:, c:",\
                ["help", "output=", "mode=", "configFile="])
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
        if o in ['-m', '--mode']:
            mode = string.upper(string.strip(a))
            if mode in ['N', 'I']: #modes
                optionDir['mode'] = mode
        if o in ['-c', '--configFile']:
            optionDir['configFile'] = string.strip(a)

    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()


