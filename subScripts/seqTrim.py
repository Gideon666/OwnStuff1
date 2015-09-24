#!/usr/bin/env python

"""
########################################################
# Special adapter for clipping should be given
# e.g. reverse seq of mirE+spacer TG
# CGCTCACTGTCAACAGCAATATACTG
# 
# or normal adapter P7-Sequence revers
# ATCTCGTATGCCGTCTTCTGCTTG is default
# P7-Sequence is clipped from sequence,
# sequence ist reverse complemented and
# first 8 Bases taken as barcode
# 
# to asure, that all Sequences start with the spacer and the mir30/E
# sequence
# usage:
    seqTrim.py [options] "fastq File"
    -o, --output        :outputfile
    -b, --barcodes      :barcodefile for barcode splitting
    -a, --adapter       :adapter for 3' adapter trimming
    -s, --silent        :silent mode, default verbose
    -a, --automatic     :called from script don't use it.
    -c, --configFile    :automatic initializing of options from configFile
    -d, --destination   :destination folder for output
"""

import sys, string, os, getopt, re
import datetime
######Level One##########

def processArgs(args, optionDir):
    """
    """
    # Process Barcode File
    if len(args) == 0:
        print __doc__
        sys.exit(1)

    if 'barcodeFile' in optionDir.keys():
        bcodes = []
        bcnames = []
        ihandle = open(optionDir['barcodeFile'], 'r')
        ibuffer = ihandle.readlines()
        ihandle.close()
        for line in ibuffer:
            if line.startswith('#'):
                continue
            bc_string = string.split(string.strip(line))
            bcodes.append(bc_string[1])
            if optionDir['verbose']:
                sys.stdout.write(bc_string[1]+"\n")
        optionDir['barcodes'] = bcodes
        sys.stdout.write('Barcodes: {0}\n'.format(str(bcodes)))
    # Process Adapter
    if 'adapter' in optionDir.keys():
        precompiledSet = []
        seq = string.upper(optionDir['adapter'])
        precompiledSet = precompAdapter(optionDir['adapter'], optionDir)
        #while(len(seq) > optionDir['adapterMinLen']):
        #    seq = seq[:len(seq)-1]
        #    precompiledSet.append(re.compile(seq))
        #print precompiledSet
        #sys.exit(0)
        optionDir['precompiledSet'] = precompiledSet
    if 'configFile' in optionDir.keys():
        barcodeInfo, commandInfo = readConfig(optionDir)
        optionDir['barcodes'] = barcodeInfo.keys()

    return args, optionDir

def manage(args, optionDir):
    """
    """
    buffer1 = None
    for a in args:
        buffer1 = readFastq(a, optionDir)
    saveBarcodedFastQ(buffer1[0], optionDir)
    return args, optionDir

######Level Two##########

def saveBarcodedFastQ(fastqDict, optionDir):
    """
    """
    print "save to files "
    print fastqDict.keys()
    for k, v in fastqDict.items():
        if optionDir['verbose']:
            sys.stdout.write(optionDir['oFile'][0]+"_BC_"+k+".fastq\n")
        if 'destination' in optionDir.keys():
            with open(optionDir['destination']+"/"+optionDir['oFile'][0]+"_BC_"+k+".fastq","w")\
                    as ohandle:
                for line in v:
                    ohandle.write(line+"\n")
        else: 
            with open(optionDir['oFile'][0]+"_BC_"+k+".fastq",'w') as ohandle:
                for line in v:
                    ohandle.write(line+"\n")
        ohandle.close()

def readFastq(fastqFile, optionDir):
    """
    """
    ihandle = open(fastqFile, 'r')
    trimBufferDict = {}
    BarcodeDict = {}
    statDict = {}	
    item = [None]*4
    item[0] = string.strip(ihandle.readline())
    counter = 1
    trimBufferDict['_OUT_'] = []

    for line in ihandle:
        if line.startswith("@M01950"):
            trimmedItem = trimItem(item, optionDir)
            #if trimmedItem[0] not in optionDir['barcodes']:
            #    print trimmedItem[0]
            if trimmedItem[0] == None:
                trimBufferDict['_OUT_'].extend(trimmedItem[1])
                counter = 0
                item = [None] * 4
                item[counter] = string.strip(line)
                counter += 1
                continue
            if trimBufferDict.has_key(trimmedItem[0]):
                trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                statDict[trimmedItem[0]] += 1
            else:
                #print "Key :"+str(trimmedItem[0])
                trimBufferDict[trimmedItem[0]] = []
                trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                statDict[trimmedItem[0]] = 1
            counter = 0
            item = [None]*4
        item[counter] = string.strip(line)
        counter += 1
    
    trimmedItem = trimItem(item, optionDir)
    if trimmedItem[0] == None:
        trimBufferDict['_OUT_'].extend(trimmedItem[1])
        counter = 0
        item = [None] *4
    else:
        if trimBufferDict.has_key(trimmedItem[0]):
            trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
            statDict[trimmedItem[0]] += 1
        else:
            trimBufferDict[trimmedItem[0]] = []
            trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
            statDict[trimmedItem[0]] = 1
    ihandle.close()
    if optionDir['verbose']:
        sys.stdout.write("------------\n")
        for k, v in trimBufferDict.items():
            sys.stdout.write("Number of seq in {0} : {1}\n".format(k, len(v)/4))
    return trimBufferDict, statDict, optionDir

def trimItem(item, optionDir, trimNumber = 24, barcodeNumber = 8):
    """
    """
               #              2+Spacer:Barcode:+4
    emPatSet = ['ACAGCAATATACTG','ACTG[\w]{8}ATCT']

    if "adapter" in optionDir.keys():
        adapterSeq = string.strip(optionDir['adapter'])
        patSet = optionDir['precompiledSet']
        #sys.stdout.write("Adapter Found! : length : {0}\n".format(len(patSet)))
    else:
        #adapterSeq = "AATATACTG"
        # mir30 + spacer TG
        # adapterSeq = "ATCTCGTATGCCGTCTTCTGCTTG"
        # careful adjust trimNumber accordingly
        
        # adapter P7
        #adapterSeq = 'ATCTCGTATGCCGTCTTCTGCTTG'
        adapterSeq =  'ATCTCGTATGCC'
        patSet = precompAdapter(adapterSeq, optionDir)

    retItem = list(item)
    seq = item[1]
    Barcode = None
    #for testing reduced patSet
    for pattern in patSet[len(patSet)-1:]:
        if pattern.search(seq):
            trimNumber = pattern.search(seq).start()
            retItem[1] = item[1][:trimNumber]
            retItem[3] = item[3][:trimNumber]
            Barcode = "Set"
            #sys.stdout.write("Pattern 1\n")
            break
    if(not(Barcode)):
        emPattern1 = re.compile(emPatSet[0])
        if emPattern1.search(seq):
            trimNumber = emPattern1.search(seq).end()+8
            retItem[1] = item[1][:trimNumber]
            retItem[3] = item[3][:trimNumber]
            Barcode = "Set"
            #sys.stdout.write("EM-Pattern 1\n")
    if(not(Barcode)):
        emPattern2 = re.compile(emPatSet[1])
        if emPattern2.search(seq):
            trimNumber = emPattern2.search(seq).end()-4
            retItem[1] = item[1][:trimNumber]
            retItem[3] = item[3][:trimNumber]
            Barcode = "Set"
            #sys.stdout.write("EM-Pattern 2\n")
    Barcode2 = reverseComplement(retItem[1][-8:])
    retItem[1] = reverseComplement(retItem[1])
    retItem[3] = retItem[3][::-1]
    item[1] = reverseComplement(item[1])
    item[3] = item[3][::-1]
    if Barcode:
        Barcode = retItem[1][:8]
        #Barcode2 = reverseComplement(retItem[1][-8:])
        #Debugger
        #if (Barcode2 != Barcode):
        #    print "------------"
        #    print retItem
        #    print "Barcode: "+Barcode
        #    print "Barcode2: "+Barcode2
        #    exit(-1)
        if 'barcodes' in optionDir.keys():
            #print "---"
            #print optionDir['barcodes']
            #print "---"
            if Barcode in optionDir['barcodes']:
                #print Barcode
                #print optionDir['barcodes']
                #print "----"
                return Barcode, retItem
            else:
                Barcode = sameKey(Barcode, optionDir['barcodes'])
                if Barcode:
                    #print Barcode
                    #print optionDir['barcodes']
                    #print Barcode in optionDir['barcodes']
                    #print "----"
                    return Barcode, retItem
    return None, item
    
######Level three####################

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


def sameKey(cKey, keys2Check):
    """
    """
    nCheck = re.compile('N')
    counter = 0
    swList = []
    if nCheck.search(cKey):
        if 1 < len(nCheck.findall(cKey)):
            return None #cKey
        else:
            
            arepl = nCheck.sub("A", cKey)
            grepl = nCheck.sub("G", cKey)
            crepl = nCheck.sub("C", cKey)
            trepl = nCheck.sub("T", cKey)
            swList = [arepl, grepl, crepl, trepl]
            if (arepl in keys2Check):
	        counter = 1
            if ((grepl in keys2Check) and (counter == 0)):
                counter = 2
            if ((crepl in keys2Check) and (counter == 0)):
                counter = 3
            if ((trepl in keys2Check) and (counter == 0)):
                counter = 4
            if (counter == 0):
                return None
            else:
                return swList[counter-1]
    else:
        return None
        
def precompAdapter(adapterSeq, optionDir, reverse=True):
    """
    precompiles the adapter into 
    regex up to the minimum length of adapterMinLen
    returns a list of the re objects
    """
    sys.stdout.write("calculating ....{0}\n".format(adapterSeq))
    retList = []
    seq = adapterSeq
    while(len(seq) > optionDir['adapterMinLen']):
        if reverse:
            seq = seq[:len(seq)-1]
            retList.append(re.compile(seq))
            sys.stdout.write("calculating ....{0}\n".format(seq))
        else:
            seq = seq[1:]
            retList.append(re.compile(seq))
    return retList

def readConfig(optionDir):
    """
    # comments
    @ special information
    @ Header SampleNum   Barcode PoolId  Description     Groups(x to ignore)
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
                #print line
                atoms = string.split(string.strip(line))
                barcodes[atoms[1]] = atoms
    return barcodes, commands 

######Defaults###########

def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['temp','txt'],\
                'verbose': True,\
                'adapter': 'ATCTCGTATGCCGTCTTCTGCTTG',\
                'adapterMinLen' : 10}
    return ValueDir

#######Main#######

def main():
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, b:, a:, s, d:, t, c:",\
                ["help", "output=", "barcodes=", "adapter=", "silent", "destination=",\
                "automatic", "configFile" ])
    except getopt.Error, msg:
        sys.stdout.write(msg + "\n")
        sys.stdout.write("For help, use -h!\n")
        sys.exit(-1)

    for o, a in opts:
        if o in ['-h', '--help']:
            print __doc__
            sys.exit(1)
        if o in ['-o', '--output']:
            optionDir['oFile'][0] = string.split(os.path.basename(string.strip(a)), sep=".")[0]
            optionDir['oFile'][0] +=\
            datetime.datetime.now().strftime('_%Hh_%d_%m_%y')
            if optionDir['verbose']:
                sys.stdout.write(optionDir['oFile'][0]+"\n")
        if o in ['-d', '--destination']:
            optionDir['destination'] = string.strip(a)
        if o in ['-b', '--barcodes']:
            optionDir['barcodeFile'] = string.strip(a)
        if o in ['-a', '--adapter']:
            optionDir['adapter'] = string.strip(a)
        if o in ['-s', '--silent']:
            #print "set silent mode?"
            optionDir['verbose'] = False
        if o in ['-t', '--automatic']:
            # if called from script
            optionDir['automatic'] = "auto"
        if o in ['-c', '--configFile']:
            optionDir['configFile'] = string.strip(a)
    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()


