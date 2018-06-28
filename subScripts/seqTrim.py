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

import sys, string, os, getopt, re, shelve
import datetime
import pdb
from aux1.recSize import total_size as getsizeof
from bMod.Sequences import SequenceCompare as SeqComp

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
            if line.startswith('#') or string.strip(line) == "":
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
        if optionDir['MODE'] == 'crspr':
            precompiledSet = precompAdapter(optionDir['adapter'], optionDir, clipFromTail=False)#, reverse=False)
        else:
            precompiledSet = precompAdapter(optionDir['adapter'], optionDir)#, reverse=False)
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
    #close shelve
    #buffer1[0].close()
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
            with open(optionDir['destination']+"/"+optionDir['oFile'][0]+"_BC_"+k+".fastq","a")\
                    as ohandle:
                for line in v:
                    ohandle.write(line+"\n")
        else: 
            with open(optionDir['oFile'][0]+"_BC_"+k+".fastq",'a') as ohandle:
                for line in v:
                    ohandle.write(line+"\n")
        ohandle.close()

def readFastq(fastqFile, optionDir):
    """
    """
    print optionDir['MODE']
    #trimBufferDict = shelve.open("trimBufferDict", writeback=True)
    trimBufferKeys = []
    trimBufferDict = {}
    BarcodeDict = {}
    #statDict = {}	
    item = [None]*4
    counter = 1
    loopCounter = 0
    MacCode = optionDir['seqMachineId']
    trimBufferDict['_OUT_'] = []
    with open(fastqFile, 'r') as ihandle:
        item[0] = string.strip(ihandle.readline())
        MacCode = mac_code_check(item[0], MacCode)
        for line in ihandle:
            #check dictionary size
            if loopCounter >= 100000:
                loopCounter = 0
                dictsize = getsizeof(trimBufferDict)
                if dictsize >= optionDir['maxChunkSize']:
                    if optionDir['verbose']:
                        sys.stdout.write("Size of trimBufferDict : {0}\n".format(str(dictsize)))
                        sys.stdout.write("Saving Data...\n")
                    saveBarcodedFastQ(trimBufferDict, optionDir)
                    trimBufferDict = {}
                    trimBufferDict['_OUT_'] = []
                    #statDict = {}
            if line.startswith(MacCode):
                loopCounter +=1
                trimBufferDict = add_to_trimBufferDict(item, trimBufferDict, optionDir)
                ## new FastQ Item, Item[0] = Barcode, Item[1] = FastQ
                #trimmedItem = trimItem(item, optionDir)
                #if trimmedItem[0] not in optionDir['barcodes']:
                #if trimmedItem[0] == None:
                #    trimBufferDict['_OUT_'].extend(trimmedItem[1])
                #    counter = 0
                #    item = [None] * 4
                #    item[counter] = string.strip(line)
                #    counter += 1
                #    continue
                # change from normal directory to shelve
                #if trimBufferDict.has_key(trimmedItem[0]):
                #    trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                #    #statDict[trimmedItem[0]] += 1
                #else:
                #    sys.stdout.write("Key :{0} added!\n".format(str(trimmedItem[0])))
                #    trimBufferDict[trimmedItem[0]] = []
                #    trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                    #statDict[trimmedItem[0]] = 1
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
            # again the change
            if trimBufferDict.has_key(trimmedItem[0]):
                trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                #statDict[trimmedItem[0]] += 1
            else:
                trimBufferDict[trimmedItem[0]] = []
                trimBufferDict[trimmedItem[0]].extend(trimmedItem[1])
                #statDict[trimmedItem[0]] = 1
                
    ##### read_in terminated! #####
    if optionDir['verbose']:
        sys.stdout.write("------------\n")
        for k, v in trimBufferDict.items():
            sys.stdout.write("Number of seq in {0} : {1}\n".format(k, len(v)/4))
    #return trimBufferDict, statDict, optionDir
    return trimBufferDict, optionDir

##################
## Level Three ##
##################

########################
#### readFastq Functions

def mac_code_check(input_line, mac_code):
    if mac_code == 'auto':
        mac_code = string.split(input_line, sep=":")[0]
        sys.stdout.write("Machine Code : {0}\n".format(mac_code))
        if not(mac_code.startswith('@')):
            sys.stdout.write("Warning! Possible FastQFile corruption!\n")
            sys.stdout.write("Sequencer Machine Code : {0}".format(mac_code))
    return mac_code
    

def add_to_trimBufferDict(item, trimBufferDict, optionDir):
    """
    """
    trimmedItem = trimItem(item, optionDir)
    if trimmedItem[0] == None:
        bc_key = '_OUT_'
    else:
        bc_key = trimmedItem[0]
    barcode_list = trimBufferDict.get(bc_key, [])
    barcode_list.extend(trimmedItem[1])
    trimBufferDict[bc_key] = barcode_list
    return trimBufferDict


########function
########endblock


def trimItem(fastqItem, optionDir):
    """
    manages different trim functions for different modes
    """
    if optionDir['MODE'] == 'default':
        return stdTrimItem(fastqItem, optionDir)
    elif optionDir['MODE'] == 'crspr':
        return crsprTrimItem(fastqItem, optionDir)
    else:
        return stdTrimItem(fasqItem, optionDir)

        

def stdTrimItem(item, optionDir, trimNumber = 24, barcodeNumber = 8):
    """
    searches for Barcode and return Fastq item with sequence reverse
    complemented and trimmed directly to first base of barcode
    Only returns sequence if barcode is found in configuration file.
    Checks barcode and allows one occurence of N with function sameKey
    """
               #              2+Spacer:Barcode:+4
    emPatSet = ['ACAGCAATATACTG','ACTG[\w]{8}ATCT']

    if "adapter" in optionDir.keys():
        adapterSeq = string.strip(optionDir['adapter'])
        patSet = optionDir['precompiledSet']
    else:
        # adapterSeq = "AATATACTG"
        # mir30 + spacer TG
        # adapterSeq = "ATCTCGTATGCCGTCTTCTGCTTG" #reverse comp
        # careful adjust trimNumber accordingly
        
        # adapter P7
        #adapterSeq = 'ATCTCGTATGCCGTCTTCTGCTTG' #reverse comp
        adapterSeq =  'ATCTCGTATGCC'
        patSet = precompAdapter(adapterSeq, optionDir)

    retItem = list(item)
    seq = item[1]
    Barcode = None
    #for testing reduced patSet
    #pdb.set_trace()
    for pattern in patSet[len(patSet)-1:]:
        if pattern.search(seq):
            trimNumber = pattern.search(seq).start()
            retItem[1] = item[1][:trimNumber]
            retItem[3] = item[3][:trimNumber]
            Barcode = "Set"
            #sys.stdout.write("Pattern 1\n")
            #sys.stdout.write("Before Trimming:\n{0}\n".format(item))
            #sys.stdout.write("After Trimming:\n{0}\n".format(retItem))
            #sys.exit()
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
        # Barcode2 = reverseComplement(retItem[1][-8:])
        # Debugger
        # if (Barcode2 != Barcode
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
   

def crsprTrimItem(fastqItem, optionDir):
    """
    """
    ### init
    allowedMM = optionDir['ALLOWEDMM']
    distance_to_barcode = 83
    #distance_to_barcode = 24
    ## init fastqItem
    Barcode = None
    retItem = list(fastqItem)
    seq = fastqItem[1]
    ## init searchmask
    if "adapter" in optionDir.keys():
        adapterSeq = optionDir['adapter']
        patSet = optionDir['precompiledSet']
    else:
        adapterSeq  = "cgtcctttccacaagatatataaagccaag".upper()
        patSet = precompAdapter(adapterSeq, optionDir, clipFromTail=False)
    ## search mask
    #for pat in patSet[len(patSet)-1:]:
    #    if pat.search(seq):
    #        trimNumber = pat.search(seq).end()+24+8
    #    else:
    #        return None, fastqItem
    ## hamming hack
    foundAdapter = findPatInHammingDist(adapterSeq.upper(), seq, allowedMM)
    if foundAdapter:
        for position in foundAdapter:
            k = min(foundAdapter.keys())
            trimNumber = foundAdapter[k][2]+distance_to_barcode+8
    else:
        #print foundAdapter
        return None, fastqItem

    ## trimm for barcode
    if trimNumber > len(seq):
        return "TooShort", fastqItem
    retItem[1] = fastqItem[1][:trimNumber]
    retItem[3] = fastqItem[3][:trimNumber]
    ## turnaround
    retItem[1] = reverseComplement(retItem[1])
    retItem[3] = retItem[3][::-1]
    Barcode = retItem[1][:8]
    #print Barcode
    #print retItem
    #print "...."
    return Barcode, retItem

######Level three####################


def findPatInHammingDist(pat, seq, maxdist):
    #pdb.set_trace()
    patlen = len(pat)
    foundPat = {}
    for i in range(len(seq)-patlen):
        chunk = seq[i:i+patlen]
        distObj = SeqComp(pat, chunk)
        dist = distObj.get_ham_dist()
        if dist <= maxdist:
            foundPat[dist] = [dist, i, i+patlen]
    return foundPat

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
        
def precompAdapter(adapterSeq, optionDir, clipFromTail=True):
    """
    precompiles the adapter into 
    regex up to the minimum length of adapterMinLen
    returns a list of the re objects
    """
    sys.stdout.write("calculating ....{0}\n".format(adapterSeq))
    retList = [re.compile(adapterSeq)]
    seq = adapterSeq
    while(len(seq) > optionDir['adapterMinLen']):
        if clipFromTail:
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
    seqMachineId is     @M01950 miseq and
                        @ST-K00207 hiseq

                'adapter': ''#'ATCTCGTATGCCGTCTTCTGCTTG',\
    """
    ValueDir = {\
                'oFile': ['temp','txt'],\
                'verbose': True,\
                'adapter': "CGTCCTTTCCACAAGATATATAAAGCCAAG",\
                'adapterMinLen' : 30,\
                'seqMachineId' : "auto",\
                'MODE' : 'default',\
                'ALLOWEDMM' : 2,\
                'DISTANCEBARCODE' : 24,\
                'maxChunkSize' : 200000000}#1073741824} # in Bytes; 1073741824 bytes = 1 Gigybyte
    return ValueDir

#######Main#######

def main():
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, b:, a:, s, d:, t, c:, m:",\
                ["help", "output=", "barcodes=", "adapter=", "silent", "destination=",\
                "automatic", "configFile", "mode=" ])
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
        if o in ['-m', '--mode']:
            optionDir['MODE'] = string.strip(a)
    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()


