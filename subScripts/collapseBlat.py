#!/usr/bin/env python

"""
collapses the pslx blat result
usage: collapseBlat.py [options} "Blat Files"

-h, --help      : displays this help doc
-o, --output    : output file
"""

import sys, string, os, getopt
import re
import numpy as np
import datetime
import matplotlib.pyplot as plt
from operator import itemgetter
from pylab import plotfile, show, gca
#try:
    #    import cPickle as pickle
#except:
import pickle

######Level One##########

def manage(args, optionDir):
    """
    """
    for a in args:
        retDict, optionDir = collapse(a, optionDir)
        if 'doStats' in optionDir.keys():
            optionDir = calcStats(a, retDict, optionDir)
            barcodeO = re.search("/(.{5,8})_Pool", a)
            name = barcodeO.group(1)
            pickle.dump(retDict, open(optionDir['destination']+"/"+name+"_pStat.b", "wb"))
        saveDict(a, retDict, optionDir)

def processArgs(args, optionDir):
    """
    """
    if 'configFile' in optionDir.keys():
        optionDir['config'] = readConfig(optionDir['configFile'])
    optionDir['statData'] = {}
    for a in args:
        optionDir['statData'][a] = []
    #if 'statsFile' in optionDir.keys():
    #    statDict = processStatFile(optionDir['statsFile'])
    #    optionDir['externalStatData'] = statDict
    return args, optionDir


######Level Two######


def processStatFile(eStatfile):
    """
    #read in File and create the correct 
    #structure
    """
    statDict = readInStatFile(eStatfile)
    statStruct = structStats(statDict)

def readInStatFile(file):
    pass

def structStats(dict):
    pass


def calcStats(filename, retDict, optionDir):
    """
    """
    readDict = optionDir['statData'][filename][0]
    sortedReads = sorted(readDict.items(), key=itemgetter(1))
    minReads = sortedReads[:10]
    maxReads = sortedReads[::-1][:10]
    add = lambda x, y : x+y
    nSplit = lambda x : string.split(x, sep="-")[0]
    try:
        readCounter = reduce(add, readDict.values())
    except TypeError:
        readCounter = 0
        sys.stderr.write("TypeError in %s \n"% filename)
        return optionDir
    index = np.arange(len(readDict.values()))
    print filename, readCounter 
    val = readDict.values()
    #val.append(15000)
    values = sorted(val, reverse=True)
    #print values
    ##########
    #drawing!!!
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(30,15))
    # overviewplots
    ax0, ax1, ax2, ax3 = axes.flat
    ax0.hist(readDict.values(),bins=100, color='springgreen')        #len(readDict.values()),range=(0,25))
    ax2.bar(index, values, color='steelblue')                      #height=10000)
    if 'config' in optionDir.keys():
        ax0.set_title(optionDir['config'][os.path.basename(filename)[:8]][3]+" Histogram")
        ax2.set_title(optionDir['config'][os.path.basename(filename)[:8]][3]+" BarPlot")
    else:
        ax0.set_title(os.path.basename(filename)+"_hist")
        ax2.set_title(os.path.basename(filename)+"_bar")

    # highes lowest
    width = 0.8 
    index2 = np.arange(len(maxReads))
    index3 = np.arange(len(minReads))
    ##
    ax1.set_title("10 Highest")
    ax1.bar(index2, map(itemgetter(1), maxReads), color="crimson")
    ax1.set_xticks(index2+width/2)
    ax1.set_xticklabels(map(nSplit, map(itemgetter(0), maxReads)))
    ##
    ax3.set_title("10 Lowest")
    ax3.bar(index3, map(itemgetter(1), minReads[::-1]), color="lightgreen")
    ax3.set_xticks(index2+width/2)
    ax3.set_xticklabels(map(nSplit, map(itemgetter(0), minReads[::-1])))
    ###
    plt.tight_layout()
    plt.savefig(optionDir['destination']+"/"+os.path.basename(filename)+"_stat.png")
    plt.close()
    return optionDir

def saveDict(orgFilename, savedDict, optionDir):
    """
    """
    ofilename=optionDir['destination']+"/"+re.split("_HF", os.path.basename(orgFilename))[0]+"_SC.pslx"
    ohandle = open(ofilename, 'w')
    for k, v in savedDict.items():
        ohandle.write(v[0])
        for item in v[1:]:
            ohandle.write("\t")
            ohandle.write(item)
        ohandle.write("\n")
    ohandle.close()

def collapse(ifile, optionDir):
    """
    """
    retDict = {}
    readDict = {}
    ihandle = open(ifile, 'r')
    for line in ihandle:
        atoms = string.split(string.strip(line))
        id = atoms[13]
        number = string.split(atoms[9], sep="-")
        if (len(number) > 2):
            id += number[2]
        if id in retDict.keys():
            oldNumber = string.split(retDict[id][9], sep="-")
            newNumber = string.atoi(oldNumber[1])+string.atoi(number[1])
            if(len(oldNumber) > 2):
                retDict[id][9] = oldNumber[0]+"-"+str(newNumber)+"-"+str(oldNumber[2])
            else:
                retDict[id][9] = oldNumber[0]+"-"+str(newNumber)
            readDict[id] += string.atoi(number[1])
                
        else:
            retDict[id] = atoms
            readDict[id] = string.atoi(number[1])
    optionDir['statData'][ifile].append(readDict)
    ihandle.close()
    return retDict, optionDir


def readConfig(file):
    """
    """
    #sys.stdout.write('Reading configuration File {}\n'.format(file))
    retDict = {}
    with open(file, 'r') as ihandle:
        for line in ihandle:
            if (line.startswith('#') or line.startswith("@")):
                continue
            atoms = string.split(string.strip(line), sep="\t")
            retDict[atoms[1]] = atoms
    return retDict

######Defaults###########

def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['col','txt'],\
                'statOut' : ['plot', 'txt'],\
                'destination' : "./"
               }
    return ValueDir

#######Main#######

def main():
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, s, d:, f:, c:",\
                ["help", "output=", "statistics", "destination=", "statsFile=", "configFile="])
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
        if o in ['-d', '--destination']:
            optionDir['destination'] = os.path.dirname(string.strip(a))
        if o in ['-s', '--statistics']:
            optionDir['doStats'] = "Yes"
        if o in ['-c', '--configFile']:
            optionDir['configFile'] = string.strip(a)
    #    if o in ['-f', '--statsFile']:
    #        optionDir['statsFile'] = string.strip(a)

    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()
