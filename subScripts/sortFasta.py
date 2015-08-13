#!/usr/bin/env python

import sys, string, os, getopt
import datetime

"""
divides the source fasta file
into the different pools
saves every pool in a different fasta file.

"""
######Level One##########

def manage(args, optionDir):
    """
    """
    for a in args:
        print a
        sortingDict = sortA(a)
        writeSoDi(sortingDict, optionDir)

def processArgs(args, optionDir):
    """
    """
    return args, optionDir


def sortA(a):
    """
    """
    switch = False
    retDict = {}
    header = ""
    Poolheader = ""
    seq = ""
    ihandle = open(a)
    for line in ihandle:
        mLine = string.strip(line)
        if mLine.startswith(">"):
            if Poolheader in retDict.keys():
                retDict[Poolheader].append([header, seq])
            else:
                retDict[Poolheader] = [[header,seq]]
            Poolheader = string.split(mLine, sep=":")[1]
            header = mLine
            seq = ""

        else:
            seq += mLine

    if Poolheader in retDict.keys():
        retDict[Poolheader].append([header, seq])
    else:
        retDict[Poolheader] = [[header,seq]]
    ihandle.close()
    return retDict

def writeSoDi(sortingDict, optionDir):
    """
    """
    for k, v in sortingDict.items():
        if 'destination' in optionDir.keys():
            ohandle = open(optionDir['destination']+"/"+k+"_pool.fasta", 'w')
        else:    
            ohandle = open(k+"_pool.fasta", "w")
        for entry in v:
            ohandle.write(entry[0]+"\n"+entry[1]+"\n")
        ohandle.close()

######Defaults###########

def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['temp','txt']\
               }
    return ValueDir

#######Main#######

def main():
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, d:",\
                               ["help", "output=", "destination="])
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
            optionDir['destination'] = os.path.abspath(string.strip(a))

    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()




