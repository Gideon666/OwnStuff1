#!/usr/bin/env python

"""
statmaker called from SeqAnalyzer.sh

usage:
    ./statMaker.py [options] -s statFile -e _pStat.b -c barcodeFile
"""

import shelve
import sys, string, os, getopt
import datetime
import re
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import plotfile, show, gca
import pickle
import operator
from os import walk
import math
import pdb

######Level One##########

def manage(args, optionDir):
    """
    """
    #pdb.set_trace()
    calcStats(args, optionDir)
    clean_up(args, optionDir)

def processArgs(args, optionDir):
    """
    """
    sys.stdout.write("Processing Config....\n")
    if 'configFile' in optionDir.keys():
        optionDir['config'] = readConfig(optionDir)
    optionDir['statData'] = {}
    sys.stdout.write("Processing Stats....\n")
    if 'statFile' in optionDir.keys():
        statDict, totalDict = processStatFile(optionDir['statFile'])
        optionDir['NumberDict'] = statDict
        optionDir['TotalsDict'] = totalDict
    return args, optionDir

######Level Two######

def calcStats(args, optionDir):
    """
    default stats calculater
    """
    sys.stdout.write("Calculating Stats....\n")
    writeCSVStats(args, optionDir)
    sys.stdout.write('CSVStats written!\n')
    args, optionDir = getData(args, optionDir)
    sys.stdout.write('read Blat data!\n')
    makePlots(args, optionDir)
    sys.stdout.write('Plots written!\n')
    writeCSVhairpins(args, optionDir)
    sys.stdout.write('CSVhairpins written!\n')

def clean_up(args, optionDir):
    #optionDir['firstBlatData_{0}'.format(optionDir['check_for_Ns'])].close()
    #os.remove(optionDir['destination']+"statMaker_shelve.tmp")
    pass

def processStatFile(eStatfile):
    """
    read in File and create the correct 
    structure
    """
    statDict = readInStatFile(eStatfile)
    statStruct = structStats(statDict)
    return statStruct[0], statStruct[1]


def structStats(statDict):
    """
    structures the statistic Data
    """
    #StatDict[Fieldcode] = [Barcode, Infos, ...]
    oneFieldKeys = ["BR", "CR", "FR", "HR", "SR", "LR"]
    twoFieldKeys = ["BC", "CU", "CO",       "SC", "LN", "MB", "VB"]
    threeFieldKeys = ["HF"]
    BarcodeDict = {}
    TotalInfo = {}
    turn2Number = lambda x : string.atoi(x)
    #### Twofields
    for info in twoFieldKeys:
        try:
            for BarcodeList in statDict[info]:
                try:
                    if BarcodeList[0] in BarcodeDict.keys():
                        BarcodeDict[BarcodeList[0]][info] = string.atoi(BarcodeList[1])
                    else:
                        BarcodeDict[BarcodeList[0]] = {}
                        BarcodeDict[BarcodeList[0]][info] = string.atoi(BarcodeList[1])
                except KeyError:
                    sys.stderr.write("KeyError: {0}\n".format(BarcodeList))
                    sys.stderr.write("BarcodeDict: {0}\n".format(BarcodeDict.keys()))
                    if BarcodeList[0] in BarcodeDict.keys():
                        BarcodeDict[BarcodeList[0]][info] = 0
                    else:
                        BarcodeDict[BarcodeList[0]] = {}
                        BarcodeDict[BarcodeList[0]][info] = 0
        except KeyError:
            sys.stderr.write("KeyError: {0}\n".format(statDict))
            sys.exit(-1)
    #### Onefields
    for tInfo in oneFieldKeys:
        try:
            TotalInfo[tInfo] = string.atoi(statDict[tInfo][0][0])
        except KeyError:
            sys.stdout.write("Key Error : {0} not found!\n".format(str(tInfo)))

    #### Threefields
    for multipleInfo in threeFieldKeys:
        for BarcodeList in statDict[multipleInfo]:
            if BarcodeList[0] in BarcodeDict.keys():
                BarcodeDict[BarcodeList[0]][multipleInfo] = map(turn2Number, BarcodeList[1:])
            else:
                BarcodeDict[BarcodeList[0]] = {} 
                BarcodeDict[BarcodeList[0]][multipleInfo] = map(turn2Number, BarcodeList[1:])
    #### compatibility reasons
    
    BarcodeDict = fillupInfos(BarcodeDict, twoFieldKeys)
    BarcodeDict = fillupInfos(BarcodeDict, threeFieldKeys, keyLen=2)
    
    #print "------------"
    #print "Stat Info"
    #print str(BarcodeDict)
    #print "Totals"
    #print str(TotalInfo)
    #print "------------"
    return BarcodeDict, TotalInfo


def readInStatFile(eStatfile):
    """
    read in external data from stats.txt file.
    
    Two char code:
    BC : Barcode Split : 2
    BR : Barcode Total Reads (high number) :1
    CU : Cutout : 2
    CR : Cutout Total Reads : 1
    CO : collapsed Sequences to fasta :2
    FR : Fasta Total Reads :1
    HF : Number of BLAT Hits in Sequences : BLAT Hits Filtered : 3
    HR : Hits Total : 1
    SC : Second collapse unique BLAT hits : 2
    SR : Second collapse total Hits : 1
    LN : Second collapse BLAT Reads : 2
    LR : Second collapse total Reads : 1
    MB : mask(loop) in barcode split files
    VB : mirsequence in barcode split files
    
    e.g
    BC  ACGTACGT        10000
    BR  100000
    Should look like this
    """
    retDict = {}
    ohandle = open(eStatfile, 'r')
    #linecrawling
    for line in ohandle:
        if string.strip(line) == '':
            continue
        if line.startswith('#') or line.startswith('@'):
            continue
        atoms = string.split(string.strip(line))
        header = atoms[0]
        info = atoms[1:]
        if len(info) > 1:
            if atoms[1] not in optionDir['config'].keys():
                continue
        if header in retDict.keys():
            retDict[header].append(info)
        else:
            retDict[header] = [info]
    ohandle.close()
    return retDict

######Level Two######

def writeCSVStats(args, optionDir):
    """
    """
    DataDict = optionDir['NumberDict']
    tDict = optionDir['TotalsDict']
    outputString = ""
    header = "#SampleID\tDescription\tBarcode\tBarcoded Reads\t\
PerfectAnchor\t5'mir Sequence\tCutout Sequences\tCollapsed to # unique Sequences \tshRNAs Blat Hits\t\
Blat Hits with total # of Reads\tHit Rate\t\tTotalBarcodeReads\n"

    for entry, values in DataDict.items():
        if not(entry in optionDir['config'].keys()):
            continue
        keys = values.keys()
        sys.stdout.write("CSVStats : "+str(entry)+" "+str(values)+"\n")
        outputString += optionDir['config'][entry][0]+"\t"+optionDir['config'][entry][3]+"\t"
        outputString += entry+"\t"
        outputString += 'BC' in keys and str(values['BC'])+"\t" or "0\t"
        outputString += 'MB' in keys and str(values['MB'])+"\t" or "0\t"
        outputString += 'VB' in keys and str(values['VB'])+"\t" or "0\t"
        outputString += 'CU' in keys and str(values['CU'])+"\t" or "0\t"
        outputString += 'CO' in keys and str(values['CO'])+"\t" or "0\t"
        outputString += 'SC' in keys and str(values['SC'])+"\t" or "0\t"
        outputString += 'LN' in keys and str(values['LN'])+"\t" or "0\t"
        try:
            outputString += str((100.00 * values['LN'])/values['BC'])+"\t\t"+str(tDict['BR'])+"\n"
        except ZeroDivisionError:
            outputString += "0%\t\t0\n"

    with open(optionDir['destination']+"/"+optionDir['oFile'][0]+"."+optionDir['oFile'][1], 'w') as ohandle:
        ohandle.write(header+outputString)
    ohandle.close()

def writeCSVhairpins(args, optionDir):
    """
    """
    if optionDir['verbose'] : sys.stdout.write("Writing Pandafile...\n")
    data = optionDir['blatResults']
    columnum = len(data.keys())
    rownum = max(map(len, data.values()))
    np.zeros((rownum, columnum))
    readDict = {}
    bcOrder = sorted(data.keys())
    row = []
    for bc in bcOrder:
        indexName = []
        column = []
        for hpId, hpData in data[bc].items():
            namesplit = string.split(hpId, sep="-")
            if len(namesplit) > 1:
            #    hpName = namesplit[0]+"-"+hpId[-7:]
                 hpName = "_".join(namesplit[:-1])+"-"+namesplit[-1]
            else:
                hpName = namesplit[0]
            reads = string.atof(string.split(hpData[9], sep="-")[1])
            #print reads
            indexName.append(hpName)
            column.append(reads)
        #readDict[bc] = pd.Series(column, index=indexName)
        try:
            readDict[bc+"_("+optionDir['config'][bc][3]+")"] = pd.Series(column, index=indexName)
        except IndexError:
            sys.stdout.write("Index Error in writeCSVhairpins statMaker.py ... {0}\n".format(bc))
            readDict[bc] = pd.Series(column, index=indexName)
    #pdb.set_trace()
    pDF = pd.DataFrame(readDict)
    pDF.to_csv(optionDir['destination']+"/"+"pandatest.tsv", sep="\t")
    #save for future work to h5
    store = pd.HDFStore(optionDir['destination']+"/"+'statStore.h5')
    store['pDF'] = pDF
    #print pDF

def getData(args, optionDir, pat=None, LC="BL"):
    """
    read first blat results
    """
    N_check = optionDir['check_for_Ns']
    bl = {}
    #bl = shelve.open(optionDir['destination']+"statMaker_shelve.tmp",writeback=True)
    list = []
    if not(pat):
        pat = '(\w{5,8})_.*%s.pslx' % LC
    for (dirpath, dirnames, filenames) in walk(optionDir['destination']):
        list.extend(filenames)
        break
    if optionDir['verbose']:
        #sys.stdout.write("Loading list for getData:\n{0}\n".format(str(list)))
        sys.stdout.write("config list:\t{0}\n".format(optionDir['config'].keys()))
    for fName in list:
        barc = re.search(pat, fName)
        if barc:
            if fName[:8] not in optionDir['config'].keys():
                continue
            if optionDir['verbose'] : sys.stdout.write("Loading : %s\n"% (dirpath+fName))
            #bl[barc.group(1)] = loadIn(dirpath+fName, args, optionDir)
            blat_data_this = loadIn(dirpath+fName, args, optionDir)
            bl.update(getNumbersFromBlat(barc.group(1), blat_data_this, N=N_check, M=0))
    optionDir['firstBlatData_{0}'.format(str(N_check))] = bl
    return args, optionDir

def makePlots(args, optionDir):
    """
    """
    #### init
    DataDict = optionDir['NumberDict']
    TotalDict= optionDir['TotalsDict']
    blatResDict = {}
    #blatResDict2 = getNumbersFromBlat(args, optionDir, N=10, M=0)
    blatResDict2 = optionDir['firstBlatData_{0}'.format(str(optionDir['check_for_Ns']))]
    allowed_Ns = optionDir['allowed_Ns_in_hit'] 
    for k, v in DataDict.items():
        pD = getPD(k, args, optionDir)
        if pD:
            blatResDict[k] = pD
        else:
            sys.stderr.write("Deleted {0} entry from Blastresults, no hits...\n".format(k))
    #### save blatResDict 
    optionDir['blatResults'] = blatResDict 
    #### plot init
    #### Pie and Bar Plots ######################################################################
    ####
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(30, 15))
    ax0, ax1, ax2, ax3 = axes.flat
    #### plotting
    
    makePie(ax0, DataDict, TotalDict, args, optionDir, title="Barcode Splitting")
    makeBoxPlot(ax1, blatResDict, args, optionDir, title="Sequence Read Quantity N<={0}".format(str(allowed_Ns)), log=True)
    makePie(ax2, DataDict, TotalDict, args, optionDir, title="Perfect Hits", LC=("LN","LR"))
    makeBoxPlot(ax3, blatResDict2, args, optionDir, title="Sequcence Read Quantity N<=10", log=True)

    #### generell things
    plt.tight_layout()
    plt.savefig(optionDir['destination']+"/"+optionDir['oFile'][0]+"_OV.svg")
    
    #### plot init 
    #### Dot Plots ##############################################################################
    ####
    groups = groupPlot(optionDir)
    plotDict = prepareShReadData(blatResDict, optionDir)
    sumFG = lambda x : (x * (x +1)) / 2
    #pdb.set_trace()
    for ind in xrange(len(groups)):
        group= groups[ind]
        grpSize = len(group[2])
        plotNum = sumFG(grpSize-1)
        #print "ind: {0}\ngrpSize: {1}\nplotNum: {2}\n".format(ind, grpSize, plotNum)
        #print groups[ind]
        fig, axes = plt.subplots(nrows=2, ncols=plotNum, figsize=(min(plotNum, 16)*10, 20), dpi=80)
        axesList = axes.flat
    #### plotting
        if (grpSize <=1):
            continue
        pairList = onePlotRow(group[2])
        print pairList, group[2]
        for i in xrange(0, plotNum):
            #split = group[2][i:i+2] 
            split = pairList[i]
            plotDotPlot(axesList[i], plotDict, split, args, optionDir,\
                    mode=['dynLim' ,'color', 'default', "annotate"])
            plotDotPlot(axesList[i+plotNum], plotDict, split, args, optionDir,\
                    mode=['dynLim',"shOnly", "annotate"]) 
            #plotDotPlot(axesList[i+grpSize-1], plotDict, split, args, optionDir, mode=["shOnly"]) 
    #### generell things
        plt.tight_layout()
        plt.savefig(optionDir['destination']+"/"+optionDir['oFile'][0]+"_grp"+str(ind)+"_dotPlot.svg")
    print groups

def makePie(ax, DataDict, TotalDict, args, optionDir, title="", LC=('BC', 'BR')):
    """
    """
    sizes = []
    labels = []
    for entry, values in DataDict.items():
        try:
            sizes.append(100 * values[LC[0]]/float(TotalDict[LC[1]]))
        except (KeyError, ZeroDivisionError) as e:
            sizes.append(0.00001)
            sys.stdout.write("Warning Error in makePie Values too small!\n{0}\n{1}\n"
                    .format(entry, e))
        labels.append(entry)
    #for i in xrange(len(sizes)):
    #    labels[i] = labels[i]+"\n"+str(sizes[i])+"%"
    ax.set_title(title, fontsize=10)
    patches, texts, autotexts = ax.pie(sizes, labels=labels, shadow=False, \
            autopct = '%.2f' , labeldistance=0.750, colors=optionDir['colorSet'][:len(sizes)])
    for t in texts:
        t.set_fontsize(8)
    ax.axis('equal')
    #plt.show()

def makeBoxPlot(ax, blatResDict, args, optionDir, title="", log=False):
    """
    receives the information from the result dicts...
    
    """
    getInt = lambda x : string.atoi(string.split(x[9], sep="-")[1])
    labels = []
    data = []
    for k in sorted(blatResDict.keys()):
        v = blatResDict[k]
        #for k, v in blatResDict.items():
        reads = map(getInt, v.values())
        data.append(reads)
        #### Load Labels from Barcode File ####
        ##
        if 'config' in optionDir.keys():
            labels.append(optionDir['config'][k][3])
        else:
            labels.append(k)
    # extends labels to the same length
    labels = adjust_labels(labels)
    ax.set_title(title, fontsize=12)
    boxDict = ax.boxplot(data, notch=0, vert=1, whis=1.5, sym='+')
    plt.setp(boxDict['boxes'], color='black')
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_ylabel('Number of Reads')
    if log:
        ax.set_yscale('log')
    ax.set_xlabel('Pools')
    xtickNames = plt.setp(ax, xticklabels=labels)
    plt.setp(xtickNames, rotation=45, fontsize=8)

    #ax.set_ylim(0,6000)

def plotDotPlot(ax, DataDict, BCkeys, args, optionDir, title="", mode=["default"]):
    """
    """
    xlabel = ""
    ylabel = ""
    softLim = [0.95, 0.95]
    annoSoftLim = [0.5, 0.5]
    limits = [2000, 2000]
    annoLim = [1800, 1800]
    getShName = lambda x : string.split(x, sep="-")[0]
    try:
        xDataDict = DataDict[BCkeys[0]]
        yDataDict = DataDict[BCkeys[1]]
    except KeyError:
        sys.stderr.write("Warning in DotPlot! Dataset : {0}-{1} not available!\n".format(BCkeys[0], BCkeys[1]))
        return
    if 'config' in optionDir.keys():
        xlabel = optionDir['config'][BCkeys[0]][3]
        ylabel = optionDir['config'][BCkeys[1]][3]
        title = xlabel+" -vs- "+ylabel
    labels = []
    plotX = []
    plotY = []
    overlapC = 0
    maxX = 10
    maxY = 10

    if "shOnly" in mode:
        xDataDict = compressiTags(xDataDict)
        yDataDict = compressiTags(yDataDict)
    
    if 'color' in mode:
        #print "ColorMode"
        shName = []
        maxXLi = []
        maxYLi = []
        for k, v in xDataDict.items():
            if k in yDataDict.keys():
                sN = getShName(k)
                if sN in shName:
                    ind = shName.index(sN)
                else:
                    shName.append(sN)
                    ind = shName.index(sN)
                    plotX.append([])
                    plotY.append([])
                    labels.append([])
                try:
                    maxXLi.append(v)
                    plotX[ind].append(v)
                    maxYLi.append(yDataDict[k])
                    plotY[ind].append(yDataDict[k])
                    labels[ind].append(k)
                    overlapC += 1
                except IndexError:
                    print ind, shName, sN, len(plotX), plotX, plotY
                    print "IndexError"
                    exit()
        maxXList = sorted(maxXLi)
        maxYList = sorted(maxYLi)
    else:
        for k, v in xDataDict.items():
            if k in yDataDict.keys():
                plotX.append(v)
                plotY.append(yDataDict[k])
                labels.append(k)
                overlapC +=1
        maxXList = sorted(plotX)
        maxYList = sorted(plotY)
    # dynamic limits
    if (("dynLim" in mode) and (len(maxXList)>0)):

        maxX = max(maxXList)
        maxY = max(maxYList)
        XLen = len(maxXList)
        YLen = len(maxYList)
        if ((maxX < 100) and (maxY < 100)):
            limits = [100, 100]
            annoLim = [0,0]
        else:
            limits[0] = maxXList[int(XLen * softLim[0])]
            limits[1] = maxYList[int(YLen * softLim[1])]
            annoLim = map(lambda x, y : x*y, limits, annoSoftLim)
    ##
    title += " ({0},{1})-{2}".format(len(xDataDict.keys()), len(yDataDict.keys()), overlapC)
    ax.set_title(title, fontsize=14)
    if 'color' in mode:
        for ind in xrange(len(plotX)):
            ax.plot(plotX[ind], plotY[ind], ".")
            if 'annotate' in mode:
                for pointInd in xrange(len(plotX[ind])):
                    if ((plotX[ind][pointInd] > annoLim[0]) or \
                            (plotY[ind][pointInd] > annoLim[1])):
                        ax.annotate('{0}'.format(labels[ind][pointInd]),\
                                xy=(plotX[ind][pointInd], plotY[ind][pointInd]), xytext=(10,0),\
                                textcoords='offset points')
                        
    else:
        ax.plot(plotX, plotY, 'r.')
        if 'annotate' in mode:
            for pointInd in xrange(len(plotX)):
                ax.annotate('{0}'.format(labels[pointInd]), xy=(plotX[pointInd], plotY[pointInd]),\
                        xytext=(10,0), textcoords='offset points')
    if "default" in mode:
        xlimDef = ax.get_xlim()
        ylimDef = ax.get_ylim()
    #dynamicsetting
    ax.set_xlim([0,limits[0]])
    ax.set_ylim([0,limits[1]])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #plt.setp(ax, xticklabels=labels)
        #print k
        #print len(v.values())
        #print len(labels)
        #break

######Level Three######

def adjust_labels(label_list):
    ret_label = ""
    ret_list = []
    max_len = max(map(len, label_list))
    f = lambda s: max_len - len(s)
    for label in label_list:
        ret_label = label+"_"*f(label)
        ret_list.append(ret_label)
    return ret_list



def fillupInfos(barcodeDict, fieldKeys, keyLen=1):
    """
    """
    retDict = {}
    for k, v in barcodeDict.items():
        retDict[k] = v
        for fK in fieldKeys:
            if not(fK in v.keys()):
                if keyLen == 1:
                    retDict[k][fK] = 0
                else:
                    retDict[k][fK] = [0]*keyLen
    return retDict

def onePlotRow(grpList):
    """
    """
    newList = grpList[:]
    if len(newList) == 2:
        return [newList]

    retList = []
    element = newList.pop()
    for i in newList:
        retList.append([element, i])
    retList.extend(onePlotRow(newList))
    return retList

def compressiTags(itagDict):
    """
    """
    getName = lambda x : string.split(x, sep="-")[0]
    retDict = {}
    for k, v in itagDict.items():
        nKey = getName(k)
        if nKey in retDict.keys():
            retDict[nKey] += v
        else:
            retDict[nKey] = v
    return retDict


def groupPlot(optionDir):
    """
    """
    groups = []
    rRows = 1
    rCols = 1
    BCkeys = []
    groupDir = {}
    counter = 0
    if 'config' in optionDir.keys():
        cfg = optionDir['config']
        if len(cfg.values()[0]) >= 5:
            for k, v in cfg.items():
                gId = v[4]
                if gId == "x":
                    sys.stderr.write("Excluded : {0}\n".format(k))
                    continue
                if gId in groupDir.keys():
                    groupDir[gId].append(v[1])
                else:
                    groupDir[gId] = [v[1]]
    for k, v in groupDir.items():
        plotNum = len(v)
        rCols = (plotNum / 2) + 1
        rRows = (plotNum / rCols) + 1
        BCkeys = v[:]
        groups.append((rRows, rCols, BCkeys))
    #print "All Groups:"
    #print groups
    return groups


def prepareShReadData(blatResDict, optionDir):
    """
    """

    getReadNumber = lambda x : string.atoi(string.split(x[9], sep="-")[1])
    getiTagName = lambda x : string.split(x[13], sep="-")[0]+"-"+string.split(x[9], sep="-")[2]
    getName = lambda x : string.split(x[13], sep="-")[0]
    retDict = {}
    for k, v in blatResDict.items():
        pool = {}
        if optionDir['iTags']:
            for entry in v.values():
                try:
                    number = getReadNumber(entry)
                    name = getiTagName(entry)
                    pool[name] = number
                except IndexError:
                    sys.stderr.write("IndexError in prepareShReadData\
                            for {1} : {0}\n".format(entry))
        else:
            for entry in v.values():
                number = getReadNumber(entry)
                name = getName(entry)
                pool[name] = number
        
        retDict[k] = pool
    return retDict


def getNumbersFromBlat(barcode, hitList, N=0, M=0):
    """
    filters data from the blat result file pslx.
    in this case the first blat file to make a comparison
    possible.
    """
    retDict = {}
    #sDict = optionDir['firstBlatData']
    retDict[barcode] = {}
    for hit in hitList:
    #for barcode, hitList in sDict.items():
        #print barcode, hitList
        #for hit in hitList:
            match = string.atoi(hit[0])
            mm = string.atoi(hit[1])
            ns = string.atoi(hit[3])
            size = string.atoi(hit[10])
            if (((mm <= M)\
                    and (ns <= N))\
                    and (match+mm+ns==size)):
                #print hit
                try:
                    info = string.split(hit[9], sep="-")
                except IndexError:
                    print hit
                    sys.exit(-1)
                hitname = hit[13]
                number = string.atoi(string.strip(info[1]))
                if (hitname in retDict[barcode].keys()):
                    #print retDict[barcode][hitname]
                    oldnum = string.atoi(string.split(retDict[barcode][hitname][9], sep="-")[1])
                    retDict[barcode][hitname][9] = "%s-%s"% (info[0], oldnum+number)
                    #print hitname, hit
                else:
                    retDict[barcode][hitname] = hit
    return retDict                

def loadIn(filename, args, optionDir):
    """
    """
    retList = []
    startData = False
    with open(filename, 'r') as ihandle:
        for line in ihandle:
            if not(startData):
                if  line.startswith("-"):
                    startData = True
                continue
            atoms = string.split(line)
            retList.append(atoms)
    return retList

def getPD(name, args, optionDir):
    """
    get Pickle Data
    """
    nameSet = name
    if 'exFile' in optionDir.keys():
        try:
            retDict = pickle.load(open(optionDir['destination']+nameSet+optionDir['exFile'], 'rb'))
        except IOError:
            sys.stderr.write("External File for {0} not found!\n".format(nameSet))
            return None
    #print retDict.keys(), retDict.values()
    return retDict
    
def readConfig(optionDir):
    """
    """
    barcodes = {}
    with open(optionDir['configFile'], 'r') as ihandle:
        for line in ihandle:
            if string.strip(line) == "":
                continue
            if not(line.startswith("#") or line.startswith("@")):
                atoms = string.split(string.strip(line), sep="\t")
                barcodes[atoms[1]] = atoms
    ihandle.close()
    return barcodes


######Defaults###########

def setDefaultValues():
    """
    standard Values settings
    """
    ValueDir = {\
                'oFile': ['Stats','tsv'],\
                'destination' : "./",\
                'iTags': False,\
                'verbose': True,\
                'check_for_Ns': 10,\
                'colorSet' : [\
                'yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'darksalmon',\
                'orangered', 'palegoldenrod','azure', 'darkorchid', 'firebrick',\
                'springgreen', 'lightgrey','steelblue', 'sandybrown','crimson',\
                'forestgreen', 'lightyellow','darkslateblue', 'moccasin', 'fuchsia',\
                'lawngreen', 'honeydew', 'lavender', 'oldlace', 'palevioletred',\
                'olive','lemonchiffon', 'indigo', 'rosybrown','indianred']
                }
    return ValueDir

#######Main#######

def main():
    global optionDir
    optionDir = setDefaultValues()
    try:
        opts, args = getopt.getopt(sys.argv[1:],\
                "h, o:, s:, d:, e:, c:, i, v, n:",\
                ["help", "output=", "statFile=", "destination=", "externalD=",\
                "configFile=", "iTags", "allowed_Ns_in_hit", "verbose"])
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
        if o in ['-s', '--statFile']:
            optionDir['statFile'] = string.strip(a)
        if o in ['-e', '--externalD']:
            optionDir['exFile'] = a
        if o in ['-d', '--destination']:
            optionDir['destination'] = a
        if o in ['-c', '--configFile']:
            optionDir['configFile'] = string.strip(a)
        if o in ['-i', '--iTags']:
            optionDir['iTags'] = True
        if o in ['-v', '--verbose']:
            optionDir['verbose'] = True
        if o in ['-n', '--allowed_Ns_in_hit']:
            optionDir['allowed_Ns_in_hit'] = int(string.strip(a))

    args, optionDir = processArgs(args, optionDir)
    manage(args, optionDir)

if __name__ == "__main__":
    main()
