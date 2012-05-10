#!/bin/env python
#
# Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)
#
# Released under the MIT license, see LICENSE.txt
#

import os
import re
import sys
import peakparser
import discordant
import argparse
import annotate

class insertionSummary:
    def __init__(self):
        self.chr      = ''
        self.start    = '' 
        self.end      = ''
        self.strand   = ''
        self.eltStart = ''
        self.eltEnd   = ''
        self.eltFam   = ''
        self.eltClass = ''
        self.index    = ''
        self.numreads = ''
        self.sources  = [] 
        self.annot    = []

    def __str__(self):
        out = ("\t".join((self.chr,self.start,self.end,self.strand,
                         self.eltStart,self.eltEnd,self.eltFam,
                         self.eltClass,self.index,self.numreads))
              + "\t" + ",".join(self.sources)
              + "\t" + "\t".join(self.annot))
        return out

    # sFS = strand from strand, sFR = strand from read positions
    def setBestStrand(self,sFS,sFR):
        if sFS == sFR and sFS != 'x':
            self.strand = sFS
            return 1
        if (sFS == 'i' or sFS == 'x') and sFR != 'x':
            self.strand = sFR
            return 1
        if sFR == 'x' and sFS != 'x' and sFS != 'i':
            self.strand = sFS
            return 1
        if sFR != 'x' and sFS !='x' and sFS != 'i' and sFR != sFS:
            self.strand = sFS
            return 1
        self.strand = '.'
        return 0

    # sets chr,start,end from chrX:NNNN-NNNN 
    def setPos(self,posString):
        if re.search(':',posString) and re.search('-',posString):
            (chrstart,self.end)   = posString.rsplit('-')
            (self.chr,self.start) = chrstart.rsplit(':')
            return 1
        else:
            return 0

    def setEltPos(self,posString):
        if re.search(',',posString):
            (self.eltStart,self.eltEnd) = posString.rsplit(',')
            return 1
        else:
            return 0

    def setSources(self,sourceString):
        if sourceString:
            for source in sourceString.rsplit(','):
                self.sources.append(source)
            return 1
        else:
            return 0

    def setClassIndex(self,ciString):
        if re.search(',',ciString):
            (self.eltClass,self.index) = ciString.rsplit(',')
        else:
            return 0

# make sure annotation information exists
def checkAnnotDir(refGenome,annotDir):
    if not os.path.exists(annotDir + "/" + refGenome):
        raise IOError("cannot find genome annotation directory for " + refGenome)
        return 0
    discordant.checkfile(annotDir + "/" + refGenome + "/names.txt")

# ins = insertionSummary object
def annotatePos(insList,refGenome,annotDir,outfileName,printout):
    namefile = open(annotDir + "/" + refGenome + "/names.txt", 'r')
    outfile  = open(outfileName, 'w')
 
    # set up each annotator (name, tabix file) present
    annDict = {}
    for nameline in namefile:
        (tname,tfile) = nameline.rstrip().split("\t")
        tfile = annotDir + "/" + refGenome + "/" + tfile
        annotator = annotate.annotator(tfile, tname, refGenome)
        annDict[tname] = annotator
    namefile.close()

    for ins in insList:
        outstring = str(ins)
        for annName,ann in annDict.iteritems():
            bedline = "\t".join((ins.chr,ins.start,ins.end))
            if ins.start > ins.end:
                bedline = "\t".join((ins.chr,ins.end,ins.start))
            annList = ann.annotate(bedline,refGenome)
            annJoin = ",".join(annList)
            ins.annot.append(annJoin)
        if printout:
            print ins
        outfile.write(str(ins) + "\n")
    outfile.close()

def main(args):
    peakparser.checkOutDir(args.outBaseName,args.outDirName)

    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    discordant.checkfile(configPath)
    configDict = peakparser.readConfig(configPath,args.outBaseName,args.outDirName)

    refGenome = configDict['refGenome']
    checkAnnotDir(refGenome, args.annotDir)

    debugFile = args.outDirName + "/" + args.outBaseName + "/" + args.outBaseName + ".debug.txt"
    discordant.checkfile(debugFile)

    # can filter on retroelement subfamilies
    discordant.checkfile(args.eltFile)
    eltDict = {}
    f = open(args.eltFile,'r')
    for elt in f:
        eltDict[elt.strip()]=1

    germlineIns = []
    cancerIns   = []
    normalIns   = []
    otherIns    = []

    f = open(debugFile, 'r')
    for line in f: # loop over lines
        if re.search('valid=1', line) and re.search('olClN=0', line):
            insdata = {}
            for c in line.rstrip().rsplit(' '): # loop over columns
               (key,value) = c.rsplit('=')
               insdata[key] = value

            if int(insdata['np']) >= int(args.minPeakSize) and eltDict.has_key(insdata['pME']): 
                insSum = insertionSummary()
                insSum.setPos(insdata['pos'])
                insSum.setEltPos(insdata['eltExt'])
                insSum.setBestStrand(insdata['eInv'],insdata['eS'])
                insSum.setSources(insdata['sources'])
                insSum.setClassIndex(insdata['index'])
                insSum.eltFam   = insdata['pME']
                insSum.numreads = insdata['np']
                if len(insSum.sources) > 1:
                    germlineIns.append(insSum)
                else:
                    if re.search('CANCER', insSum.sources[0]):
                        cancerIns.append(insSum)
                    elif re.search('NORMAL', insSum.sources[0]):
                        normalIns.append(insSum)
                    else:
                        otherIns.append(insSum)
    f.close()
    annotatePos(germlineIns,refGenome,args.annotDir,args.outDirName + "/" + args.outBaseName + "/germline.tab.txt",args.printout)
    annotatePos(cancerIns,refGenome,args.annotDir,args.outDirName + "/" + args.outBaseName + "/canceronly.tab.txt",args.printout)
    annotatePos(normalIns,refGenome,args.annotDir,args.outDirName + "/" + args.outBaseName + "/normalonly.tab.txt",args.printout)
    annotatePos(otherIns,refGenome,args.annotDir,args.outDirName + "/" + args.outBaseName + "/other.tab.txt",args.printout)

if __name__ == '__main__':
    # commandline args
    parser = argparse.ArgumentParser(description='parse the output of discordant.py')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='config file left by discordant.py')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for output files')
    parser.add_argument('-e', '--eltfile', dest='eltFile', default='sumEltList.txt',
                        help='list of element families to include')
    parser.add_argument('-m', '--minpeak', dest='minPeakSize', default=8,
                        help='minimum reads per peak (default 8)')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-a', '--annotdir', dest='annotDir', required=True,
                        help='base directory for annotations (annotation/refgenome should exist too)')
    parser.add_argument('--printout', action="store_true")
    args = parser.parse_args()

    main(args)

