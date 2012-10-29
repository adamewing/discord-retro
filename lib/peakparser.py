#!/bin/env python

#
# Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)
#
# Released under the MIT license, see LICENSE.txt
#

import pysam, argparse, pairinfo, discordant, sys, os, math, logging

# facilitates keeping reads from different element classes seperate
class PeakBuilder:
    def __init__(self,eltClass):
        self.eltClass = eltClass
        self.peak     = pairinfo.Peak()
        self.peaks    = []
        self.lastPeakNum = 0

def makeEltDict(fname):
    eltfile = open(fname,'r')
    eltdict = {}

    for eltline in eltfile:
        if eltline[0] != '#':
            (eltclass,eltname) = eltline.strip().split()
            eltdict[eltname] = eltclass
    return eltdict

# function for checking coordinates against repeatmasker tabix files
# if eltClass is set fo 'all', all entries in the tabix fill that overlap
# will generate a match. Returns True if an element overlapped
def checkRmskDict(tabixFile,eltDict,eltClass,chr,start,end):
    if start > end: # can happen in messy invalid peaks
       tmp   = start
       start = end
       end   = tmp
    for tabixRecord in tabixFile.fetch(reference=chr, start=start, end=end):
        t = tabixRecord.strip().split()
        eltName = t[3]
        if eltClass == 'All': # match anything
            return True
        else:
            if eltDict.has_key(eltName) and eltDict[eltName] == eltClass:
                return True
    return False

def longOutput(eltPeakDict,outBaseName,outDirName,tabixFile,eltDict):
    fname = outBaseName + ".debug.txt"
    f = open(outDirName + "/" + outBaseName + "/" + fname, 'w')
    
    pindex = {}
    for eltClass in eltPeakDict.keys():
        pindex[eltClass] = 0
        for peak in eltPeakDict[eltClass].peaks:
            if peak.hasPairs(): # an empty peak would be a bad sign, but ignore that for now...
                breakpoint = peak.getBreakpoint()

                # annotate overlaps with repeatmasker elements
                overlapDictEltWide   = checkRmskDict(tabixFile,eltDict,eltClass,peak.getChr(),
                                                     breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                overlapAnyEltWide    = checkRmskDict(tabixFile,eltDict,'All',peak.getChr(),
                                                     breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                overlapDictEltNarrow = checkRmskDict(tabixFile,eltDict,eltClass,peak.getChr(),
                                                     breakpoint.leftpos,breakpoint.rightpos)
                overlapAnyEltNarrow  = checkRmskDict(tabixFile,eltDict,'All',peak.getChr(),
                                                     breakpoint.leftpos,breakpoint.rightpos)

                outstr = ("pos=" + peak.getChr() + ":%d-%d" % (breakpoint.leftpos,breakpoint.rightpos)
                       + " str=" + breakpoint.leftstrand + "," + breakpoint.rightstrand + " valid=%d" % breakpoint.isValid()
                       + " lmd=%d" % breakpoint.leftMeanDist() + " rmd=%d" % breakpoint.rightMeanDist()
                       + " lws=%d" % breakpoint.leftWrongStrand() + " rws=%d" % breakpoint.rightWrongStrand()
                       + " lMws=%d" % breakpoint.leftMateWrongStrand() + " rMws=%d" % breakpoint.rightMateWrongStrand()
                       + " np=%d" % breakpoint.npairs + " intD=%d" % (breakpoint.rightpos-breakpoint.leftpos)
                       + " eltExt=%d,%d" % (breakpoint.eltMinPos(),breakpoint.eltMaxPos())
                       + " pME=" + peak.primaryMateElt() + " olClW=%d" % overlapDictEltWide + " olAnyW=%d" % overlapAnyEltWide
                       + " olClN=%d" % overlapDictEltNarrow + " olAnyN=%d" % overlapAnyEltNarrow 
                       + " eS=" + breakpoint.eltStrandFromReads() + " sources=" + ','.join(breakpoint.sourceList())
                       + " eInv=" + breakpoint.eltStrandFromStrand() + " index=" + eltClass + ",%d" % pindex[eltClass] + "\n")
                f.write(outstr)
            pindex[eltClass] += 1
    f.close()

def bedOutput(eltPeakDict,outBaseName,outDirName,refGenome,tabixFile,eltDict):
    callsFileName = outBaseName + ".calls.bed"
    callsFile = open(outDirName + "/" + outBaseName + "/" + callsFileName, 'w')

    # write headers
    callsHeader = "track\tname=" + outBaseName + "_Calls\tvisibility=2\tuseScore=1\titemRgb=On\tdb=" + refGenome + "\n"
    callsFile.write(callsHeader)

    pindex = {}
    for eltClass in eltPeakDict.keys():
        pindex[eltClass] = 0
        for peak in eltPeakDict[eltClass].peaks:
            if peak.hasPairs(): # an empty peak would be a bad sign, but ignore that for now...
                breakpoint = peak.getBreakpoint()

                # annotate overlaps with repeatmasker elements
                overlapDictEltWide   = checkRmskDict(tabixFile,eltDict,eltClass,peak.getChr(),
                                                     breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                overlapAnyEltWide    = checkRmskDict(tabixFile,eltDict,'All',peak.getChr(),
                                                     breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                overlapDictEltNarrow = checkRmskDict(tabixFile,eltDict,eltClass,peak.getChr(),
                                                     breakpoint.leftpos,breakpoint.rightpos)
                overlapAnyEltNarrow  = checkRmskDict(tabixFile,eltDict,'All',peak.getChr(),
                                                     breakpoint.leftpos,breakpoint.rightpos)

                # deduct from the bed 'score' if 
                bedscore = 1000
                if overlapAnyEltWide:
                    bedscore = bedscore - 200
                if overlapDictEltWide:
                    bedscore = bedscore - 200
                if overlapAnyEltNarrow:
                    bedscore = bedscore - 200
                if overlapDictEltNarrow:
                    bedscore = bedscore - 200
                if breakpoint.eltStrandFromReads() == 'x':
                    bedscore = bedscore - 100
                if breakpoint.eltStrandFromStrand() == 'x':
                    bedscore = bedscore - 100 

                bedname =  eltClass + ",%d" % pindex[eltClass]

                bedline = (peak.getChr() + "\t%d\t%d" % (breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                        + "\t" + bedname + "\t%d" % bedscore + "\t" + breakpoint.bestStrand()
                        + "\t%d\t%d" % (breakpoint.leftpos,breakpoint.rightpos) + "\n")

                # using bedscore here is debateable, but in general calls with score < 200 are shaky
                if breakpoint.isValid() and bedscore > 200:
                    callsFile.write(bedline)
            pindex[eltClass] = pindex[eltClass] + 1
    callsFile.close()

def vcfOutput(eltPeakDict,outBaseName):
    pass

def checkOutDir(outBaseName,outDirName):
    if not os.path.exists(outDirName + "/" + outBaseName):
        raise IOError('cannot find output directory (' + outDirName + ')')
    if not os.path.exists(outDirName + "/" + outBaseName):
        raise IOError("cannot find sample directory (" + outDirName + "/" + outBaseName + ")")
    if not os.path.exists(outDirName + "/" + outBaseName + '/logs'):
        os.mkdir(outDirName + "/" + outBaseName + '/logs')

def readConfig(configPath,outBaseName, outDirName):
    f = open(configPath, 'r')
    configDict = dict()
    for configline in f:
        (param, value) = configline.rstrip().rsplit("=")
        configDict[param] = value

    # check for minimum necessary parameters and some sanity checks
    if configDict['outBaseName']:
        if configDict['outBaseName'] != outBaseName:
            raise ValueError("baseName mismatch: " + outBaseName + "=/="
                             + configDict['outBaseName'])
    else:
        raise ValueError("outBaseName is not set in " + configPath)
    if not configDict['tabixFileName']:
        raise ValueError("tabixFileName is not set in " + configPath)
    if not configDict['eltFileName']:
        raise ValueError("eltFileName is not set in " + configPath)

    configDict['readFileName'] = outDirName + "/" + outBaseName + "/" + outBaseName + ".readpairs.txt"

    return configDict

def main(args):
    checkOutDir(args.outBaseName,args.outDirName)
    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    discordant.checkfile(configPath)
    configDict = readConfig(configPath,args.outBaseName,args.outDirName)

    # debugging info
    logfile = args.outDirName + "/" + args.outBaseName + "/logs/%d" % os.getpid() + "." + args.outBaseName +".peakparser.log"
    logging.basicConfig(format='%(asctime)s %(message)s',filename=logfile,level=logging.DEBUG)

    # log parameters
    logging.info("\noutBaseName=%s\nbamFileName=%s\ntabixFileName=%s\neltFileName=%s"
                 % (args.outBaseName,configDict['readFileName'],configDict['tabixFileName'],
                    configDict['eltFileName']))

    discordant.checkfile(configDict['readFileName'])
    discordant.checkfile(configDict['tabixFileName'])
    discordant.checkfile(configDict['eltFileName'])

    tabixFile = pysam.Tabixfile(configDict['tabixFileName'], 'r')
    eltDict = makeEltDict(configDict['eltFileName'])

    # structure for keeping element classes seperate
    eltPeakDict = {}
    for (eltName,eltClass) in eltDict.iteritems():
        eltPeakDict[eltClass] = PeakBuilder(eltClass)

    f = open(configDict['readFileName'])
    for line in f:
        fields = line.rsplit()
        eltClass = ''
        try: 
            eltClass = fields[13].rstrip()
        except:
            print "bad line: " + line
        chrom      = fields[0]
        readPos    = int(fields[1])
        readStrand = fields[2]
        mateElt    = fields[3]
        mateChrom  = fields[4]
        mateStrand = fields[5]
        matePos    = int(fields[6])
        eltStart   = int(fields[7])
        eltEnd     = int(fields[8])
        eltStrand  = fields[9]
        eltFullLen = fields[10]
        genomeName = fields[11]
        peakIndex  = fields[12]

        pair = pairinfo.Pair(chrom,readPos,readStrand,mateElt,mateChrom,
                             mateStrand,matePos,eltStart,eltEnd,eltStrand,
                             eltFullLen,genomeName,peakIndex,eltClass)
        if eltPeakDict[eltClass].lastPeakNum != peakIndex:
            eltPeakDict[eltClass].peaks.append(eltPeakDict[eltClass].peak)
            eltPeakDict[eltClass].peak = pairinfo.Peak()
            eltPeakDict[eltClass].peak.addpair(pair)
            eltPeakDict[eltClass].lastPeakNum = peakIndex
        else:
            eltPeakDict[eltClass].peak.addpair(pair)

    logging.info("starting long output...")
    longOutput(eltPeakDict,args.outBaseName,args.outDirName,tabixFile,eltDict)

    logging.info("starting bed output...")
    bedOutput(eltPeakDict,args.outBaseName,args.outDirName,configDict['refGenome'],tabixFile,eltDict)

if __name__ == '__main__':
    # commandline args
    parser = argparse.ArgumentParser(description='parse the output of discordant.py')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='config file left by discordant.py')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for output files')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('--permissive', action="store_true")
    args = parser.parse_args()

    main(args)

