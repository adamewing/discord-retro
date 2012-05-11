#!/bin/env python

#
# Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)
#
# Released under the MIT license, see LICENSE.txt
#

import pysam, argparse, sys, pairinfo, logging, os, re

def checkfile(fname):
    try:
        open(fname)
    except IOError as e:
        print "can't find file: " + fname
        sys.exit()

# reads the element dictionary that specifies which retroelement annotations belong
# to a common class of element. e.g. L1PA2, L1P1, and L1MA5 are all LINEs whereas
# AluYa5, AluYb8, and AluY are all SINE annotations. This allows fine-grained control
# of which subfamilies are considered

def makeEltDict(fname):
    eltfile = open(fname,'r')
    eltdict = dict()

    for eltline in eltfile:
        if eltline.strip()[0] != '#':
            (eltclass,eltname) = eltline.strip().split()
            eltdict[eltname] = eltclass
    eltfile.close()
    return eltdict

def makeEltLenDict(fname):
    eltfile = open(fname,'r')
    eltdict = dict()

    for eltline in eltfile:
        if eltline.strip()[0] != '#':
            (eltname,eltlen) = eltline.strip().split()
            eltdict[eltname] = int(eltlen)
    eltfile.close()
    return eltdict

def prepOutDir(outBaseName,outDirName,overwrite):
    if not os.path.exists(outDirName):
        os.mkdir(outDirName)
    if not os.path.exists(outDirName + "/" + outBaseName):
        os.mkdir(outDirName + "/" + outBaseName)
    else:
        if not overwrite:
            sys.exit(outDirName + " exists, can be overridden with --overwrite")
    if not os.path.exists(outDirName + "/" + outBaseName + "/logs"):
        os.mkdir(outDirName + "/" + outBaseName + "/logs")

def writeConfig(args):
    f = open(args.outDirName + "/" + args.outBaseName + "/config.txt", 'w')
    for k,v in vars(args).iteritems():
         f.write(k + "=" + str(v) + "\n")
    f.close()

# fix for the most common variation on chromosome names (leaving out the 'chr')
def fixChrName(name):
    if name[0] != 'c':
        return "chr" + name
    else:
        return name

class clusterBuild:
    """Stores placeholder data for building clusters"""
    def __init__(self,eltClass):
        self.lastPos   = 0
        self.lastChr   = ''
        self.peakIndex = 0

# used for generating new index numbers
def reindexReads(readlist, maxDist, eltLenDict):
    # used for keeping track of the last readpair seen for a given element class
    eltClustDict = {}
    for eltName in eltLenDict.keys():
        eltClustDict[eltName] = clusterBuild(eltName)

    for readline in readlist:
        (readChrName,readPos,readStrand,mateElt,mateChrom,
         mateStrand,matePos,eltStart,eltEnd,eltStrand,
         eltFullLen,source,oldIndex,eltClass) = readline

        readPos = int(readPos)

        if (readChrName == eltClustDict[eltClass].lastChr):
            if (readPos - eltClustDict[eltClass].lastPos < 0): # sanity check
                raise IndexError('.bam file does not appear to be properly sorted')
            if (readPos - eltClustDict[eltClass].lastPos < maxDist): # keep the same peak number
                eltClustDict[eltClass].lastPos = readPos
            else: # create a new peak and add the pair to it if pair is > maxDist from last one
                eltClustDict[eltClass].peakIndex = eltClustDict[eltClass].peakIndex + 1
                eltClustDict[eltClass].lastPos = readPos
        else: # create a new peak and add the pair to it if pair is on a new chromosome
            eltClustDict[eltClass].peakIndex = eltClustDict[eltClass].peakIndex + 1
            eltClustDict[eltClass].lastChr = readChrName
            eltClustDict[eltClass].lastPos = readPos

        newlist.append((readChrName,str(readPos),readStrand,mateElt,mateChrom,
                        mateStrand,matePos,eltStart,eltEnd,eltStrand,eltFullLen,
                        source,str(eltClustDict[eltClass].peakIndex),eltClass))
    return newlist

def main(args):
    prepOutDir(args.outBaseName,args.outDirName,args.overwrite)
    writeConfig(args) # save arguments in output directory

    # debugging info
    logfile = args.outDirName + "/" + args.outBaseName + "/logs/%d" % os.getpid() + "." + args.outBaseName +".discordant.log"
    logging.basicConfig(format='%(asctime)s %(message)s',filename=logfile,level=logging.DEBUG)

    # log parameters
    logging.info("\noutBaseName=%s\nbamFileName=%s\ntabixFileName=%s\neltFileName=%s\nreadLength=%s\ninsertSize=%s" 
             % (args.outBaseName,args.bamFileName,args.tabixFileName,args.eltFileName,args.readLength,args.insertSize))

    chromFile = 'chromlist.txt'

    checkfile(args.bamFileName)
    checkfile(args.tabixFileName)
    checkfile(args.eltFileName)

    bamFile = pysam.Samfile(args.bamFileName, 'rb') # rb = read, binary
    tabixFile = pysam.Tabixfile(args.tabixFileName, 'r')
    tabixRefs = {}
    for ref in tabixFile.contigs:
        tabixRefs[ref] = 1

    eltDict = makeEltDict(args.eltFileName)
    chrDict = {}
    for ref in bamFile.references:
        chrDict[ref] = 1 

    # Dictionary of element lengths (hard coded right now)
    # FIXME allow this to be specified externally, maybe use this as default
    # eltLenDict = {'LINE': 6019, 'SINE': 282, 'Other':2000, 'LTR':950}
    eltLenDict = makeEltLenDict(args.eltLenFileName)

    # used for keeping track of the last readpair seen for a given element class
    eltClustDict = {}
    for eltName in eltLenDict.keys():
        eltClustDict[eltName] = clusterBuild(eltName)

    # parameters for reads
    readLength = int(args.readLength)
    insertSize = int(args.insertSize)
    maxDist = readLength + (2*insertSize)

    bamIndex = 0
    discordIndex = 0

    # output file
    readsFileName = args.outBaseName + ".readpairs.txt"
    bedFileName   = args.outBaseName + ".reads.bed"
    readsFile = open(args.outDirName + "/" + args.outBaseName + "/" + readsFileName, 'w')
    bedFile   = open(args.outDirName + "/" + args.outBaseName + "/" + bedFileName, 'w')

    bedHeader = "track\tname=" + args.outBaseName + "_Reads\tvisibility=2\titemRgb=On\tdb=" + args.refGenome + "\n"
    bedFile.write(bedHeader)

    for read in bamFile.fetch():
    #if maxLines > 0 and maxLines < bamIndex:
        #break
    # avoid putting anything else before this 'if statement' as it will be evaluated for _every_ read
    # in the .bam file (and there are a lot of reads...)

        if (not read.is_proper_pair and not read.is_duplicate 
            and not read.mate_is_unmapped and not read.is_unmapped
            and read.tid > -1 and read.mrnm > -1
            and chrDict.has_key(bamFile.getrname(read.tid))
            and chrDict.has_key(bamFile.getrname(read.mrnm))):

            discordIndex += 1
            readElt = ''
            mateElt = ''

            readChrName = fixChrName(bamFile.getrname(read.tid))
            mateChrName = fixChrName(bamFile.getrname(read.mrnm))

            tabixTupleParse = None
            tabixMateTupleParse = None

            # fetch() using the pysam.asTuple parser returns tabix results as a tuple
            if tabixRefs.has_key(readChrName) and tabixRefs.has_key(mateChrName):
                tabixTupleParse = tabixFile.fetch(reference=readChrName, 
                                                  start=read.pos, 
                                                  end=read.pos+1, 
                                                  parser=pysam.asTuple())

                tabixMateTupleParse = tabixFile.fetch(reference=mateChrName, 
                                                      start=read.mpos, 
                                                      end=read.mpos+1, 
                                                      parser=pysam.asTuple())

            if tabixTupleParse:
                for tabixTuple in tabixTupleParse:
                    readElt = tabixTuple[3]
                for tabixMateTuple in tabixMateTupleParse:
                    mateElt = tabixMateTuple[3]

                # For a given read pair, add it to the set if a read's mate matches an 
                # annotation in the element dictionary, and the read itself is either 
                # not a repeat or is not a repeat in the dictionary or is a repeat in 
                # the dictionary but not the same type as the mate read. Mapping quality
                # must be >= 20 (this seems to filter out a significant portion of the noise)

                if (eltDict.has_key(mateElt) and (not readElt or not eltDict.has_key(readElt) or 
                   (eltDict.has_key(readElt) and eltDict[mateElt] != eltDict[readElt])) and
                   (read.mapq >= int(args.minMapQ))):

                    eltClass   = eltDict[mateElt]
                    readStrand = '+'
                    mateStrand = '+'
                    if read.is_reverse:      readStrand = '-'
                    if read.mate_is_reverse: mateStrand = '-'

                    eltStart   = int(tabixMateTuple[1])
                    eltEnd     = int(tabixMateTuple[2])
                    eltDiverge = tabixMateTuple[4]
                    eltStrand  = tabixMateTuple[5]

                    if (readChrName == eltClustDict[eltClass].lastChr):
                        if (read.pos - eltClustDict[eltClass].lastPos < 0): # sanity check
                            raise IndexError('.bam file does not appear to be properly sorted')
                        if (read.pos - eltClustDict[eltClass].lastPos < maxDist): # keep the same peak number
                            eltClustDict[eltClass].lastPos = read.pos
                        else: # create a new peak and add the pair to it if pair is > maxDist from last one
                            eltClustDict[eltClass].peakIndex = eltClustDict[eltClass].peakIndex + 1
                            eltClustDict[eltClass].lastPos = read.pos
                    else: # create a new peak and add the pair to it if pair is on a new chromosome
                        eltClustDict[eltClass].peakIndex = eltClustDict[eltClass].peakIndex + 1
                        eltClustDict[eltClass].lastChr = readChrName
                        eltClustDict[eltClass].lastPos = read.pos

                    readline = (readChrName + "\t%d" % read.pos
                             + "\t" + readStrand + "\t" + mateElt + "\t"
                             + mateChrName + "\t"
                             + mateStrand + "\t%d" % read.mpos + "\t%d"
                             % eltStart + "\t%d" % eltEnd + "\t" + eltStrand
                             +"\t%d" % eltLenDict[eltClass] + "\t" 
                             + args.outBaseName + "\t%d" % eltClustDict[eltClass].peakIndex 
                             + "\t" + eltClass + "\n")

                    # make a pair object to compute things like relative position in element, etc.
                    pair = pairinfo.Pair(readChrName,read.pos,readStrand,mateElt,mateChrName,mateStrand,
                                         read.mpos,eltStart,eltEnd,eltStrand,eltLenDict[eltClass],
                                         args.outBaseName,eltClustDict[eltClass].peakIndex,eltClass)

                    bedcolor = "0,0,0"
                    if re.search("CANCER", args.outBaseName):
                        bedcolor = "255,0,0"
                    if re.search("NORMAL", args.outBaseName):
                        bedcolor = "0,0,255"
 
                    bedname = eltClass + "_%d" % pair.relPos()
                    readEnd = read.pos + int(args.readLength)
                    bedline = (readChrName + "\t%d" % read.pos + "\t%d" % readEnd + "\t" + bedname
                            + "\t1000\t" + readStrand + "\t%d" % read.pos + "\t%d" % readEnd
                            + "\t" + bedcolor + "\n")

                    readsFile.write(readline)
                    bedFile.write(bedline)

        bamIndex = bamIndex + 1
    readsFile.close()
    bedFile.close()
    logging.info("%s finished, (%d reads, %d discordant)" % (args.bamFileName,bamIndex,discordIndex))

# commandline args
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse .bam file and output discordant reads where one end is in an annotation specified in a tabix index')
    parser.add_argument('-b', '--bam', dest='bamFileName', required=True,
                        help='name of .bam file')
    parser.add_argument('-t', '--tabix', dest='tabixFileName', required=True,
                        help='name of RepeatMasker tabix index (.gz.tbi)')
    parser.add_argument('-e', '--eltfile', dest='eltFileName', required=True,
                        help='name of file containing class/family list')
    parser.add_argument('-l', '--eltlen', dest='eltLenFileName', required=True,
                        help='list of element maximum lengths')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', default='defaultName',
                        help='basename for output files')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-r', '--readlen', dest='readLength', default='100',
                        help='read length in basepairs (default  100 bp)')
    parser.add_argument('-i', '--insertsize', dest='insertSize', default='300',
                        help='expected insert size in bp (default 300 bp)')
    parser.add_argument('-g', '--refgenome', dest='refGenome', default='hg19',
                        help='reference genome assembly (defaults to hg19)')
    parser.add_argument('-q', '--minqual', dest='minMapQ', default='20',
                        help='minimum mapping quality')

    parser.add_argument('--overwrite', action="store_true")
    args = parser.parse_args()

    main(args)
