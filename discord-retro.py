#!/bin/env python
#
# Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)
#
# Released under the MIT license, see LICENSE.txt
#

import sys
import re
import time
from os import path as path
import pp
import argparse
import subprocess
import pysam
import ConfigParser

import lib.configtest as configtest
# other imports are done via pp

def asmFromPhrase(phrase):
    if (re.search('build36', phrase) or
        re.search('HG36', phrase) or
        re.search('hg36', phrase) or
        re.search('assembly18', phrase) or
        re.search('HG18', phrase) or
        re.search('hg18', phrase) or
        re.search('GRCh36', phrase)):
        return 'hg18'

    if (re.search('build37', phrase) or
        re.search('HG37', phrase) or
        re.search('hg37', phrase) or
        re.search('assembly19', phrase) or
        re.search('HG19', phrase) or
        re.search('hg19', phrase) or
        re.search('GRCh37', phrase)):
        return 'hg19'

class PairedSample:
    def __init__(self,name,config,outdir,tcga=False, pinonly=False):
        self.name = name
        self.tcga = tcga      # if true, changes validation, naming
        self.pinonly   = pinonly
        self.config    = config  # ConfigParser
        self.outdir    = outdir
        self.cancerBam = None
        self.normalBam = None
        self.cancerIdx = None
        self.normalIdx = None 
        self.assembly  = None

    def __str__(self):
        return "\t".join((str(self.cancerBam),str(self.normalBam),
                          str(self.cancerIdx),str(self.normalIdx)))

    def validate(self):
        if self.tcga:
            if (self.cancerBam and path.exists(self.cancerBam) and
                self.normalBam and path.exists(self.normalBam) and
                self.cancerIdx and path.exists(self.cancerIdx) and
                self.normalIdx and path.exists(self.normalIdx)):
                return True
            else:
                return False
        else:
            if (self.normalBam and path.exists(self.normalBam) and
                self.normalIdx and path.exists(self.normalIdx)):
                return True
            else:
                return False

    def addFile(self,sampleType,extension,filePath):
        if sampleType == 'NORMAL' and extension == 'bam':
            self.normalBam = filePath
            self.normalIdx = filePath + ".bai"
        if sampleType == 'CANCER' and extension == 'bam':
            self.cancerBam = filePath
            self.cancerIdx = filePath + ".bai"

    def getAssembly(self,default):
        """
        use pysam to try and determine whether this is hg18 or hg19
        must be run before runDiscordant
        """
        try:
            bamFile = pysam.Samfile(self.normalBam, 'rb')
            header = bamFile.header
            asm = ''
            asmphrase = ''
            if header.has_key('SQ'):
                sqdict = header['SQ'][0]
                if sqdict.has_key('AS'):
                    asm = asmFromPhrase(sqdict['AS'])
                elif sqdict.has_key('UR'):
                    asm = asmFromPhrase(sqdict['UR'])
                else:
                    asm = default
            else:
                asm = default

            self.assembly = asm
            return asm
        except:
            sys.stderr.write("unrecognized entry in header, assuming " + default + "\n")
            self.assembly = default
            return default

    def runDiscordant(self):
        """Runs discordant.py, peakparser.py, summarize,py, pinpoint.py"""
        tabixLoc = self.config.get('discord', 'annotDir') + "/" + self.assembly + "/" + self.assembly + ".rmsk.bed.gz"
        lib.discordant.checkfile(tabixLoc)
        usechr = False
        if self.config.get('discord','usechr') == "True":
            usechr = True
        args = argparse.Namespace(bamFileName    = self.normalBam,  # discordant.py
                                  outBaseName    = self.name,       # discordant.py, peakparser.py, summarize.py, pinpoint.py
                                  tabixFileName  = tabixLoc,        # discordant.py
                                  refGenome      = self.assembly,   # discordant.py
                                  overwrite      = True,            # discordant.py
                                  printout       = False,           # summarize.py
                                  inDir1         = None,            # mergepairs.py
                                  inDir2         = None,            # mergepairs.py
                                  outDirName     = self.outdir,     # everything
                                  usechr         = usechr,          # pinpoint.py
                                  eltFileName    = self.config.get('discord','eltFileName'),
                                  eltLenFileName = self.config.get('discord','eltLenFileName'),
                                  readLength     = self.config.get('discord','readLength'),
                                  insertSize     = self.config.get('discord','insertSize'),
                                  minMapQ        = self.config.get('discord','minMapQ'),
                                  configFileName = self.config.get('discord','configFileName'),
                                  eltFile        = self.config.get('discord','eltFile'),
                                  minPeakSize    = self.config.get('discord','minPeakSize'),
                                  maxReadLen     = self.config.get('discord','maxReadLen'),
                                  zeroChar       = self.config.get('discord','zeroChar'),
                                  minClipQual    = self.config.get('discord','minClipQual'),
                                  refFastaDir    = self.config.get('discord','refFastaDir'),
                                  annotDir       = self.config.get('discord','annotDir'),
                                  refGenomeFile  = self.config.get('discord',self.assembly))
        if self.tcga:
            normalName = self.name + "-NORMAL"
            cancerName = self.name + "-CANCER"
            mergeName  = self.name + "-MERGE"

            args.inDir1 = normalName
            arge.inDir2 = cancerName
            args.outBaseName = normalName
            # normal .bam
            lib.discordant.main(args)

            # cancer .bam
            args.bamFileName = self.cancerBam
            args.outBaseName = cancerName
            lib.discordant.main(args)

            # merged results
            args.outBaseName = mergeName
            lib.mergepairs.main(args)

            if not self.pinonly:
                lib.peakparser.main(args)
                lib.summarize.main(args)
            lib.pinpoint.main(args)

        else:
            if not self.pinonly:
                lib.discordant.main(args)
                lib.peakparser.main(args)
                lib.summarize.main(args)
            lib.pinpoint.main(args)

def main(args):
    configtest.check(args.configFile)
    
    pairedSamples = {}

    for line in open(args.sampleFile, 'r'):
        if not re.search("^#", line):
            (filePath,label) = line.strip().split()
            base = path.basename(filePath)

            sampleType = 'NORMAL'
            sampleName = label
            fileparts = base.split('.')
            extension = fileparts[-1]

            if args.tcga: # see https://wiki.nci.nih.gov/display/TCGA/Working+with+TCGA+Data
                baseparts = base.split('-')
                participant = baseparts[2]
                samplenum = baseparts[3]
                if re.search('01', samplenum):
                    sampleType = 'CANCER'
                sampleName = '-'.join((label,participant))

            config = ConfigParser.ConfigParser()
            config.read(args.configFile)

            if sampleName not in pairedSamples:
                pairedSamples[sampleName] = PairedSample(sampleName,config, 
                                                         args.outDirName,
                                                         tcga=args.tcga, 
                                                         pinonly=args.pinpointonly)

            pairedSamples[sampleName].addFile(sampleType,extension,filePath)

    # assign assemblies from header information
    print "trying to determine assemblies from .bam file headers..."
    for sampleName in pairedSamples.keys():
        if pairedSamples[sampleName].validate():
            asm = pairedSamples[sampleName].getAssembly(args.defaultRef)
            print sampleName + " is " + asm

    # parallel python stuff
    ncpus = 16
    jobServer = pp.Server(ncpus,ppservers=())

    sampleJobs = {}
    for sampleName in pairedSamples.keys():
        if pairedSamples[sampleName].validate():
            print sampleName + " validated"
            job = jobServer.submit(pairedSamples[sampleName].runDiscordant,(),(),
                  ("lib.discordant","lib.mergepairs","lib.peakparser","lib.pinpoint","lib.summarize","argparse"))
            sampleJobs[sampleName] = job
            print "started job: " + sampleName
            jobServer.print_stats()
        else:
            print sampleName + " did not validate, make sure the .bam exists and is indexed (has a .bam and a .bam.bai file)"

    for sampleName in sampleJobs:
        print "=== stdout of job: " + sampleName + " ==="
        print sampleJobs[sampleName]()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pipeline that aims to identify TE insertions from .bam alignments")
    parser.add_argument('-s', '--sampleFile', dest="sampleFile", required=True,
                    help="List of filenames and sample names. Should include both .bam files (.bam) and their indexes (.bai)")
    parser.add_argument('-r', '--defaultReference', dest="defaultRef", default='hg19',
                    help="Default reference assembly. If assembly can be determined from .bam header that will trump the default.")
    parser.add_argument('-c', '--config', dest='configFile', required=True)
    parser.add_argument('-p', '--processors', dest='numCPUs', default='1', help="Number of CPUs to use (1 job per CPU)")
    parser.add_argument('-o', '--outdir', dest='outDirName', required=True, help="output base directory")
    parser.add_argument('--pinpointonly', action='store_true', help="only run pinpoint.py")
    parser.add_argument('--tcga', action="store_true", help="Use TCGA filename format to determine cancer/normal pair")
    args = parser.parse_args()
    main(args)
