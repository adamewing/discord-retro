#!/bin/env python

#
# Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)
#
# Released under the MIT license, see LICENSE.txt
#

import math

class Pair:
    """Stores information for a genome/element read pair"""
    def __init__(self,chrom,readPos,readStrand,mateElt,mateChrom,mateStrand,matePos,eltStart,eltEnd,eltStrand,eltFullLen,source,peakIndex,eltClass):
        self.chrom      = chrom
        self.readPos    = readPos
        self.readStrand = readStrand
        self.mateElt    = mateElt
        self.mateChrom  = mateChrom
        self.mateStrand = mateStrand
        self.matePos    = matePos
        self.eltStart   = eltStart
        self.eltEnd     = eltEnd
        self.eltStrand  = eltStrand
        self.eltFullLen = eltFullLen
        self.source     = source # keep track of genome
        self.peakIndex  = peakIndex
        self.eltClass   = eltClass
        if matePos < eltStart-500 or matePos > eltEnd+500:
            raise ValueError('matePos outside of element boundaries: matepos=' + str(matePos) + " eltStart=" + str(eltStart) + " eltEnd=" + str(eltEnd))

    def __str__(self):
        output = "\t".join((self.chrom, str(self.readPos), self.readStrand, self.mateElt, self.mateChrom, self.mateStrand, str(self.matePos),
                           str(self.eltStart), str(self.eltEnd), self.eltStrand, str(self.eltFullLen), self.source, str(self.peakIndex), self.eltClass))
        return output

    # return the strand of the ref element relative to the alignment
    def relStrand(self):
        if self.eltStrand == self.mateStrand:
            return '+'
        else:
            return '-'

    # return the position of the mate read relative to the ref element coords
    def relPos(self):
        refLen   = self.eltEnd - self.eltStart
        refLeft  = self.matePos - self.eltStart
        refRight = self.eltEnd - self.matePos
        if self.eltStrand == '+':
            if (int(self.eltFullLen) - refRight) < 0:
                return 0
            else:
                return int(self.eltFullLen) - refRight
        else: # - strand
            if (int(self.eltFullLen) - refLeft) < 0:
                return 0
            else:
                return int(self.eltFullLen) - refLeft

class Peak:
    """Stores objects of class Pair and provides functions for analysis"""
    def __init__(self):
        self.pairs = []
        self.npairs = 0 
    def addpair(self,pair):
        self.pairs.append(pair)
        self.npairs = self.npairs + 1
    def __iter__(self):
        return self
    def next(self):
        if self.npairs == 0:
            self.index = len(self.pairs)
            raise StopIteration
        self.npairs = self.npairs - 1
        return self.pairs[self.npairs]

    # find breakpoint (where ref strands switch):
    def getBreakpoint(self):
        plusScore  = [0] * len(self.pairs) 
        minusScore = [0] * len(self.pairs)
        index = 0
        for pair in self.pairs:
            if pair.readStrand == '+':
                if index == 0:
                    plusScore[0]  = 1
                    minusScore[0] = 0
                else:
                    plusScore[index]  = plusScore[index-1] + 1
                    minusScore[index] = minusScore[index-1] - 1
                if minusScore[index] < 0:
                    minusScore[index] = 0

            if pair.readStrand == '-':
                if index == 0:
                    plusScore[0]  = 0
                    minusScore[0] = 1
                else:
                    minusScore[index] = minusScore[index-1] + 1
                    plusScore[index]  = plusScore[index-1] - 1
                if plusScore[index] < 0:
                    plusScore[index] = 0

            index = index + 1

        leftHiScorePlus  = 0
        leftHiScoreMinus = 0
        leftHiIndexPlus  = 0
        leftHiIndexMinus = 0
        
#        print "forward result (going left to right across the peak):"
        for i in range(0,index):
            if plusScore[i] > leftHiScorePlus:
                leftHiScorePlus = plusScore[i]
                leftHiIndexPlus = i
            if minusScore[i] > leftHiScoreMinus:
                leftHiScoreMinus = minusScore[i]
                leftHiIndexMinus = i
            
#            print "%d %d %d" % (i, plusScore[i], minusScore[i])
#        print "leftHiScorePlus: %d" % leftHiScorePlus + " index: %d" % leftHiIndexPlus
#        print "leftHiScoreMinus: %d" % leftHiScoreMinus + " index: %d" % leftHiIndexMinus

        # reset and run in the other direction across the peak
        plusScore  = [0] * len(self.pairs)
        minusScore = [0] * len(self.pairs)
        index = len(self.pairs) - 1
        for pair in reversed(self.pairs):
            if pair.readStrand == '+':
                if index == len(self.pairs) - 1:
                    plusScore[index]  = 1
                    minusScore[index] = 0
                else:
                    plusScore[index]  = plusScore[index+1] + 1
                    minusScore[index] = minusScore[index+1] - 1
                if minusScore[index] < 0:
                    minusScore[index] = 0

            if pair.readStrand == '-':
                if index == len(self.pairs) - 1:
                    plusScore[index]  = 0
                    minusScore[index] = 1
                else:
                    minusScore[index] = minusScore[index+1] + 1
                    plusScore[index] = plusScore[index+1] - 1
                if plusScore[index] < 0:
                    plusScore[index] = 0

            index = index - 1

        rightHiScorePlus  = 0
        rightHiScoreMinus = 0
        rightHiIndexPlus  = 0
        rightHiIndexMinus = 0

#        print "reverse result (going right to left across the peak):"
        for i in range(0,len(self.pairs)):
            if plusScore[i] > rightHiScorePlus:
                rightHiScorePlus = plusScore[i]
                rightHiIndexPlus = i
            if minusScore[i] > rightHiScoreMinus:
                rightHiScoreMinus = minusScore[i]
                rightHiIndexMinus = i

#            print "%d %d %d" % (i, plusScore[i], minusScore[i])
#        print "rightHiScorePlus: %d" % rightHiScorePlus + " index: %d" % rightHiIndexPlus
#        print "rightHiScoreMinus: %d" % rightHiScoreMinus + " index: %d" % rightHiIndexMinus

        if (leftHiIndexPlus-rightHiIndexMinus < leftHiIndexMinus-rightHiIndexPlus):
            myBreakPoint = BreakPoint(leftHiIndexPlus,
                                      rightHiIndexMinus,
                                      self.pairs[leftHiIndexPlus].readPos,
                                      self.pairs[rightHiIndexMinus].readPos,
                                      '+', '-',self.npairs,self)
            return myBreakPoint
	else: # will need to filter out bad breakpoints elsewhere since they're all reported
            myBreakPoint = BreakPoint(leftHiIndexMinus,
                                      rightHiIndexPlus,
                                      self.pairs[leftHiIndexMinus].readPos,
                                      self.pairs[rightHiIndexPlus].readPos,
                                      '-','+',self.npairs,self)
            return myBreakPoint

    # return the number of reads that are not on the expected strand between indexStart
    # and indexEnd, inclusive.
    def countWrongStrand(self,indexStart,indexEnd,expectedStrand):
        wsCount = 0
        for index in range(indexStart,indexEnd+1):
            if (self.pairs[index].readStrand != expectedStrand):
                wsCount = wsCount + 1
        return wsCount

    def countWrongMateStrand(self,indexStart,indexEnd,expectedStrand):
        wsCount = 0
        for index in range(indexStart,indexEnd+1):
            if (self.pairs[index].relStrand() != expectedStrand):
                wsCount = wsCount + 1
        return wsCount

    # compute the mean distance between reads across a range indices. Note that since a range
    # statement range(indexStart,indexEnd+1) is used, indexStart and indexEnd will both be 
    # included in the mean. The distance for each index is the difference in position between
    # the read with that index and the read with index-1 (so it doesn't make sense to call this
    # with index 0)
    def meanDist(self,indexStart,indexEnd):
        distSum = 0
        for index in range(indexStart,indexEnd+1):
            distSum = distSum + (self.pairs[index].readPos - self.pairs[index-1].readPos)
        if (indexEnd + 1 - indexStart) < 1:
            return 0
        else:
            return distSum/(indexEnd + 1 - indexStart)

    # get mean position of reads within the insertion for a range of indices (inclusive)
    # needs the expected length of a full-length insertion to determine the relative position
    # inside the insertion from the genomic coords
    def meanMatePos(self,indexStart,indexEnd):
        posSum = 0
        for index in range(indexStart,indexEnd+1):
            posSum = posSum + self.pairs[index].relPos()
        if (indexEnd + 1 - indexStart) < 1:
            return 0
        else:
            return posSum/(indexEnd + 1 - indexStart)

    # as meanMatePos, but returns maximum position in element
    def maxMatePos(self,indexStart,indexEnd):
        maxPos = self.pairs[indexStart].relPos()
        for index in range(indexStart,indexEnd+1):
            myPos = self.pairs[index].relPos()
            if (myPos > maxPos):
                maxPos = myPos
        return maxPos

    # as meanMatePos, but returns minimum position in element
    def minMatePos(self,indexStart,indexEnd):
        minPos = self.pairs[indexStart].relPos()
        for index in range(indexStart,indexEnd+1):
            myPos = self.pairs[index].relPos()
            if (myPos < minPos):
                minPos = myPos
        return minPos

    # return true if this peak has anything in it. Good for sanity checks.
    def hasPairs(self):
        if len(self.pairs) < 1:
            return 0
        else:
            return 1

    # return chromosome name.
    def getChr(self):
        if (self.pairs[0]):
            return self.pairs[0].chrom
        else:
            return 0

    # returns a dictionary of the sources of reads for this peak and the number
    # of reads each source was responsible for
    def getSources(self):
        sourceDict = dict()
        for pair in self.pairs:
            if sourceDict.has_key(pair.source):
                sourceDict[pair.source] = sourceDict[pair.source] + 1
            else :
                sourceDict[pair.source] = 1

        return sourceDict

    # returns the element family to which the majority of mates are mapped
    def primaryMateElt(self):
        eltDict = {}
        for pair in self.pairs:
            if eltDict.has_key(pair.mateElt):
                eltDict[pair.mateElt] = eltDict[pair.mateElt] + 1
            else:
                eltDict[pair.mateElt] = 1
        maxEltCount = 0
        maxEltName = ''
        for elt,count in eltDict.iteritems():
            if count > maxEltCount:
                maxEltCount = count
                maxEltName  = elt
        return maxEltName

    # infers the element strand given the indices of the reads flanking the breakpoint
    def eltStrand(self, indexLeft, indexRight):
        leftStrand  = self.pairs[indexLeft].relStrand()
        rightStrand = self.pairs[indexRight].relStrand()
        if leftStrand == '-' and rightStrand == '+':
            return '+'
        if leftStrand == '+' and rightStrand == '-':
            return '-'
        if leftStrand == rightStrand:
            if leftStrand == '+': # inversion due to twin priming (Ostertag et al.)
                return 'i'
            if leftStrand == '-': # ??
                return 'x'

class BreakPoint:
    """Stores information about a breakpoint consisting of a 5' side and a 3' side"""
    def __init__(self,leftindex,rightindex,leftpos,rightpos,leftstrand,rightstrand,npairs,peak):
        self.leftindex   = leftindex
        self.rightindex  = rightindex
        self.leftpos     = leftpos
        self.rightpos    = rightpos
        self.leftstrand  = leftstrand
        self.rightstrand = rightstrand
        self.npairs      = npairs
        self.peak        = peak

    # functions to find the mean distance between reads on the left
    # and right sides of the breakpoint. Useful in validation.
    def leftMeanDist(self):
        return self.peak.meanDist(1,self.leftindex)
    def rightMeanDist(self):
        return self.peak.meanDist(self.rightindex+1,self.npairs-1)
    
    # useful for validation
    def leftWrongStrand(self):
        return self.peak.countWrongStrand(0,self.leftindex,self.leftstrand)

    def rightWrongStrand(self):
        return self.peak.countWrongStrand(self.rightindex,self.npairs-1,self.rightstrand)

    def leftMateWrongStrand(self):
        leftMateStrand = self.peak.pairs[self.leftindex].relStrand()
        return self.peak.countWrongMateStrand(0,self.leftindex,leftMateStrand)

    def rightMateWrongStrand(self):
        rightMateStrand = self.peak.pairs[self.rightindex].relStrand()
        return self.peak.countWrongMateStrand(self.rightindex,self.npairs-1,rightMateStrand)

    # leftmost/rightmost peak coords
    def minPeakPos(self):
        return self.peak.pairs[0].readPos
    def maxPeakPos(self):
        return self.peak.pairs[self.npairs-1].readPos

    # mean positions of reads in elements
    def leftMateMeanPos(self):
        return self.peak.meanMatePos(0,self.leftindex)
    def rightMateMeanPos(self):
        return self.peak.meanMatePos(self.rightindex,self.npairs-1)

    # min positions of reads in elements
    def leftMateMinPos(self):
        return self.peak.minMatePos(0,self.leftindex)
    def rightMateMinPos(self):
        return self.peak.minMatePos(self.rightindex,self.npairs-1)

    # max positions of reads in elements
    def leftMateMaxPos(self):
        return self.peak.maxMatePos(0,self.leftindex)
    def rightMateMaxPos(self):
        return self.peak.maxMatePos(self.rightindex,self.npairs-1)

    # element strand from read strands (can be 'i' if inverted or 'x' if FUBAR)
    def eltStrandFromStrand(self):
        return self.peak.eltStrand(self.leftindex,self.rightindex)

    # gets element strand from where the reads are mapped with in the element
    def eltStrandFromReads(self):
        minside='left'
        maxside='right'
        if (self.rightMateMinPos() < self.leftMateMinPos()):
            minside='right'
        if (self.leftMateMaxPos() > self.rightMateMaxPos()):
            maxside='left'

        if (minside == maxside):
            return 'x'
        else:
            if (minside=='left' and maxside=='right'):
                return '+'
            if (minside=='right' and maxside=='left'):
                return '-'

    def sourceList(self):
        return self.peak.getSources().keys()

    # uses the other two strand calling function to make a best guess as to the orientation
    # of the element, returns '.' by default
    def bestStrand(self):
        sFS = self.eltStrandFromStrand()
        sFR = self.eltStrandFromReads()
        if sFS == sFR and sFS != 'x':
            return sFS
        if (sFS == 'i' or sFS == 'x') and sFR != 'x':
            return sFR
        if sFR == 'x' and sFS != 'x' and sFS != 'i':
            return sFS
        if sFR != 'x' and sFS !='x' and sFS != 'i' and sFR != sFS:
            return sFS
        return '.'
        
    # extrema for positions within the element
    def eltMinPos(self):
        return min(self.leftMateMinPos(),self.rightMateMinPos())
    def eltMaxPos(self):
        return max(self.leftMateMaxPos(),self.rightMateMaxPos())
    # element strand from read positions, can try to resolve orientations of inverted elements
    def eltStrandFromPos(self):
        pass

    # returns whether a breakpoint meets the following conditions:
    # genomic 'anchor' alignments are either +/- or -/+ (5'/3')
    # reads are within one position of each other
    # peak size cutoff
    def isValid(self):
        if math.fabs(self.leftindex-self.rightindex) > 1:
            return False
        if self.leftstrand != '+' or self.rightstrand != '-':
            return False
        # require one read on each junction
        if self.leftindex < 1 or self.npairs-self.rightindex-1 < 1:
            return False
        # requirements for mean inter-read distances inside of a peak
        if self.leftMeanDist() == 0 or self.rightMeanDist() == 0 or self.leftMeanDist() > 100 or self.rightMeanDist() > 100:
            return False
        # less than 20% of reads on either side can have an inconsistent strand
        if ((float(self.leftWrongStrand()) / float(self.leftindex+1) > .2) or 
            (float(self.rightWrongStrand()) / float(self.rightindex-self.leftindex+1) > .2)):
            return False
        if float(self.leftWrongStrand() + self.rightWrongStrand())/float(self.npairs) > .4:
            return False
        # less than 10% of mates on either side can have an inconsistent strand
        if ((float(self.leftMateWrongStrand()) / float(self.leftindex+1) > .1) or
            (float(self.rightMateWrongStrand()) / float(self.rightindex-self.leftindex+1) > .1)):
            return False
        if float(self.leftMateWrongStrand() + self.rightMateWrongStrand())/float(self.npairs) > .1:
            return False
        return True



