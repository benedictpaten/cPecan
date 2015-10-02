#!/usr/bin/env python

from __future__ import print_function
from itertools import product


class EmissionMatrix(object):
    """base class for code-generating emission matricies for cPecan"""
    def __init__(self, alpahbet, kmerLength):
        #super(emissionMatrix, self).__init__()
        self.alphabet = alpahbet
        self.kmerLength = kmerLength
        self.makeXaYKmers()

    def makeXaYKmers(self):
        self.Xkmers = []
        self.Ykmers = []

        for _ in product(self.alphabet, repeat=self.kmerLength):
            self.Xkmers.append(''.join(_))
            self.Ykmers.append(''.join(_))


    def matchKmers(self, kmerX, kmerY):
        prob = ''
        kmerX, kmerY = kmerX, kmerY
        for x, y in zip(kmerX, kmerY):
            prob += self.matchNucleotide(x, y)
        return prob

    def matchNucleotide(self, x, y):
        transitions = [
            ('A', 'G'), ('G', 'A'),
            ('C', 'T'), ('T', 'C')
        ]
        if (x == 'N') or (y == 'N'):
            return 'N' # for NN match not gap
        if (x == y) and (x != 'N'):
            return 'M'
        if (x, y) in transitions:
            return 'S'
        else:
            return 'V'

    def generateMatrix(self):
        probs = []

        Xkmers = self.Xkmers
        Ykmers = self.Ykmers

        for kmerX in Xkmers:
            lastX = kmerX[-1]
            for kmerY in Ykmers:
                prob = ''
                lastY = kmerY[-1]

                if lastX == lastY:
                    if lastX == 'N':
                        prob += 'N'
                    else:
                        prob += 'M'

                if (kmerX[:-1] == kmerY[:-1]) and (lastX != lastY):
                    prob = self.matchNucleotide(lastX, lastY)

                if (kmerX[:-1] != kmerY[:-1]) and (lastX != lastY):
                    prob = self.matchKmers(kmerX, kmerY)

                probs.append(prob)
                yield prob

    def generateMatrix2(self):
        probs = []
        Xkmers = self.Xkmers
        Ykmers = self.Ykmers

        for kmerX in Xkmers:
            for kmerY in Ykmers:
                prob = ''
                for n in range(len(kmerX)):
                    prob += self.matchNucleotide(kmerX[n], kmerY[n])
                yield prob
