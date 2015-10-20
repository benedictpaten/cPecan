#!/usr/bin/env python

from __future__ import division, print_function
from itertools import product

# Functions
    # A = 0
    # C = 1
    # G = 2
    # T = 3
    # N = 4

def nucleotideToIndex(base):
    if base == 'A':
        return 0
    if base == 'C':
        return 1
    if base == 'G':
        return 2
    if base == 'T':
        return 3
    if base == 'N':
        return 4

def getKmerIndex(kmer):
    """This is the algorithm for finding the rank (index) of a kmer)
    """
    axisLength = len(alphabet)**len(kmer)
    l = axisLength/len(alphabet)
    i = 0
    index = 0
    while l > 1:
        index += l*nucleotideToIndex(kmer[i])
        i += 1
        l = l/len(alphabet)
    index += nucleotideToIndex(kmer[-1])
    return int(index)

# Main program
kmers = []
alphabet = 'ACGTN'

for kmer in product(alphabet, repeat=6): #change repeat for longer kmers
    kmers.append(''.join(kmer))

#print("[", end="")
for _ in kmers:
    request = _
    x = getKmerIndex(request)
    k = kmers[x]
#    print("\"",_,"\"",",",sep='',end=" ")
#    print("\"",x,"\"",",",sep='',end=" ")
    print("testing:", x, kmers[x])
#    print(kmers[x])
#    if request == k:
#        print('All good')
    assert request == k
#print("]", end="\n")
