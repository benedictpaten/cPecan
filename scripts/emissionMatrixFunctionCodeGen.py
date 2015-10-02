#!/usr/bin/env python
from __future__ import print_function
from argparse import ArgumentParser
import sys
import os
import emissionMatrixLib as eMLib

def parse_args():
    """Parses command line arguments.
    """
    parser = ArgumentParser (description=__doc__)
    parser.add_argument ('--outFile', '-o', action='store', dest='outFile',
                         default=sys.stdout, required=False, type=str,
                         help='''output .c file, does not make header file'''
                        )
    parser.add_argument ('--alphabet', '-a', action='store', dest='alpha',
                         default='ATGCN', required=False, type=str,
                         help='''alphabet to use'''
                        )
    parser.add_argument ('--kmerLength', '-l', action='store', dest='kLength',
                         default=2, required=False, type=int,
                         help='''kmer length to use'''
                        )
    args = parser.parse_args()
    return args

header = """#include "{}.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>"""

emissionsMatchfunction_start = """
void emissionsKmers_setMatchProbsToDefaults(double *emissionMatchProbs) {
    /*
     * This sets the match probabilities to default values for matching kmers
     */

    const double M=-2.1149196655034745; //log(0.12064298095701059);
    const double V=-4.5691014376830479; //log(0.010367271172731285);
    const double S=-3.9833860032220842; //log(0.01862247669752685);
    const double N=-2.772588722;        //log(0.25**2)

    //Symmetric matrix of emission probabilities.

    const double i[MATRIX_SIZE] = {
        """
emissionsMatchfunction_end = """};

    memcpy(emissionMatchProbs, i, sizeof(double)*MATRIX_SIZE);
}
\n"""

emissionsGapFunction_start = """
void emissionsKmers_setGapProbsToDefaults(double *emissionGapProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double G = -1.6094379124341003; //log(0.2)
    const double i[%d] = {
        """
emissionsGapFunction_end = """};

    memcpy(emissionGapProbs, i, sizeof(double)*%d);
}
\n"""



def writeLines(generator, numOfCols, outFile):
    outFile = outFile
    lineBreak = 0
    for entry in generator:
        if len(entry) == 1:
            print(entry, end=', ', file=outFile)
        if len(entry) > 1:
            entry = '+'.join(list(entry))
            print(entry, end=', ', file=outFile)

        lineBreak += 1
        if lineBreak == numOfCols:
            print('', file=outFile)
            print('        ', end='', file=outFile)
            lineBreak = 0


def main():

    # Do the set up
    args = parse_args()
    eMatrix = eMLib.EmissionMatrix(args.alpha, args.kLength)
    num_of_kmers = len(eMatrix.Xkmers)

    # hande output
    if args.outFile != sys.stdout:
        outFile = open(args.outFile+'.c', 'a')
    if args.outFile == sys.stdout:
        outFile = sys.stdout

    # write out header and emission match default function
    print(header.format(args.outFile), file=outFile) # header
    print(emissionsMatchfunction_start, end='', file=outFile) # everything up
                                                              # to the matrix
    #writeLines(eMatrix.generateMatrix(), 25, outFile) # the matrix
    writeLines(eMatrix.generateMatrix2(), 25, outFile) # the matrix
    print(emissionsMatchfunction_end, file=outFile) # the end of the function
    print('\n', file=outFile) # to seperate functions

    print(emissionsGapFunction_start % num_of_kmers, end='',
         file=outFile) # start of gap emissions function

    # print out the gap emissions
    writeLines(['G' for k in eMatrix.Xkmers], 5, outFile)

    print(emissionsGapFunction_end % num_of_kmers,
          file=outFile) # end of the function


if __name__ == '__main__':
    main()
