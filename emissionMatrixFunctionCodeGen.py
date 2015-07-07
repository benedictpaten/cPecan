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

Cfunction_start = """
void emissions_setMatchProbsToDefaults(double *emissionMatchProbs) {
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
Cfunction_end = """};

    memcpy(emissionMatchProbs, i, sizeof(double)*MATRIX_SIZE);
}
\n"""


def main():
    args = parse_args()

    matrix = eMLib.emissionMatrix(args.alpha, args.kLength)

    if args.outFile != sys.stdout:
        outFile = open(args.outFile+'.c', 'a')
    if args.outFile == sys.stdout:
        outFile = sys.stdout

    print(header.format(args.outFile), file=outFile)
    print(Cfunction_start, end='', file=outFile)
    lineBreak = 0

    for prob in matrix.generateMatrix():
        if len(prob) == 1:
            print(prob, end=', ', file=outFile)
        if len(prob) > 1:
            prob = '*'.join(list(prob))
            print(prob, end=', ', file=outFile)

        lineBreak += 1
        if lineBreak == 10:
            print('', file=outFile)
            print('        ', end='', file=outFile)
            lineBreak = 0

    print(Cfunction_end, file=outFile)
    #print('\n', file=outFile)

if __name__ == '__main__':
    main()
