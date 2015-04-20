#!/usr/bin/env python
from cPecan.cPecanEm import Hmm, SYMBOL_NUMBER
import sys
import numpy as np
from optparse import OptionParser

"""Modifies an HMM model trained with margin align.
"""

toMatrix = lambda e : map(lambda i : e[SYMBOL_NUMBER*i:SYMBOL_NUMBER*(i+1)], xrange(SYMBOL_NUMBER))
fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)
    
def normaliseHmmByReferenceGCContent(hmm, gcContent):
    #Normalise background emission frequencies, if requested to GC% given
    for state in range(hmm.stateNumber):
        if state not in (2, 4): #Don't normalise GC content of insert states (as they don't have any ref bases!)
            n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)])
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), range(SYMBOL_NUMBER))) #Normalise

def modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate):
    #Normalise background emission frequencies, if requested to GC% given
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), xrange(SYMBOL_NUMBER**2)))
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

def setHmmIndelEmissionsToBeFlat(hmm):
    #Set indel emissions to all be flat
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2
    

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputModel outputModel [options]", 
                          version="%prog 0.1")

    parser.add_option("--substitutionRate", dest="substitutionRate", 
                      help="The probability per base of a difference between \
                      the sequenced reference and the reference the reads are aligned to. \
                      Value must be between 0 and 1.", 
                      default=0.00, type=float)
    
    parser.add_option("--gcContent", dest="gcContent", 
                      help="The desired GC content of the model. \
                      By default no adjustment is made; value must be between 0 and 1.", 
                      default=None, type=float)
    
    parser.add_option("--setFlatIndelEmissions", dest="setFlatIndelEmissions", 
                      help="Set all indel emissions probability to be equal regardless of base.", 
                      default=False, action="store_true")

    #Parse the options/arguments
    options, args = parser.parse_args()
    
    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    #Exit if the arguments are not what we expect
    if len(args) != 2:
        raise RuntimeError("Expected two arguments, got: %s" % " ".join(args))

    #Load HMM
    hmm = Hmm.loadHmm(sys.argv[1])

    #Normalise background emission frequencies, if requested to GC% given
    if options.gcContent != None:
        if options.gcContent < 0 or options.gcContent > 1.0:
            raise RuntimeError("Substitution rate is not a value between 0 and 1, got: %s" % options.gcContent)
        normaliseHmmByReferenceGCContent(hmm, options.gcContent)
    
    #Modify match emissions by proposed substitution rate
    if options.substitutionRate < 0 or options.substitutionRate > 1.0:
        raise RuntimeError("Substitution rate is not a value between 0 and 1, got: %s" % options.substitutionRate)
    modifyHmmEmissionsByExpectedVariationRate(hmm, options.substitutionRate)
    
    if options.setFlatIndelEmissions:
        setHmmIndelEmissionsToBeFlat(hmm)
    
    #Write out HMM
    hmm.write(sys.argv[2])

if __name__ == '__main__':
    main()

