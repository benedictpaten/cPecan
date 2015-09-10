#ifndef EMISSIONS_MATRIX_H
#define EMISSIONS_MATRIX_H

#define KMER_LENGTH 2
#define NUM_OF_KMERS 25 // for alphabet 'AGCTN', may not need 'N' in there
#define MATRIX_SIZE 625 // assuming 25 kmers, 25*25
                        // may also want to make this file code-generated.

void emissions_kmer_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_kmer_setGapProbsToDefaults(double *emissionGapProbs);

#endif
