#ifndef EMISSIONS_MATRIX_H
#define EMISSIONS_MATRIX_H

#define MATRIX_SIZE 625
#define NUM_OF_KMERS 25

void emissionsKmers_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissionsKmers_setGapProbsToDefaults(double *emissionGapProbs);

#endif
