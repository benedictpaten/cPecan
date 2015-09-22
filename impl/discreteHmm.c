#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "discreteHmm.h"

#include "../inc/emissionMatrix.h"
#include "../inc/stateMachine.h"
#include "../../sonLib/lib/sonLibCommon.h"
#include "../../sonLib/lib/sonLibList.h"
#include "../../sonLib/lib/sonLibString.h"
#include "../../sonLib/lib/sonLibFile.h"
#include "../inc/pairwiseAligner.h"
#include "../../sonLib/lib/sonLibRandom.h"
#include "../inc/discreteHmm.h"

// Construct
Hmm *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                        void (*addToTransitionExpFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to, double p),
                                        void (*setTransitionFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to, double p),
                                        double (*getTransitionsExpFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to),
                                        void (*addEmissionsExpFcn)(double *emissions, int64_t nSymbols, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        void (*setEmissionExpFcn)(double *emissions, int64_t nSymbols, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        double (*getEmissionExpFcn)(double *emissions, int64_t nSymbols, int64_t matrixSize, int64_t state, int64_t x, int64_t y)) {
    // malloc
    Hmm *hmmD = st_malloc(sizeof(Hmm));

    // Set up constants
    hmmD->stateNumber = stateNumber;
    hmmD->symbolSetSize = symbolSetSize;
    hmmD->matrixSize = symbolSetSize*symbolSetSize; // working with symmetric matrices
    hmmD->type = type;

    // Set up transitions matrix
    hmmD->transitions = st_malloc(hmmD->stateNumber * hmmD->stateNumber * sizeof(double));
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->stateNumber; i++) {
        hmmD->transitions[i] = pseudocount;
    }

    // Set up emissions matrix
    hmmD->emissions = st_malloc(hmmD->stateNumber * hmmD->matrixSize * sizeof(double));
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->matrixSize; i++) {
        hmmD->emissions[i] = pseudocount;
    }

    // Initialize likelihood
    hmmD->likelihood = 0.0;

    // Set up functions
    // transitions
    hmmD->addToTransitionExpectationFcn = addToTransitionExpFcn; //add
    hmmD->setTransitionFcn = setTransitionFcn; // set
    hmmD->getTransitionsExpFcn = getTransitionsExpFcn; // get
    // emissions
    hmmD->addToEmissionExpectationFcn = addEmissionsExpFcn; // add
    hmmD->setEmissionExpectationFcn = setEmissionExpFcn; // set
    hmmD->getEmissionExpFcn = getEmissionExpFcn; // get

    return hmmD;
}
// Transitions
void hmmDiscrete_addToTransitionExpectation(double *transitions, int64_t stateNumber,
                                            int64_t from, int64_t to, double p) {
    transitions[from * stateNumber + to] += p;
}

void hmmDiscrete_setTransitionExpectation(double *transitions, int64_t stateNumber, int64_t from, int64_t to, double p) {
    transitions[from * stateNumber + to] = p;
}

double hmmDiscrete_getTransitionExpectation(double *transitions, int64_t stateNumber, int64_t from, int64_t to) {
    return transitions[from * stateNumber + to];
}

// Emissions
void hmmDiscrete_addToEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                          int64_t state, int64_t x, int64_t y, double p) {
    int64_t tableIndex = x * symbolSetSize + y;
    emissions[(state * matrixSize) + tableIndex] += p;
}

void hmmDiscrete_setEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                        int64_t state, int64_t x, int64_t y, double p) {
    int64_t tableIndex = x * symbolSetSize + y;
    emissions[(state * matrixSize) + tableIndex] = p;
}

double hmmDiscrete_getEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                          int64_t state, int64_t x, int64_t y) {
    int64_t tableIndex = x * symbolSetSize + y;
    return emissions[(state * matrixSize) + tableIndex];
}

// Randomize/Normalize
void hmmDiscrete_randomize(Hmm *hmmD) {
    // Transitions
    for (int64_t from = 0; from < hmmD->stateNumber; from++) {
        for (int64_t to = 0; to < hmmD->stateNumber; to++) {
            hmmDiscrete_setTransitionExpectation(hmmD->transitions, hmmD->stateNumber, from, to, st_random());
        }
    }
    // Emissions
    for (int64_t state = 0; state < hmmD->stateNumber; state++) {
        for (int64_t x = 0; x < hmmD->symbolSetSize; x++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                hmmDiscrete_setEmissionExpectation(hmmD->emissions, hmmD->symbolSetSize, hmmD->matrixSize,
                                                   state, x, y, st_random());
            }
        }
    }
    hmmDiscrete_normalize(hmmD);
}

void hmmDiscrete_normalize(Hmm *hmmD) {
    // Transitions
    for (int64_t from = 0; from < hmmD->stateNumber; from++) {
        double total = 0.0;
        for (int64_t to = 0; to < hmmD->stateNumber; to++) {
            total += hmmDiscrete_getTransitionExpectation(hmmD->transitions, hmmD->stateNumber, from, to);
        }
        for (int64_t to = 0; to < hmmD->stateNumber; to++) {
            double newProb = hmmDiscrete_getTransitionExpectation(hmmD->transitions, hmmD->stateNumber,
                                                                  from, to) / total;
            hmmDiscrete_setTransitionExpectation(hmmD->transitions, hmmD->stateNumber, from, to, newProb);
        }
    }
    for (int64_t state = 0; state < hmmD->stateNumber; state++) {
        double total = 0.0;
        for (int64_t x = 0; x < hmmD->symbolSetSize; x++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                total += hmmDiscrete_getEmissionExpectation(hmmD->emissions, hmmD->symbolSetSize, hmmD->matrixSize,
                                                            state, x, y);
            }
        }
        for (int64_t x = 0; x < hmmD->symbolSetSize; x ++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                double newProb = hmmDiscrete_getEmissionExpectation(hmmD->emissions, hmmD->symbolSetSize, hmmD->matrixSize,
                                                                    state, x, y) / total;
                hmmDiscrete_setEmissionExpectation(hmmD->emissions, hmmD->symbolSetSize, hmmD->matrixSize,
                                                   state, x, y, newProb);
            }
        }
    }
}

// Loaders
//static void hmmDiscrete_loadSymmetric()

// Housekeeping
void hmmDiscrete_destruct(HmmDiscrete *hmmD) {
    free(hmmD->transitions);
    free(hmmD->emissions);
    free(hmmD);
}

// stateMachine interface
StateMachineFunctions *stateMachineFunctions_construct(double (*gapXProbFcn)(const double *, void *),
                                                       double (*gapYProbFcn)(const double *, void *),
                                                       double (*matchProbFcn)(const double *, void *, void *)) {
    StateMachineFunctions *sMfs = malloc(sizeof(StateMachineFunctions));
    sMfs->gapXProbFcn = gapXProbFcn;
    sMfs->gapYProbFcn = gapYProbFcn;
    sMfs->matchProbFcn = matchProbFcn;
    return sMfs;
}

//StateMachineFunctions *stateMachineFunctions_constructFromType(int64_t stateMachineType) {
//
//}