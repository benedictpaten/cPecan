#ifndef DISCRETE_HMM_H_
#define DISCRETE_HMM_H_

#include <stdio.h>
#include <stdint.h>
#include "stateMachine.h"

typedef struct _hmmDiscrete {
    double likelihood;
    int64_t stateNumber;
    int64_t symbolSetSize;
    int64_t matrixSize
    double *transitions;
    double *emissions;

    void (*randomizeFcn)(double *transitions, double *emissions);
    void (*normalizeFcn)(double *transitions, double *emissions);

    //void (*writeFcn)(struct hmm, FILE *fileHandle); TODO figure out how this is going to work

    void (*addToTransitionExpectationFcn)(double *transitions, int64_t stateNumber,
                                       int64_t from, int64_t to, double p);

    void (*setTransitionFcn)(double *transitions, int64_t from, int64_t to, double p);

    double (*getTransitionsExpFcn)(double *transitions, int64_t stateNumber, int64_t from, int64_t to);

    void (*addToEmissionExpectationFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                        int64_t state, int64_t x, int64_t y, double p);

    void (*setEmissionExpectationFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                      int64_t state, int64_t x, int64_t y, double p);

    double (*getEmissionExpFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                int64_t state, int64_t x, int64_t y);

    //void (*loadSymmetric)(struct stateMachine, struct hmm);
    //void (*loadAsymmetric)(struct stateMachine, struct hmm);

} HmmDiscrete;

// Construct
HmmDiscrete *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize,
                                        void (*addToTransitionExpFcn)(double *transitions, int64_t stateNumber, int64_t from, int64_t to, double p),
                                        void (*setTransitionFcn)(double *transitions, int64_t stateNumber, int64_t from, int64_t to, double p),
                                        double (*getTransitionsExpFcn)(double *transitions, int64_t stateNumber, int64_t from, int64_t to),
                                        void (*addEmissionsExpFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        void (*setEmissionExpFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        double (*getEmissionExpFcn)(double *emissions, int64_t symbolSetSize, int64_t matrixSize, int64_t state, int64_t x, int64_t y));

// Transitions
void hmmDiscrete_addToTransitionExpectation(double *transitions, int64_t stateNumber,
                                            int64_t from, int64_t to, double p);
void hmmDiscrete_setTransitionExpectation(double *transitions, int64_t stateNumber,
                                          int64_t from, int64_t to, double p);
double hmmDiscrete_getTransitionExpectation(double *transitions, int64_t stateNumber, int64_t from, int64_t to);

// Emissions
void hmmDiscrete_addToEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                          int64_t state, int64_t x, int64_t y, double p);
void hmmDiscrete_setEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                        int64_t state, int64_t x, int64_t y, double p);
double hmmDiscrete_getEmissionExpectation(double *emissions, int64_t symbolSetSize, int64_t matrixSize,
                                          int64_t state, int64_t x, int64_t y);

// Randomize/Normalize
void hmmDiscrete_randomize(HmmDiscrete *hmmD);
void hmmDiscrete_normalize(HmmDiscrete *hmmD);


#endif