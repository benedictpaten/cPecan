#ifndef DISCRETE_HMM_H_
#define DISCRETE_HMM_H_

#include <stdio.h>
#include <stdint.h>
#include "stateMachine.h"

// Construct
Hmm *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                        void (*addToTransitionExpFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to, double p),
                                        void (*setTransitionFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to, double p),
                                        double (*getTransitionsExpFcn)(double *transitions, int64_t nStates, int64_t from, int64_t to),
                                        void (*addEmissionsExpFcn)(double *emissions, int64_t nStates, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        void (*setEmissionExpFcn)(double *emissions, int64_t nSymbols, int64_t matrixSize, int64_t state, int64_t x, int64_t y, double p),
                                        double (*getEmissionExpFcn)(double *emissions, int64_t nSymbols, int64_t matrixSize, int64_t state, int64_t x, int64_t y));

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
void hmmDiscrete_randomize(Hmm *hmmD);

void hmmDiscrete_normalize(Hmm *hmmD);

// Housekeeping
void hmmDiscrete_destruct(Hmm *hmmD);

// stateMachine interface
StateMachineFunctions *stateMachineFunctions_construct(double (*gapXProbFcn)(const double *, void *),
                                                       double (*gapYProbFcn)(const double *, void *),
                                                       double (*matchProbFcn)(const double *, void *, void *));

// Not yet implemented
//StateMachineFunctions *stateMachineFunctions_constructFromType(int64_t stateMachineType);

#endif
