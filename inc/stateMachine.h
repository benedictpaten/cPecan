/*
 * stateMachine.h
 *
 *  Created on: Aug 8, 2014
 *      Author: benedictpaten
 */

#ifndef STATEMACHINE_H_
#define STATEMACHINE_H_

#include "sonLib.h"

#define SYMBOL_NUMBER 5
#define SYMBOL_NUMBER_NO_N 4

/*
 * The statemachine object for computing pairwise alignments with.
 */

// TODO expand this?
typedef enum {
    fiveState=0,
    fiveStateAsymmetric=1,
    threeState=2,
    threeStateAsymmetric=3
} StateMachineType;

typedef struct _stateMachine StateMachine;
typedef struct _hmm Hmm;

/*
 * Hmm for loading/unloading HMMs and storing expectations.
 * Maybe move these definitions to stateMachine.c to clean this up?
 */

typedef struct _hmm {
    double likelihood;
    StateMachineType type;
    int64_t stateNumber;
    int64_t symbolSetSize;
    int64_t matrixSize;
    double *transitions;
    double *emissions;

    //void (*writeFcn)(struct hmm, FILE *fileHandle); TODO figure out how this is going to work

    void (*addToTransitionExpectationFcn)(Hmm *hmm, int64_t from, int64_t to, double p);

    void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p);

    double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to);

    void (*addToEmissionExpectationFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);

    void (*setEmissionExpectationFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);

    double (*getEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y);

    //void (*loadSymmetric)(StateMachine sM, Hmm hmm);
    //void (*loadAsymmetric)(struct stateMachine, struct hmm);
};

struct _stateMachine {
    StateMachineType type; // TODO get rid of this
    int64_t stateNumber;
    int64_t matchState;

    double (*startStateProb)(StateMachine *sM, int64_t state);

    double (*endStateProb)(StateMachine *sM, int64_t state);

    double (*raggedEndStateProb)(StateMachine *sM, int64_t state);

    double (*raggedStartStateProb)(StateMachine *sM, int64_t state);

    void (*setMatchDefaultsFcn)(double* emissionMatchProbs);

    void (*setXGapDefaultsFcn)(double* emissionXGapProbs);

    void (*setYGapDefaultsFcn)(double* emissionYGapProbs);

    double (*getXGapProbFcn)(const double* emissionGapProbs, void *i);

    double (*getYGapProbFcn)(const double* emissionGapProbs, void *i);

    double (*getMatchProbFcn)(const double* emissionMatchProbs, void *x, void *y);

    //void (*updateExpectationsFcn)(Hmm *hmmExpectations, int64_t from, int64_t to, double p);

    //Cells (states at a given coordinate)
    void (*cellCalculate)(StateMachine *sM, double *current, double *lower, double *middle, double *upper,
                          void* cX, void* cY,
                          void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
                          void *extraArgs);
};

typedef struct _stateMachineFunctions {
    double (*gapXProbFcn)(const double *, void *);
    double (*gapYProbFcn)(const double *, void *);
    double (*matchProbFcn)(const double *, void *, void *);
} StateMachineFunctions;


// constructers
Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type);
Hmm *hmm_kmer_constructEmpty(double pseudoExpectation, StateMachineType type);

// randomizers
void hmm_randomise(Hmm *hmm); //Creates normalised HMM with parameters set to small random values.
void hmm_kmer_randomise(Hmm *hmm);

// destruct
void hmm_destruct(Hmm *hmmExpectations);

// writers
void hmm_write(Hmm *hmmExpectations, FILE *fileHandle);
void hmm_kmer_write(Hmm *hmm, FILE *fileHandle);

// transitions
void hmm_addToTransitionExpectation(Hmm *hmmExpectations, int64_t from, int64_t to, double p);
double hmm_getTransition(Hmm *hmmExpectations, int64_t from, int64_t to);
void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p);

// emissions
void hmm_addToEmissionsExpectation(Hmm *hmmExpectations, int64_t state, int64_t x, int64_t y, double p);
void hmm_kmer_addToEmissionsExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p); //kmer addition

double hmm_getEmissionsExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y);
double hmm_kmer_getEmissionsExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y); //kmer addition

void hmm_setEmissionsExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);
void hmm_kmer_setEmissionsExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p);

// loaders
Hmm *hmm_loadFromFile(const char *fileName);
Hmm *hmm_kmer_loadFromFile(const char *fileName);

// normalizers
void hmm_normalise(Hmm *hmm);
void hmm_kmer_normalise(Hmm *hmm);

//////////////////
// stateMachine //
//////////////////

// indexing
int64_t emissions_getKmerIndex(void *kmer);

int64_t emissions_getBaseIndex(void *base);

// defaults
void emissions_symbol_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_symbol_setGapProbsToDefaults(double *emissionGapProbs);

void emissions_initMatchProbsToZero(double *emissionMatchProbs, int64_t symbolSetSize);

void emissions_initGapProbsToZero(double *emissionGapProbs, int64_t symbolSetSize);

// get prob
double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base);

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

double emissions_kmer_getGapProb(const double *emissionGapProbs, void *kmer);

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

// EM
StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs);

StateMachine *hmm_getStateMachine(Hmm *hmm);

StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs);

StateMachine *hmm_kmer_getStateMachine(Hmm *hmm);

StateMachine *stateMachine5_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setXGapDefaultsFcn)(double *),
                                      void (*setYGapDefaultsFcn)(double *),
                                      void (*setMatchDefaultsFcn)(double *),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *));

StateMachine *stateMachine3_construct(StateMachineType type); //the type is to specify symmetric/asymmetric

void stateMachine_destruct(StateMachine *stateMachine);

#endif /* STATEMACHINE_H_ */
