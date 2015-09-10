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

typedef enum {
    fiveState=0,
    fiveStateAsymmetric=1,
    threeState=2,
    threeStateAsymmetric=3
} StateMachineType;

typedef struct _stateMachine StateMachine;

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

    void (*updateExpectationsFcn)(double *fromCells, double *toCells, int64_t from, int64_t to,
                                  double eP, double tP, void *extraArgs);

    //Cells (states at a given coordinate)
    void (*cellCalculate)(StateMachine *sM, double *current, double *lower, double *middle, double *upper,
                          void* cX, void* cY,
                          void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
                          void *extraArgs);
};

/*
 * Hmm for loading/unloading HMMs and storing expectations.
 */

typedef struct _stateMachineFunctions {
    void (*setMatchDefaultsFcn)(double* emissionMatchProbs);
    void (*setXGapDefaultsFcn)(double* emissionXGapProbs);
    void (*setYGapDefaultsFcn)(double* emissionYGapProbs);
    double (*getXGapProbFcn)(const double* emissionGapProbs, void *i);
    double (*getYGapProbFcn)(const double* emissionGapProbs, void *i);
    double (*getMatchProbFcn)(const double* emissionMatchProbs, void *x, void *y);
    void (*updateExpectationsFcn)(double *fromCells, double *toCells, int64_t from, int64_t to,
                                  double eP, double tP, void *extraArgs);
} StateMachineFunctions;

typedef struct _hmm {
    StateMachineType type;
    double likelihood;
    int64_t stateNumber;
    StateMachineFunctions sM_functions;

    void (*randomizeFcn)(struct hmm);
    void (*normalizeFcn)(struct hmm);

    void (*writeFcn)(struct hmm, FILE *fileHandle);

    void (*addToTransitionExpectation)(struct hmm, int64_t from, int64_t to, double p);
    double (*getTransitions)(struct hmmExpectations, int64_t from, int64_t to);
    void (*setTransition)(struct hmm, int64_t from, int64_t to, double p);

    void (*addToEmissionExpectation);
    double (*getEmissionExpectation);
    void (*setEmissionExpectation);

    void (*loadSymmetric)(struct stateMachine, struct hmm);
    void (*loadAsymmetric)(struct stateMachine, struct hmm);

} Hmm;

typedef struct _hmm_discrete {
    Hmm baseHMM;
    int64_t matrixSize;
    double *transitions;
    double *emissions;
} Hmm_discrete;


// constructers
Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type,
                        void (*setEmissionsToZero)(Hmm*, int64_t, double*, double));

void hmm_symbol_setEmissionsToZero(Hmm *hmm, int64_t matrixSize, double* emissions, double pseudoExpectation);

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
int64_t emissions_getKmerIndex(void *kmer);

int64_t emissions_getBaseIndex(void *base);

void emissions_symbol_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_symbol_setGapProbsToDefaults(double *emissionGapProbs);

double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base);

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

double emissions_kmer_getGapProb(const double *emissionGapProbs, void *kmer);

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

StateMachine *hmm_getStateMachine(Hmm *hmm);

StateMachine *hmm_kmer_getStateMachine(Hmm *hmm);

StateMachine *stateMachine5_construct(StateMachineType type,
                                      void (*setMatchDefaultsFcn)(double *),
                                      void (*setXGapDefaultsFcn)(double *),
                                      void (*setYGapDefaultsFcn)(double *),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *),
                                      void (*updateExpectations)(double*, double*, int64_t, int64_t, double, double, void*));

StateMachine *stateMachine3_construct(StateMachineType type); //the type is to specify symmetric/asymmetric

void stateMachine_destruct(StateMachine *stateMachine);

#endif /* STATEMACHINE_H_ */
