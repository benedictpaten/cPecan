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
#define MODEL_PARAMS 4 // we record the the mean and standard deviation for the level and the noise for each kmer
#define SQRT_TWO 1.4142135623730951
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

    int64_t (*getElementIndexFcn)(void *);

    //void (*loadSymmetric)(StateMachine sM, Hmm hmm);
    //void (*loadAsymmetric)(struct stateMachine, struct hmm);
};

struct _stateMachine {
    StateMachineType type; // TODO get rid of this
    int64_t stateNumber;
    int64_t matchState;
    int64_t parameterSetSize;
    double *EMISSION_MATCH_PROBS; //Match emission probs
    double *EMISSION_GAP_X_PROBS; //Gap emission probs
    double *EMISSION_GAP_Y_PROBS; //Gap emission probs

    double (*startStateProb)(StateMachine *sM, int64_t state);

    double (*endStateProb)(StateMachine *sM, int64_t state);

    double (*raggedEndStateProb)(StateMachine *sM, int64_t state);

    double (*raggedStartStateProb)(StateMachine *sM, int64_t state);

    //void (*setMatchDefaultsFcn)(double* emissionMatchProbs);
    //void (*setXGapDefaultsFcn)(double* emissionXGapProbs);
    //void (*setYGapDefaultsFcn)(double* emissionYGapProbs);

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

typedef struct _StateMachine5 StateMachine5;

struct _StateMachine5 {
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_SHORT_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_X; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_X; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_X; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_MATCH_FROM_SHORT_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_Y; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_Y; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_Y; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_Y; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_Y; //0.0073673675173412815f;
    //double *EMISSION_MATCH_PROBS; //Match emission probs
    //double *EMISSION_GAP_X_PROBS; //Gap emission probs
    //double *EMISSION_GAP_Y_PROBS; //Gap emission probs

};

typedef struct _StateMachine3 StateMachine3;

struct _StateMachine3 {
    //3 state state machine, allowing for symmetry in x and y.
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_GAP_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_SWITCH_TO_Y; //0.0073673675173412815f;
    //double *EMISSION_MATCH_PROBS; //Match emission probs
    //double *EMISSION_GAP_X_PROBS; //Gap X emission probs
    //double *EMISSION_GAP_Y_PROBS; //Gap Y emission probs
};

typedef struct _stateMachineFunctions {
    double (*gapXProbFcn)(const double *, void *);
    double (*gapYProbFcn)(const double *, void *);
    double (*matchProbFcn)(const double *, void *, void *);
} StateMachineFunctions;

//////////////////
// stateMachine //
//////////////////

// StateMachine constructors
StateMachine *stateMachine5_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setEmissionsDefaults)(StateMachine *sM),
        //void (*setXGapDefaultsFcn)(double *),
        //void (*setYGapDefaultsFcn)(double *),
        //void (*setMatchDefaultsFcn)(double *),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *));

StateMachine *stateMachine3_construct(StateMachineType type, int64_t parameterSetSize,
                                      void (*setEmissionsDefaults)(StateMachine *sM),
        //void (*setXGapDefaultsFcn)(double *),
        //void (*setYGapDefaultsFcn)(double *),
        //void (*setMatchDefaultsFcn)(double *),
                                      double (*gapXProbFcn)(const double *, void *),
                                      double (*gapYProbFcn)(const double *, void *),
                                      double (*matchProbFcn)(const double *, void *, void *));

// indexing
int64_t emissions_getKmerIndex(void *kmer);

int64_t emissions_getBaseIndex(void *base);

// defaults
void emissions_symbol_setMatchProbsToDefaults(double *emissionMatchProbs);

void emissions_symbol_setGapProbsToDefaults(double *emissionGapProbs);

void emissions_discrete_initEmissionsToZero(StateMachine *sM);

void emissions_symbol_setEmissionsToDefaults(StateMachine *sM);

void emissions_signal_initEmissionsToZero(StateMachine *sM);

// get prob
double emissions_symbol_getGapProb(const double *emissionGapProbs, void *base);

double emissions_symbol_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

double emissions_kmer_getGapProb(const double *emissionGapProbs, void *kmer);

double emissions_kmer_getMatchProb(const double *emissionMatchProbs, void *x, void *y);

double emissions_signal_getKmerGapProb(const double *kmerGapModel, void *kmer);

double emissions_signal_getEventGapProb(const double *eventGapModel, void *event);

double emissions_signal_getlogGaussPDFMatchProb(const double *eventModel, void *kmer, void *event);

double emissions_signal_getBivariateGaussPdfMatchProb(const double *eventModel, void *kmer, void *event);

// signal
//void emissions_signal_loadPoreModel(StateMachine *sM, const char *modelFile);

double emissions_signal_getModelEntry(const double *model, void *kmer);

void emissions_signal_scaleModel(StateMachine *sM, double scale, double shift, double var,
                                 double scale_sd, double var_sd);

StateMachine *getSignalStateMachine3(const char *modelFile, StateMachineFunctions *sMfs);

// EM
StateMachine *getStateMachine5(Hmm *hmmD, StateMachineFunctions *sMfs);

void stateMachine_destruct(StateMachine *stateMachine);

#endif /* STATEMACHINE_H_ */
