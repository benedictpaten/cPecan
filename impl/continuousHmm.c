#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "discreteHmm.h"
#include "stateMachine.h"
#include "continuousHmm.h"

// maybe make this static
static HmmContinuous *hmmContinuous_constructEmpty(
        int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
        void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
        void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
        double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
        void (*addEmissionsExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
        void (*setEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
        double (*getEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y),
        int64_t (*getElementIndexFcn)(void *)) {
    // malloc
    HmmContinuous *hmmC = st_malloc(sizeof(HmmContinuous));

    // setup base class
    hmmC->baseHmm.type = type;
    hmmC->baseHmm.stateNumber = stateNumber;
    hmmC->baseHmm.symbolSetSize = symbolSetSize;
    hmmC->baseHmm.matrixSize = MODEL_PARAMS;

    // initialize match models, for storage in between iterations
    hmmC->matchModel = st_malloc(hmmC->baseHmm.matrixSize * hmmC->baseHmm.symbolSetSize * sizeof(double));
    hmmC->extraEventMatchModel = st_malloc(hmmC->baseHmm.matrixSize * hmmC->baseHmm.symbolSetSize
                                           * sizeof(double));

    // Set up functions
    // transitions
    hmmC->baseHmm.addToTransitionExpectationFcn = addToTransitionExpFcn; // add
    hmmC->baseHmm.setTransitionFcn = setTransitionFcn;                   // set
    hmmC->baseHmm.getTransitionsExpFcn = getTransitionsExpFcn;           // get
    // emissions
    hmmC->baseHmm.addToEmissionExpectationFcn = addEmissionsExpFcn;      // add
    hmmC->baseHmm.setEmissionExpectationFcn = setEmissionExpFcn;         // set
    hmmC->baseHmm.getEmissionExpFcn = getEmissionExpFcn;                 // get
    // indexing
    hmmC->baseHmm.getElementIndexFcn = getElementIndexFcn;               // indexing

    return hmmC;
}

Hmm *continuousPairHmm_constructEmpty(
        double pseudocount, int64_t stateNumber,
        int64_t symbolSetSize, StateMachineType type,
        void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
        void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
        double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
        void (*addToKmerGapExpFcn)(Hmm *hmm, int64_t state, int64_t ki, int64_t ignore, double p),
        void (*setKmerGapExpFcn)(Hmm *hmm, int64_t state, int64_t ki, int64_t ignore, double p),
        double (*getKmerGapExpFcn)(Hmm *hmm, int64_t state, int64_t ki, int64_t ignore),
        int64_t (*getElementIndexFcn)(void *)) {
    // malloc
    ContinuousPairHmm *cpHmm = st_malloc(sizeof(ContinuousPairHmm));
    cpHmm->baseContinuousHmm =  *hmmContinuous_constructEmpty(stateNumber, symbolSetSize, type,
                                                              addToTransitionExpFcn,
                                                              setTransitionFcn,
                                                              getTransitionsExpFcn,
                                                              addToKmerGapExpFcn,
                                                              setKmerGapExpFcn,
                                                              getKmerGapExpFcn,
                                                              getElementIndexFcn);
    // transitions
    int64_t nb_states = cpHmm->baseContinuousHmm.baseHmm.stateNumber;
    cpHmm->transitions = st_malloc(nb_states * nb_states * sizeof(double));
    for (int64_t i = 0; i < (nb_states * nb_states); i++) {

        cpHmm->transitions[i] = pseudocount;
    }

    // individual kmer skip probs
    cpHmm->individualKmerGapProbs = st_malloc(cpHmm->baseContinuousHmm.baseHmm.symbolSetSize * sizeof(double));
    for (int64_t i = 0; i < cpHmm->baseContinuousHmm.baseHmm.symbolSetSize; i++) {
        cpHmm->individualKmerGapProbs[i] = pseudocount;
    }
    return (Hmm *) cpHmm;
}

// transitions
void continuousPairHmm_addToTransitionsExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    cpHmm->transitions[from * cpHmm->baseContinuousHmm.baseHmm.stateNumber + to] += p;
}

void continuousPairHmm_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    cpHmm->transitions[from * cpHmm->baseContinuousHmm.baseHmm.stateNumber + to] = p;
}

double continuousPairHmm_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    return cpHmm->transitions[from * cpHmm->baseContinuousHmm.baseHmm.stateNumber + to];
}

// kmer/gap emissions
void continuousPairHmm_addToKmerGapExpectation(Hmm *hmm, int64_t state, int64_t kmerIndex, int64_t ignore, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    (void) ignore;
    (void) state;
    cpHmm->individualKmerGapProbs[kmerIndex] += p;
}

void continuousPairHmm_setKmerGapExpectation(Hmm *hmm, int64_t state, int64_t kmerIndex, int64_t ignore, double p) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    (void) ignore;
    (void) state;
    // need a check for in-bounds kmer index?
    cpHmm->individualKmerGapProbs[kmerIndex] = p;
}

double continuousPairHmm_getKmerGapExpectation(Hmm *hmm, int64_t ignore, int64_t kmerIndex, int64_t ignore2) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    (void) ignore;
    (void) ignore2;
    return cpHmm->individualKmerGapProbs[kmerIndex];
}

// destructor
void continuousPairHmm_destruct(Hmm *hmm) {
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *) hmm;
    free(cpHmm->transitions);
    free(cpHmm->individualKmerGapProbs);
    free(cpHmm);
}

// normalizers/randomizers
void continuousPairHmm_normalize(Hmm *hmm) {
    // normalize transitions
    hmmDiscrete_normalize2(hmm, 0);
    // tally up the total
    double total = 0.0;
    for (int64_t i = 0; i < hmm->symbolSetSize; i++) {
        total += hmm->getEmissionExpFcn(hmm, 0, i, 0);
    }
    // normalize
    for (int64_t i = 0; i < hmm->symbolSetSize; i++) {
        double newProb = hmm->getEmissionExpFcn(hmm, 0, i, 0) / total;
        hmm->setEmissionExpectationFcn(hmm, 0, i, 0, newProb);
    }
}

void continuousPairHmm_randomize(Hmm *hmm) {
    // set all the transitions to random numbers
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm->setTransitionFcn(hmm, from, to, st_random());
        }
    }
    for (int64_t i = 0; i < hmm->symbolSetSize; i++) {
        hmm->setEmissionExpectationFcn(hmm, 0, i, 0, st_random());
    }
    continuousPairHmm_normalize(hmm);
}

void continuousPairHmm_loadTransitionsAndKmerGapProbs(StateMachine *sM, Hmm *hmm) {
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    // load transitions
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm->getTransitionsExpFcn(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, match));
    sM3->TRANSITION_MATCH_FROM_GAP_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_OPEN_X = log(hmm->getTransitionsExpFcn(hmm, match, shortGapX));
    sM3->TRANSITION_GAP_OPEN_Y = log(hmm->getTransitionsExpFcn(hmm, match, shortGapY));
    sM3->TRANSITION_GAP_EXTEND_X = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_EXTEND_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(hmm->getTransitionsExpFcn(hmm, shortGapY, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y = log(hmm->getTransitionsExpFcn(hmm, shortGapX, shortGapY));
    // load kmer gap probs
    for (int64_t i = 0; i < hmm->symbolSetSize; i++) {
        sM3->model.EMISSION_GAP_X_PROBS[i] = hmm->getEmissionExpFcn(hmm, 0, i, 0);
    }
}

Hmm *vanillaHmm_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                               void (*addToKmerBinExpFcn)(Hmm *hmm, int64_t bin, int64_t ignore, double p),
                               void (*setKmerBinFcn)(Hmm *hmm, int64_t bin, int64_t ignore, double p),
                               double (*getKmerBinExpFcn)(Hmm *hmm, int64_t bin, int64_t ignore),
                               int64_t (*getElementIndexFcn)(void *)) {
    VanillaHmm *vHmm = st_malloc(sizeof(VanillaHmm));

    vHmm->baseContinuousHmm =  *hmmContinuous_constructEmpty(stateNumber, symbolSetSize, type,
                                                             addToKmerBinExpFcn,
                                                             setKmerBinFcn,
                                                             getKmerBinExpFcn,
                                                             NULL,
                                                             NULL,
                                                             NULL,
                                                             getElementIndexFcn);
    vHmm->kmerSkipBins = st_malloc(30 * sizeof(double));
    for (int64_t i = 0; i < 30; i++) {
        vHmm->kmerSkipBins[i] = pseudocount;
    }

    vHmm->getKmerSkipBin = emissions_signal_getKmerSkipBin;

    return (Hmm *) vHmm;
}

void vanillaHmm_addToKmerSkipBinExpectation(Hmm *hmm, int64_t bin, int64_t ignore, double p) {
    VanillaHmm *vHmm = (VanillaHmm *) hmm;
    (void) ignore;
    vHmm->kmerSkipBins[bin] += p;
}

void vanillaHmm_setKmerSkipBinExpectation(Hmm *hmm, int64_t bin, int64_t ignore, double p) {
    VanillaHmm *vHmm = (VanillaHmm *) hmm;
    (void) ignore;
    vHmm->kmerSkipBins[bin] = p;
}

double vanillaHmm_getKmerSkipBinExpectation(Hmm *hmm, int64_t bin, int64_t ignore) {
    VanillaHmm *vHmm = (VanillaHmm *) hmm;
    (void) ignore;
    return vHmm->kmerSkipBins[bin];
}

void vanillaHmm_normalize(Hmm *hmm) {
    double total = 0.0;
    for (int64_t i = 0; i < 30; i++) {
        total += hmm->getTransitionsExpFcn(hmm, i, 0);
    }
    for (int64_t i = 0; i < 30; i++) {
        double newProb = hmm->getTransitionsExpFcn(hmm, i, 0) / total;
        hmm->setTransitionFcn(hmm, i, 0, newProb);
    }
}
