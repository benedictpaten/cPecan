#include <stdio.h>
#include <stdint.h>
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

Hmm *continuousPairHmm_constructEmpty(double pseudocount, int64_t stateNumber,
                                      int64_t symbolSetSize, StateMachineType type,
                                      void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                      void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                      double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
                                      int64_t (*getElementIndexFcn)(void *)) {
    // malloc
    ContinuousPairHmm *cpHmm = st_malloc(sizeof(ContinuousPairHmm));

    cpHmm->baseContinuousHmm =  *hmmContinuous_constructEmpty(stateNumber, symbolSetSize, type,
                                                              addToTransitionExpFcn,
                                                              setTransitionFcn,
                                                              getTransitionsExpFcn,
                                                              NULL,
                                                              NULL,
                                                              NULL,
                                                              getElementIndexFcn);
    // transitions
    int64_t nb_states = cpHmm->baseContinuousHmm.baseHmm.stateNumber;
    cpHmm->transitions = st_malloc(nb_states * nb_states * sizeof(double));
    for (int64_t i = 0; i < (nb_states * nb_states); i++) {
        cpHmm->transitions[i] = pseudocount;
    }

    // individual kmer skip probs
    for (int64_t i = 0; i < cpHmm->baseContinuousHmm.baseHmm.symbolSetSize; i++) {
        cpHmm->individualKmerSkipProbs[i] = pseudocount;
    }

    return (Hmm *) cpHmm;
}

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

Hmm *vanillaHmm_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                               void (*addToKmerBinExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                               void (*setKmerBinFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                               double (*getKmerBinExpFcn)(Hmm *hmm, int64_t from, int64_t to),
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
    (void) ignore; //TODO try this out
    vHmm->kmerSkipBins[bin] += p;
}

void vanillaHmm_setKmerSkipBinExpectation(Hmm *hmm, int64_t bin, int64_t ignore, double p) {
    VanillaHmm *vHmm = (VanillaHmm *) hmm;
    (void) ignore; //TODO try this out
    vHmm->kmerSkipBins[bin] = p;
}

double vanillaHmm_getKmerSkipBinExpectation(Hmm *hmm, int64_t bin, int64_t ignore) {
    VanillaHmm *vHmm = (VanillaHmm *) hmm;
    (void) ignore; //TODO try this out
    return vHmm->kmerSkipBins[bin];
}

