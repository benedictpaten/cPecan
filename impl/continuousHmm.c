#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "discreteHmm.h"
#include "emissionMatrix.h"
#include "stateMachine.h"
#include "pairwiseAligner.h"
#include "continuousHmm.h"


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
    hmmC->baseHmm.likelihood = 0.0;

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
///////////////////////////////////////////// Continuous Pair HMM /////////////////////////////////////////////////////
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
    if (type != threeState) {
        st_errAbort("ContinuousPair HMM construct: Wrong HMM type for this function got: %i", type);
    }
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
    //sM3->TRANSITION_GAP_SWITCH_TO_X = LOG_ZERO;
    //sM3->TRANSITION_GAP_SWITCH_TO_Y = LOG_ZERO;

    /// / load kmer gap probs
    //for (int64_t i = 0; i < hmm->symbolSetSize; i++) {
        //sM3->model.EMISSION_GAP_X_PROBS[i] = hmm->getEmissionExpFcn(hmm, 0, i, 0);
    //    sM3->model.EMISSION_GAP_X_PROBS[i] = log(hmm->getEmissionExpFcn(hmm, 0, i, 0));
    //}
}

void continuousPairHmm_writeToFile(Hmm *hmm, FILE *fileHandle) {
    /*
     * Format:
     * type \t stateNumber \t symbolSetSize \n
     * [transitions... \t] likelihood \n
     * [kmer skip probs ... \t] \n
     */
    // downcast
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    // write the basic stuff to disk (positions are line:item#, not line:col)
    fprintf(fileHandle, "%i\t", cpHmm->baseContinuousHmm.baseHmm.type); // type 0:0
    fprintf(fileHandle, "%lld\t", cpHmm->baseContinuousHmm.baseHmm.stateNumber); // stateNumber 0:1
    fprintf(fileHandle, "%lld\t", cpHmm->baseContinuousHmm.baseHmm.symbolSetSize); // symbolSetSize 0:2
    fprintf(fileHandle, "\n"); // newLine

    // write the transitions to disk
    int64_t nb_transitions = (cpHmm->baseContinuousHmm.baseHmm.stateNumber
                              * cpHmm->baseContinuousHmm.baseHmm.stateNumber);
    fprintf(stdout, "exps:");
    for (int64_t i = 0; i < nb_transitions; i++) {
        fprintf(stdout, "%f ", cpHmm->transitions[i]);
        fprintf(fileHandle, "%f\t", cpHmm->transitions[i]); // transitions 1:(0-9)
    }
    fprintf(stdout, "\n");

    // write the likelihood
    fprintf(fileHandle, "%f\n", cpHmm->baseContinuousHmm.baseHmm.likelihood); // likelihood 1:10, newLine

    // write the individual kmer skip probs to disk
    for (int64_t i = 0; i < cpHmm->baseContinuousHmm.baseHmm.symbolSetSize; i++) {
        fprintf(fileHandle, "%f\t", cpHmm->individualKmerGapProbs[i]); // indiv kmer skip probs 2:(0-4096)
    }
    fprintf(fileHandle, "\n"); // newLine
}

Hmm *continuousPairHmm_loadFromFile(const char *fileName) {
    // open file
    FILE *fH = fopen(fileName, "r");

    // line 0
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    int type;
    int64_t stateNumber, symbolSetSize;
    int64_t j = sscanf(stList_get(tokens, 0), "%i", &type); // type
    if (j != 1) {
        st_errAbort("Failed to parse type (int) from string: %s\n", string);
    }
    int64_t s = sscanf(stList_get(tokens, 1), "%lld", &stateNumber); // stateNumber
    if (s != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    int64_t n = sscanf(stList_get(tokens, 2), "%lld", &symbolSetSize); // symbolSetSize
    if (n != 1) {
        st_errAbort("Failed to parse symbol set size (int) from string: %s\n", string);
    }

    // make empty cpHMM
    Hmm *hmm = continuousPairHmm_constructEmpty(0.0, stateNumber, symbolSetSize, type,
                                                  continuousPairHmm_addToTransitionsExpectation,
                                                  continuousPairHmm_setTransitionExpectation,
                                                  continuousPairHmm_getTransitionExpectation,
                                                  continuousPairHmm_addToKmerGapExpectation,
                                                  continuousPairHmm_setKmerGapExpectation,
                                                  continuousPairHmm_getKmerGapExpectation,
                                                  emissions_discrete_getKmerIndexFromKmer);

    // Downcast
    ContinuousPairHmm *cpHmm = (ContinuousPairHmm *)hmm;

    // cleanup
    free(string);
    stList_destruct(tokens);

    // Transitions
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    int64_t nb_transitions = (cpHmm->baseContinuousHmm.baseHmm.stateNumber
                              * cpHmm->baseContinuousHmm.baseHmm.stateNumber);

    // check for the correct number of transitions
    if (stList_length(tokens) != nb_transitions + 1) { // + 1 bc. likelihood is also on that line
        st_errAbort(
                "Incorrect number of transitions in the input HMM file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), nb_transitions + 1);
    }
    // load them
    for (int64_t i = 0; i < nb_transitions; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(cpHmm->transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(cpHmm->baseContinuousHmm.baseHmm.likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }
    // Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    // Emissions
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check
    if (stList_length(tokens) != cpHmm->baseContinuousHmm.baseHmm.symbolSetSize) {
        st_errAbort(
                "Incorrect number of emissions in the input HMM file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), cpHmm->baseContinuousHmm.baseHmm.symbolSetSize);
    }
    // load them
    for (int64_t i = 0; i < cpHmm->baseContinuousHmm.baseHmm.symbolSetSize; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(cpHmm->individualKmerGapProbs[i]));
        if (j != 1) {
            st_errAbort("Failed to parse the individual kmer skip probs from string %s\n", string);
        }
    }

    // Cleanup emissions line
    free(string);
    stList_destruct(tokens);
    // close file
    fclose(fH);

    return (Hmm *)cpHmm;

}

////////////////////////////////////////////////// Vanilla HMM ///////////////////////////////////////////////////////
Hmm *vanillaHmm_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                               void (*addToKmerBinExpFcn)(Hmm *hmm, int64_t bin, int64_t ignore, double p),
                               void (*setKmerBinFcn)(Hmm *hmm, int64_t bin, int64_t ignore, double p),
                               double (*getKmerBinExpFcn)(Hmm *hmm, int64_t bin, int64_t ignore)) {
    if (type != vanilla) {
        st_errAbort("Vanilla HMM construct: Wrong HMM type for this function got: %i", type);
    }

    VanillaHmm *vHmm = st_malloc(sizeof(VanillaHmm));

    vHmm->baseContinuousHmm =  *hmmContinuous_constructEmpty(stateNumber, symbolSetSize, type,
                                                             addToKmerBinExpFcn,
                                                             setKmerBinFcn,
                                                             getKmerBinExpFcn,
                                                             NULL,
                                                             NULL,
                                                             NULL,
                                                             NULL);
    vHmm->kmerSkipBins = st_malloc(60 * sizeof(double));
    for (int64_t i = 0; i < 60; i++) {
        vHmm->kmerSkipBins[i] = pseudocount;
    }
    // make a variable to make life easier, add 1 because of the correlation param
    int64_t nb_matchModelBuckets = 1 + (vHmm->baseContinuousHmm.baseHmm.symbolSetSize * MODEL_PARAMS);
    vHmm->matchModel = st_malloc(nb_matchModelBuckets * sizeof(double));
    vHmm->scaledMatchModel = st_malloc(nb_matchModelBuckets * sizeof(double));
    vHmm->getKmerSkipBin = emissions_signal_getKmerSkipBin;

    return (Hmm *) vHmm;
}
// transitions (kmer skip bins)
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

// normalize/randomize
void vanillaHmm_normalizeKmerSkipBins(Hmm *hmm) {
    double total = 0.0;
    for (int64_t i = 0; i < 60; i++) { // this is wrong, want to normalize alpha and beta seperately
        total += hmm->getTransitionsExpFcn(hmm, i, 0);
    }
    for (int64_t i = 0; i < 60; i++) {
        double newProb = hmm->getTransitionsExpFcn(hmm, i, 0) / total;
        hmm->setTransitionFcn(hmm, i, 0, newProb);
    }
}

void vanillaHmm_randomizeKmerSkipBins(Hmm *hmm) {
    for (int64_t i = 0; i < 60; i++) {
        hmm->setTransitionFcn(hmm, i, 0, st_random());
    }
    vanillaHmm_normalizeKmerSkipBins(hmm);
}

// load pore model
void vanillaHmm_implantMatchModelsintoHmm(StateMachine *sM, Hmm *hmm) {
    // down cast
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *)sM;
    VanillaHmm *vHmm = (VanillaHmm *)hmm;

    // go through the match and scaled match models and load them into the hmm
    int64_t nb_matchModelBuckets = 1 + (sM3v->model.parameterSetSize * MODEL_PARAMS);
    for (int64_t i = 0; i < nb_matchModelBuckets; i++) {
        vHmm->matchModel[i] = sM3v->model.EMISSION_MATCH_PROBS[i];
        vHmm->scaledMatchModel[i] = sM3v->model.EMISSION_GAP_Y_PROBS[i];
    }
}

// load into stateMachine
void vanillaHmm_loadKmerSkipBinExpectations(StateMachine *sM, Hmm *hmm) {
    if (hmm->type != vanilla) {
        st_errAbort("you gave me the wrong type of HMM");
    }
    StateMachine3Vanilla *sM3v = (StateMachine3Vanilla *)sM; // might be able to make this vanilla/echelon agnostic
    // made this 60 so that it loads alpha and beta probs
    for (int64_t i = 0; i < 60; i++) {
        sM3v->model.EMISSION_GAP_X_PROBS[i] = hmm->getTransitionsExpFcn(hmm, i, 0);
    }
}

// destructor
void vanillaHmm_destruct(Hmm *hmm) {
    VanillaHmm *vHmm = (VanillaHmm *)hmm;
    free(vHmm->matchModel);
    free(vHmm->scaledMatchModel);
    free(vHmm->kmerSkipBins);
    free(vHmm);
}

void vanillaHmm_writeToFile(Hmm *hmm, FILE *fileHandle) {
    /*
     * Format:
     * line 0: type \t stateNumber \t symbolSetSize \n
     * line 1: skip bins (alpha and beta) \t likelihood \n
     * line 2: [correlation coeff] \t [match model .. \t]  \n
     * line 3: [correlation coeff] [extra event matchModel]
     * See emissions_signal_loadPoreModel for description of matchModel
     * TODO might want to make poremodel and this more similar?
     */
    VanillaHmm *vHmm = (VanillaHmm *)hmm;

    // Line 0 - write the basic stuff to disk (positions are line:item#, not line:col)
    fprintf(fileHandle, "%i\t", vHmm->baseContinuousHmm.baseHmm.type); // type 0:0
    fprintf(fileHandle, "%lld\t", vHmm->baseContinuousHmm.baseHmm.stateNumber); // stateNumber 0:1
    fprintf(fileHandle, "%lld\t", vHmm->baseContinuousHmm.baseHmm.symbolSetSize); // symbolSetSize 0:2
    fprintf(fileHandle, "\n"); // newLine


    // Line 1 - write kmer skip bins to disk
    for (int64_t i = 0; i < 60; i++) {
        fprintf(fileHandle, "%f\t", vHmm->kmerSkipBins[i]); // kmer skip bins
    }
    fprintf(fileHandle, "%f\n", vHmm->baseContinuousHmm.baseHmm.likelihood); // likelihood, newline


    // Line 2 - write matchModel to disk
    int64_t nb_matchModelBuckets = 1 + (vHmm->baseContinuousHmm.baseHmm.symbolSetSize * MODEL_PARAMS);
    for (int64_t i = 0; i < nb_matchModelBuckets; i++) {
        fprintf(fileHandle, "%f\t", vHmm->matchModel[i]); // correlation coeff, matchModel
    }
    fprintf(fileHandle, "\n"); // newLine

    // Line 3 - write extra event model to disk
    for (int64_t i = 0; i < nb_matchModelBuckets; i++) {
        fprintf(fileHandle, "%f\t", vHmm->scaledMatchModel[i]); // correlation coeff, extra event matchModel
    }
    fprintf(fileHandle, "\n"); // newLine
}

Hmm *vanillaHmm_loadFromFile(const char *fileName) {

    // open file
    FILE *fH = fopen(fileName, "r");

    // line 0
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    int type;
    int64_t stateNumber, symbolSetSize;
    int64_t j = sscanf(stList_get(tokens, 0), "%i", &type); // type
    if (j != 1) {
        st_errAbort("Failed to parse type (int) from string: %s\n", string);
    }
    int64_t s = sscanf(stList_get(tokens, 1), "%lld", &stateNumber); // stateNumber
    if (s != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    int64_t n = sscanf(stList_get(tokens, 2), "%lld", &symbolSetSize); // symbolSetSize
    if (n != 1) {
        st_errAbort("Failed to parse symbol set size (int) from string: %s\n", string);
    }

    // make empty vanillaHmm
    Hmm *hmm = vanillaHmm_constructEmpty(0.0, stateNumber, symbolSetSize, type,
                                          vanillaHmm_addToKmerSkipBinExpectation,
                                          vanillaHmm_setKmerSkipBinExpectation,
                                          vanillaHmm_getKmerSkipBinExpectation);
    // Downcast
    VanillaHmm *vHmm = (VanillaHmm *)hmm;

    // cleanup
    free(string);
    stList_destruct(tokens);

    // kmer skip bins
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check
    if (stList_length(tokens) != 61) {
        st_errAbort("Did not find the correct number of kmer skip bins and/or likelihood\n");
    }
    // load
    for (int64_t i = 0; i < 60; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(vHmm->kmerSkipBins[i]));
        if (j != 1) {
            st_errAbort("Error parsing kmer skip bins from string %s\n", string);
        }
    }
    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(vHmm->baseContinuousHmm.baseHmm.likelihood));
    if (j != 1) {
        st_errAbort("error parsing likelihood from string %s\n", string);
    }

    // cleanup
    free(string);
    stList_destruct(tokens);

    // match model
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check
    int64_t nb_matchModelBuckets = 1 + (vHmm->baseContinuousHmm.baseHmm.symbolSetSize * MODEL_PARAMS);
    if (stList_length(tokens) != nb_matchModelBuckets) {
        st_errAbort("incorrect number of members for match model in HMM %s got %lld instead of %lld",
                    fileName, stList_length(tokens), nb_matchModelBuckets);

    }
    // load
    for (int64_t i = 0; i < nb_matchModelBuckets; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(vHmm->matchModel[i]));
        if (j != 1) {
            st_errAbort("error parsing match model for string %s", string);
        }
    }

    // cleanup
    free(string);
    stList_destruct(tokens);

    // and finally, the extra event match model
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check
    if (stList_length(tokens) != nb_matchModelBuckets) {
        st_errAbort("incorrect number of members for extra event match model in HMM %s got %lld instead of %lld",
                    fileName, stList_length(tokens), nb_matchModelBuckets);

    }
    // load
    for (int64_t i = 0; i < nb_matchModelBuckets; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(vHmm->scaledMatchModel[i]));
        if (j != 1) {
            st_errAbort("error parsing extra event match model for string %s", string);
        }
    }

    // cleanup
    free(string);
    stList_destruct(tokens);
    fclose(fH);

    return (Hmm *)vHmm;
}

///////////////////////////////////////////////// CORE FUNCTIONS //////////////////////////////////////////////////////
Hmm *hmmContinuous_loadSignalHmm(const char *fileName, StateMachineType type) {
    assert((type == vanilla) || (type == threeState));
    if (type == vanilla) {
        Hmm *hmm = vanillaHmm_loadFromFile(fileName);
        return hmm;
    }
    if (type == threeState) {
        Hmm *hmm = continuousPairHmm_loadFromFile(fileName);
        return hmm;
    }
    return 0;
}

void hmmContinuous_loadExpectations(StateMachine *sM, Hmm *hmm, StateMachineType type) { // todo rename this function
    assert((type == vanilla) || (type == threeState));
    if (type == vanilla) {
        vanillaHmm_loadKmerSkipBinExpectations(sM, hmm);
    }
    if (type == threeState) {
        continuousPairHmm_loadTransitionsAndKmerGapProbs(sM, hmm);
    }
}

void hmmContinuous_destruct(Hmm *hmm, StateMachineType type) {
    assert((type == vanilla) || (type == threeState));
    if (type == vanilla) {
        vanillaHmm_destruct(hmm);
    }
    if (type == threeState) {
        continuousPairHmm_destruct(hmm);
    }
}

Hmm *hmmContinuous_getEmptyHmm(StateMachineType type, double pseudocount) {
    assert((type == vanilla) || (type == threeState));
    if (type == vanilla) {
        Hmm *hmm = vanillaHmm_constructEmpty(pseudocount, 3, NUM_OF_KMERS, vanilla,
                                             vanillaHmm_addToKmerSkipBinExpectation,
                                             vanillaHmm_setKmerSkipBinExpectation,
                                             vanillaHmm_getKmerSkipBinExpectation);
        return hmm;
    }
    if (type == threeState) {
        Hmm *hmm = continuousPairHmm_constructEmpty(pseudocount, 3, NUM_OF_KMERS, threeState,
                                                    continuousPairHmm_addToTransitionsExpectation,
                                                    continuousPairHmm_setTransitionExpectation,
                                                    continuousPairHmm_getTransitionExpectation,
                                                    continuousPairHmm_addToKmerGapExpectation,
                                                    continuousPairHmm_setKmerGapExpectation,
                                                    continuousPairHmm_getKmerGapExpectation,
                                                    emissions_discrete_getKmerIndexFromKmer);
        return hmm;
    }
    return 0;
}

void hmmContinuous_normalize(Hmm *hmm, StateMachineType type) {
    assert((type == vanilla) || (type == threeState));
    if (type == vanilla) {
        vanillaHmm_normalizeKmerSkipBins(hmm);
    }
    if (type == threeState) {
        continuousPairHmm_normalize(hmm);
    }
}

void hmmContinuous_writeToFile(const char *outFile, Hmm *hmm, StateMachineType type) {
    assert((type == vanilla) || (type == threeState));
    FILE *fH = fopen(outFile, "w");
    if (type == vanilla) {
        vanillaHmm_writeToFile(hmm, fH);
    }
    if (type == threeState) {
        continuousPairHmm_writeToFile(hmm, fH);
    }
    fclose(fH);
}
