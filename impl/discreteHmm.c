#include <stdlib.h>
#include <stdio.h>
#include "bioioC.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "discreteHmm.h"

// Construct
Hmm *hmmDiscrete_constructEmpty(double pseudocount, int64_t stateNumber, int64_t symbolSetSize, StateMachineType type,
                                void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
                                void (*addEmissionsExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                void (*setEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p),
                                double (*getEmissionExpFcn)(Hmm *hmm, int64_t state, int64_t x, int64_t y),
                                int64_t (*getElementIndexFcn)(void *)) {
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
    // indexing
    hmmD->getElementIndexFcn = getElementIndexFcn;

    return hmmD;
}
// Transitions
void hmmDiscrete_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] += p;
}

void hmmDiscrete_setTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] = p;
}

double hmmDiscrete_getTransitionExpectation(Hmm *hmm, int64_t from, int64_t to) {
    return hmm->transitions[from * hmm->stateNumber + to];
}

// Emissions
void hmmDiscrete_addToEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p) {
    int64_t tableIndex = x * hmm->symbolSetSize + y;
    hmm->emissions[(state * hmm->matrixSize) + tableIndex] += p;
}

void hmmDiscrete_setEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y, double p) {
    int64_t tableIndex = x * hmm->symbolSetSize + y;
    hmm->emissions[(state * hmm->matrixSize) + tableIndex] = p;
}

double hmmDiscrete_getEmissionExpectation(Hmm *hmm, int64_t state, int64_t x, int64_t y) {
    int64_t tableIndex = x * hmm->symbolSetSize + y;
    return hmm->emissions[(state * hmm->matrixSize) + tableIndex];
}

// Randomize/Normalize
void hmmDiscrete_randomize(Hmm *hmmD) {
    // Transitions
    for (int64_t from = 0; from < hmmD->stateNumber; from++) {
        for (int64_t to = 0; to < hmmD->stateNumber; to++) {
            hmmDiscrete_setTransitionExpectation(hmmD, from, to, st_random());
        }
    }
    // Emissions
    for (int64_t state = 0; state < hmmD->stateNumber; state++) {
        for (int64_t x = 0; x < hmmD->symbolSetSize; x++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                hmmDiscrete_setEmissionExpectation(hmmD, state, x, y, st_random());
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
            total += hmmDiscrete_getTransitionExpectation(hmmD, from, to);
        }
        for (int64_t to = 0; to < hmmD->stateNumber; to++) {
            double newProb = hmmDiscrete_getTransitionExpectation(hmmD, from, to) / total;
            hmmDiscrete_setTransitionExpectation(hmmD, from, to, newProb);
        }
    }
    for (int64_t state = 0; state < hmmD->stateNumber; state++) {
        double total = 0.0;
        for (int64_t x = 0; x < hmmD->symbolSetSize; x++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                total += hmmDiscrete_getEmissionExpectation(hmmD, state, x, y);
            }
        }
        for (int64_t x = 0; x < hmmD->symbolSetSize; x ++) {
            for (int64_t y = 0; y < hmmD->symbolSetSize; y++) {
                double newProb = hmmDiscrete_getEmissionExpectation(hmmD, state, x, y) / total;
                hmmDiscrete_setEmissionExpectation(hmmD, state, x, y, newProb);
            }
        }
    }
}

// writers
void hmmDiscrete_write(Hmm *hmmD, FILE *fileHandle) {
    /*
     * Format:
     * type \t stateNumber \t symbolSetSize \n
     * transition \t transition ... likelihood \n
     * emission \t emission \t ... \n
     */
    // basics
    fprintf(fileHandle, "%i\t", hmmD->type); // 0
    fprintf(fileHandle, "%lld\t", hmmD->stateNumber); // 1
    fprintf(fileHandle, "%lld\t", hmmD->symbolSetSize); // 2
    fprintf(fileHandle, "\n");
    // transitions
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->stateNumber; i++) {
        fprintf(fileHandle, "%f\t", hmmD->transitions[i]);
    }
    // likelihood
    fprintf(fileHandle, "%f\n", hmmD->likelihood);
    // emissions
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->matrixSize; i++) {
        fprintf(fileHandle, "%f\t", hmmD->emissions[i]);
    }
    fprintf(fileHandle, "\n");
}

// Loaders
Hmm *hmmDiscrete_loadFromFile(const char *fileName) {
    // open the file
    FILE *fH = fopen(fileName, "r");
    // read the first line of the file and split at whitespace
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);

    // improper input check
    if (stList_length(tokens) < 2) {
        st_errAbort("Got an empty line in the input state machine file %s\n", fileName); // this is the wrong error message?
    }
    // setup
    int type;
    int64_t stateNumber, parameterSetSize;
    // get StateMachineType
    int64_t j = sscanf(stList_get(tokens, 0), "%i", &type);
    if (j != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    // get stateNumber
    int64_t s = sscanf(stList_get(tokens, 1), "%lld", &stateNumber);
    if (s != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    // get parameterSetSize
    int64_t n = sscanf(stList_get(tokens, 2), "%lld", &parameterSetSize);
    if (n != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }

    // make empty Hmm
    Hmm *hmmD = hmmDiscrete_constructEmpty(0.0, stateNumber, parameterSetSize, type,
                                          hmmDiscrete_addToTransitionExpectation,
                                          hmmDiscrete_setTransitionExpectation,
                                          hmmDiscrete_getTransitionExpectation,
                                          hmmDiscrete_addToEmissionExpectation,
                                          hmmDiscrete_setEmissionExpectation,
                                          hmmDiscrete_getEmissionExpectation,
                                          emissions_getBaseIndex);
    // cleanup setup line
    free(string);
    stList_destruct(tokens);

    // Transitions
    // parse transitions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check for the correct number of transitions
    if (stList_length(tokens) != hmmD->stateNumber * hmmD->stateNumber + 1) { // + 1 bc. likelihood is also on that line
        st_errAbort(
                "Got the wrong number of transitions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmmD->stateNumber * hmmD->stateNumber + 1);
    }
    // load transitions
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->stateNumber; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hmmD->transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    // load likelihood
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(hmmD->likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }
    // Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    // Now parse the emissions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // check for the correct number of emissions
    if (stList_length(tokens) != hmmD->stateNumber * hmmD->matrixSize) {
        st_errAbort(
                "Got the wrong number of emissions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmmD->matrixSize);
    }

    // load emissions
    for (int64_t i = 0; i < hmmD->stateNumber * hmmD->matrixSize; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hmmD->emissions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse emission prob (float) from string: %s\n", string);
        }
    }
    // Final cleanup
    free(string);
    stList_destruct(tokens);
    fclose(fH);

    return hmmD;
}


//static void hmmDiscrete_loadSymmetric()

// Housekeeping
void hmmDiscrete_destruct(Hmm *hmmD) {
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

