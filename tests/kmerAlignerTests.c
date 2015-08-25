// Art Rand

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "randomSequences.h"
#include "shim.h"
#include "../inc/pairwiseAligner.h"
#include "../../sonLib/lib/sonLibCommon.h"
#include "../../sonLib/lib/sonLibRandom.h"
#include "../../sonLib/lib/CuTest.h"
#include "../../sonLib/lib/sonLibExcept.h"
#include "../../sonLib/lib/sonLibList.h"
#include "../inc/multipleAligner.h"
#include "../inc/emissionMatrix.h"



static void test_Kmers_cell(CuTest *testCase) {
    // construct the state machine
    StateMachine *sM = stateMachine5_kmer_construct(fiveState);

    // make arrays
    double lowerF[sM->stateNumber], middleF[sM->stateNumber], upperF[sM->stateNumber], currentF[sM->stateNumber];
    double lowerB[sM->stateNumber], middleB[sM->stateNumber], upperB[sM->stateNumber], currentB[sM->stateNumber];

    // set array values to start values
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        middleF[i] = sM->startStateProb(sM, i);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = sM->endStateProb(sM, i);
    }

    // Make some test sequences, these are the same sequences as in test_cell
    const char *testXseq = "AGCG";
    const char *testYseq = "AGTTCG";

    // Alternate test sequences
    //char* testXseq = "ATTA";
    //char* testYseq = "TGCT";
    //const char *testXseq = "AATT";
    //const char *testYseq = "AATT";

    // Initialize the sequence objects
    Sequence* xSeq = sequenceConstruct(4, testXseq, kmer);
    Sequence* ySeq = sequenceConstruct(6, testYseq, kmer);
    // Use the get function to retrieve one element of the sequence
    char* cX = xSeq->get(xSeq->elements, 1);
    char* cY = ySeq->get(ySeq->elements, 1);

    //Do forward
    cell_calculateForward(sM, lowerF, NULL, NULL, middleF, cX, cY, NULL);
    cell_calculateForward(sM, upperF, middleF, NULL, NULL, cX, cY, NULL);
    cell_calculateForward(sM, currentF, lowerF, middleF, upperF, cX, cY, NULL);
    //Do backward
    cell_calculateBackward(sM, currentB, lowerB, middleB, upperB, cX, cY, NULL);
    cell_calculateBackward(sM, upperB, middleB, NULL, NULL, cX, cY, NULL);
    cell_calculateBackward(sM, lowerB, NULL, NULL, middleB, cX, cY, NULL);
    double totalProbForward = cell_dotProduct2(currentF, sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, sM, sM->startStateProb);
    st_logInfo("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal

    // clean up
    sequenceDestroy(xSeq);
    sequenceDestroy(ySeq);
}

static void test_Kmers_diagonalDPCalculations(CuTest *testCase) {
    // These are the originals:
    const char *sX = "AGCG";
    const char *sY = "AGTTCG";

    // These are alternative sequences for testing
    //const char *sX = getRandomSequence(5);
    //const char *sY = evolveSequence(sX);

    // set lX and lY to the lengths of those sequences
    int64_t slX = strlen(sX);
    int64_t slY = strlen(sY);

    // construct a sequence struct from those sequences
    Sequence* sX2 = sequenceConstruct(slX, sX, kmer);
    Sequence* sY2 = sequenceConstruct(slY, sY, kmer);
    int64_t lX = sX2->length;
    int64_t lY = sY2->length;

    // construct a 5-state Kmer state machine, the forward and reverse DP Matrices, the band, the band
    // iterators and the anchor pairs
    StateMachine *sM = stateMachine5_kmer_construct(fiveState);
    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Initialise matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //initialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);

    //Forward algorithm
    for (int64_t i = 1; i <= lX + lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(sM, i, dpMatrixForward, sX2, sY2);
    }
    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, sX2, sY2);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    //Check the forward and back probabilities are about equal
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001);

    // Test the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       sX2, sY2);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01);
    }

    // Now do the posterior probabilities, get aligned pairs with posterior match probs above threshold
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    // Perform alignment
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.4; // used to be 0.2 needs to be doubled for kmers;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, sX2, sY2,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    // Make a list of the correct anchor points
    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                          (void (*)(void *)) stIntTuple_destruct);

    // these are correct for AGCG and AGTTCG
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 4));

    // make sure alignedPairs is correct
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        //printf("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, 3, (int) stList_length(alignedPairs));

    // clean up
    sequenceDestroy(sX2);
    sequenceDestroy(sY2);
}

// repeat of getRandomAnchorPairs for Kmer aligner tests
static stList *getRandomAnchorPairs_KmerTests(int64_t lX, int64_t lY) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    int64_t x = -1;
    int64_t y = -1;
    while (1) {
        x += st_randomInt(1, 20);
        y += st_randomInt(1, 20);
        if (x >= lX || y >= lY) {
            break;
        }
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY);
        stList_append(anchorPairs, stIntTuple_construct2(x, y));
    }
    return anchorPairs;
}

static void checkAlignedPairs_kmers(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                (void (*)(void *)) stIntTuple_destruct
                                                );
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 3);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX+1); // remove offset for aligned kmers
        CuAssertTrue(testCase, y < lY+1); // at the end of the sequences

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

static void test_Kmers_getAlignedPairsWithBanding(CuTest *testCase) {
    for (int64_t test = 0; test < 3; test++) {
        // Get two random sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);

        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);
        // Make sequence objects
        Sequence* sX2 = sequenceConstruct(lX, sX, kmer);
        Sequence* sY2 = sequenceConstruct(lY, sY, kmer);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->traceBackDiagonals = st_randomInt(1, 10);
        p->minDiagsBetweenTraceBack = p->traceBackDiagonals + st_randomInt(2, 10);
        p->diagonalExpansion = st_randomInt(0, 10) * 2;
        StateMachine *sM = stateMachine5_kmer_construct(fiveState);
        stList *anchorPairs = getRandomAnchorPairs_KmerTests(lX, lY);
        stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        void *extraArgs[1] = { alignedPairs };
        getPosteriorProbsWithBanding(sM,                                     //state machine
                                     anchorPairs,                            //(random) anchor pairs
                                     sX2, sY2,                               //sequence objects
                                     p,                                      //params
                                     0, 0,                                   //ragged left and right end booleans
                                     diagonalCalculationPosteriorMatchProbs, //posterior probability function
                                     extraArgs);                             //bin for aligned pairs

        //Check the aligned pairs.
        checkAlignedPairs_kmers(testCase, alignedPairs, lX, lY);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        sequenceDestroy(sX2);
        sequenceDestroy(sY2);
        stList_destruct(alignedPairs);
    }
}

static void test_kmers_getAlignedPairs(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(1, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        Sequence* SsX = sequenceConstruct(lX, sX, kmer);
        Sequence* SsY = sequenceConstruct(lY, sY, kmer);
        //st_logInfo("Sequence X to align: %s END\n", sX);
        //st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold=0.02;
        StateMachine *sM = stateMachine5_kmer_construct(fiveState);
        stList *alignedPairs = getAlignedPairs(sM, SsX, SsY, kmer, p, 0, 0); // this is where the magic happens

        //Check the aligned pairs.
        checkAlignedPairs_kmers(testCase, alignedPairs, SsX->length, SsY->length);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        sequenceDestroy(SsX);
        sequenceDestroy(SsY);
        stList_destruct(alignedPairs);
    }
}

static void test_kmers_getAlignedPairsWithRaggedEnds(CuTest *testCase) {
    for (int64_t test = 0; test < 1000; test++) {
        //Make a core sequence and then add a random portion to either end of it
        int64_t coreLength = 100, randomPortionLength = 100;
        char *sX = getRandomSequence(coreLength);
        char *sY = stString_print("%s%s%s", getRandomSequence(randomPortionLength), sX,
                                  getRandomSequence(randomPortionLength)); //x with an extra bit at the end.
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        Sequence* SsX = sequenceConstruct(lX, sX, kmer);
        Sequence* SsY = sequenceConstruct(lY, sY, kmer);

        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        // Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine5_kmer_construct(fiveState);
        stList *alignedPairs = getAlignedPairs(sM, SsX, SsY, kmer, p, 1, 1);
        alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, SsX->repr, SsY->repr, 0.2);

        // Check the aligned pairs.
        checkAlignedPairs_kmers(testCase, alignedPairs, strlen(SsX->repr), strlen(SsY->repr));
        CuAssertIntEquals(testCase, stList_length(alignedPairs), coreLength-1);
        for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
            stIntTuple *j = stList_get(alignedPairs, i);
            CuAssertTrue(testCase, stIntTuple_length(j) == 3);

            int64_t x = stIntTuple_get(j, 1);
            int64_t y = stIntTuple_get(j, 2);
            CuAssertIntEquals(testCase, x + randomPortionLength, y);
        }

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        sequenceDestroy(SsX);
        sequenceDestroy(SsY);
        stList_destruct(alignedPairs);
    }
}

static void test_Kmer_hmm(CuTest *testCase, StateMachineType stateMachineType) {
    //Expectation object
    Hmm *hmm = hmm_Kmer_constructEmpty(0.0, stateMachineType);
    //Add some transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_addToTransitionExpectation(hmm, from, to, from * hmm->stateNumber + to);
        }
    }

    //Add some emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < NUM_OF_KMERS; x++) {
            for (int64_t y = 0; y < NUM_OF_KMERS; y++) {
                double dummy = state * SYMBOL_NUMBER * SYMBOL_NUMBER + x * SYMBOL_NUMBER + y;
                hmm_Kmer_addToEmissionsExpectation(hmm, state, x, y, dummy);
            }
        }
    }

    //Write to a file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile)); //Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    hmm_Kmer_write(hmm, fH);
    fclose(fH);
    hmm_destruct(hmm);

    //Load from a file
    hmm = hmm_Kmer_loadFromFile(tempFile);
    stFile_rmrf(tempFile);

    //Check the transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            CuAssertTrue(testCase, hmm_getTransition(hmm, from, to) == from * hmm->stateNumber + to);
        }
    }

    //Check the emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < NUM_OF_KMERS; x++) {
            for (int64_t y = 0; y < NUM_OF_KMERS; y++) {
                double expected = state * SYMBOL_NUMBER * SYMBOL_NUMBER + x * SYMBOL_NUMBER + y;
                double actual = hmm_Kmer_getEmissionsExpectation(hmm, state, x, y);
                CuAssertDblEquals(testCase, expected, actual, 0.001);
                //CuAssertTrue(testCase,
                //             hmm_Kmer_getEmissionsExpectation(hmm, state, x, y) == state * SYMBOL_NUMBER * SYMBOL_NUMBER + x * SYMBOL_NUMBER + y);
            }
        }
    }

    //Normalise
    hmm_Kmer_normalise(hmm);

    //Recheck transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            double z = from * hmm->stateNumber * hmm->stateNumber + (hmm->stateNumber * (hmm->stateNumber - 1)) / 2;
            CuAssertDblEquals(testCase, (from * hmm->stateNumber + to) / z, hmm_getTransition(hmm, from, to), 0.0);
        }
    }

    //Recheck the emissions
    double z[5] = {45000, 60625, 76250, 91875, 107500}; // correct totals for each state
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < NUM_OF_KMERS; x++) {
            for (int64_t y = 0; y < NUM_OF_KMERS; y++) {
                CuAssertTrue(testCase,
                             hmm_Kmer_getEmissionsExpectation(hmm, state, x, y) == (state * SYMBOL_NUMBER * SYMBOL_NUMBER + x * SYMBOL_NUMBER + y)/z[state]);
            }
        }
    }
    //Clean up
    hmm_destruct(hmm);
}

static void test_Kmer_em(CuTest *testCase, StateMachineType stateMachineType) {
    for (int64_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(10, 200));
        char *sY = evolveSequence(sX); //stString_copy(seqX);

        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);

        Sequence* SsX = sequenceConstruct(lX, sX, kmer);
        Sequence* SsY = sequenceConstruct(lY, sY, kmer);

        //Currently starts from random model and iterates.
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        double pLikelihood = -INFINITY;
        Hmm *hmm = hmm_Kmer_constructEmpty(0.0, stateMachineType);
        hmm_Kmer_randomise(hmm);
        StateMachine *sM = hmm_Kmer_getStateMachine(hmm);
        hmm_destruct(hmm);

        for (int64_t iteration = 0; iteration < 10; iteration++) {
            hmm = hmm_Kmer_constructEmpty(0.000000000001, stateMachineType); //The tiny pseudo count prevents overflow 0.000000000001
            getExpectations(sM, hmm, SsX, SsY, p, 0, 0);
            hmm_Kmer_normalise(hmm);
            //Log stuff
            // Transitions
            for (int64_t from = 0; from < sM->stateNumber; from++) {
                for (int64_t to = 0; to < sM->stateNumber; to++) {
                    //printf("Transition from %" PRIi64 " to %" PRIi64 " has expectation %f\n", from, to,
                    //           hmm_getTransition(hmm, from, to));
                }
            }
            // Emissions
            for (int64_t x = 0; x < NUM_OF_KMERS; x++) {
                for (int64_t y = 0; y < NUM_OF_KMERS; y++) {
                    //printf("Emission x %" PRIi64 " y %" PRIi64 " has expectation %f\n", x, y,
                    //           hmm_Kmer_getEmissionsExpectation(hmm, sM->matchState, x, y));
                }
            }

            st_logInfo("->->-> Got expected likelihood %f for trial %" PRIi64 " and  iteration %" PRIi64 "\n",
                       hmm->likelihood, test, iteration);
            //printf("->->-> Got expected likelihood %f for trial %" PRIi64 " and  iteration %" PRIi64 "\n",
            //           hmm->likelihood, test, iteration);
            assert(pLikelihood <= hmm->likelihood * 0.85); //Used to be 0.95 for nucleotide model
            CuAssertTrue(testCase, pLikelihood <= hmm->likelihood * 0.85); // same here
            pLikelihood = hmm->likelihood;
            stateMachine_destruct(sM);
            sM = hmm_getStateMachine(hmm);
            hmm_destruct(hmm);
        }

        //Cleanup
        pairwiseAlignmentBandingParameters_destruct(p);
        sequenceDestroy(SsX);
        sequenceDestroy(SsY);
        free(sX);
        free(sY);
    }
}

static void test_Kmer_hmm_5State(CuTest *testCase) {
    test_Kmer_hmm(testCase, fiveState);
}

static void test_Kmer_em_5State(CuTest *testCase) {
    test_Kmer_em(testCase, fiveState);
}

CuSuite* kmerTestSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_Kmers_cell);
    SUITE_ADD_TEST(suite, test_Kmers_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_Kmers_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_kmers_getAlignedPairs);
    SUITE_ADD_TEST(suite, test_kmers_getAlignedPairsWithRaggedEnds);
    SUITE_ADD_TEST(suite, test_Kmer_hmm_5State);
    SUITE_ADD_TEST(suite, test_Kmer_em_5State);

    return suite;
}
