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

    // These are the same sequences as in test_cell
    const char *testXseq = "AGCG";
    const char *testYseq = "AGTTCG";

    // Alternate test sequences
    //char* testXseq = "ATTA";
    //char* testYseq = "TGCT";
    //const char *testXseq = "AATT";
    //const char *testYseq = "AATT";

    Sequence* xSeq = sequenceConstruct(4, testXseq, kmer);
    Sequence* ySeq = sequenceConstruct(6, testYseq, kmer);

    //char* xkmer = xSeq->get(xSeq->elements, 0);
    //char* ykmer = ySeq->get(ySeq->elements, 0);
    //printf("%s is the kmer at position %d\n", xkmer, 0);
    //printf("%s is the kmer at position %d\n", ykmer, 0);

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
    //printf("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal
}

static void test_Kmers_diagonalDPCalculations(CuTest *testCase) {
    // Sets up a complete matrix for the following example and checks the total
    // marginal probability and the posterior probabilities of the matches

    // make some simple DNA sequences
    // These are the originals:
    const char *sX = "AGCG";
    const char *sY = "AGTTCG";

    // These are not the originals
    //const char *sX = getRandomSequence(5);
    //const char *sY = evolveSequence(sX);

    // set lX and lY to the lengths of those sequences
    // NOTE The length of kmers is the length of bases-1
    int64_t slX = strlen(sX);
    int64_t slY = strlen(sY);


    // construct a sequence struct from those sequences and assign the get function as get base
    Sequence* sX2 = sequenceConstruct(slX, sX, kmer);
    Sequence* sY2 = sequenceConstruct(slY, sY, kmer);
    int64_t lX = sX2->length;
    int64_t lY = sY2->length;


    // construct a 5-state state machine, the forward and reverse DP Matrices, the band, the band
    // iterators and the anchor pairs
    StateMachine *sM = stateMachine5_kmer_construct(fiveState);
    //printf("just finished running stateMachine5 kmer construct, states: %lld\n", sM->stateNumber);

    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    //printf("just finished running dpMatrix_construct, twice\n");

    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    //printf("just finished running band_construct\n");
    BandIterator *bandIt = bandIterator_construct(band);
    //printf("just finished running bandIterator_construct\n");

    //Initialise matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //initialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    //printf("just finished initializing matrices\n");

    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);
    //printf("just finished initializing values\n");

    //Forward algorithm
    //printf("\n-->At forward algorithm\n");
    for (int64_t i = 1; i <= lX + lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(sM, i, dpMatrixForward, sX2, sY2);
    }
    //Backward algorithm
    //printf("\n-->At Backward algorithm\n");
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, sX2, sY2);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    //printf("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    //Check the forward and back probabilities are about equal

    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001);

    // Test calculating the posterior probabilities along the diagonals of the
    // matrix.
    //printf("\n-->Calculating posterior probabilities\n");
    for (int64_t i = 0; i <= lX + lY; i++) {
        //Calculate the total probs
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       sX2, sY2);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01);
    }

    //Now do the posterior probabilities
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    // aligned pairs has length 0 here, just constructed
    void *extraArgs[1] = { alignedPairs };
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

    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 4));


    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        //printf("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }

    CuAssertIntEquals(testCase, 3, (int) stList_length(alignedPairs));

}

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
    //st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    //printf("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
            (void (*)(void *)) stIntTuple_destruct);
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
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);


        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);

        //st_logInfo("Sequence X to align: %s END\n", sX);
        //st_logInfo("Sequence Y to align: %s END\n", sY);
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
        getPosteriorProbsWithBanding(sM,
                                     anchorPairs,
                                     sX2, sY2,
                                     p,
                                     0, 0,
                                     diagonalCalculationPosteriorMatchProbs, extraArgs);
        //Check the aligned pairs.
        checkAlignedPairs_kmers(testCase, alignedPairs, lX, lY);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        // TODO need sequence destruct function
        //free(sX2->elements);
        //free(sY2->elements);
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
        //printf("Sequence X to align: %s END\n", sX);
        //printf("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold=0.02;
        StateMachine *sM = stateMachine5_kmer_construct(fiveState);

        stList *alignedPairs = getAlignedPairs(sM, SsX, SsY, kmer, p, 0, 0);

        //Check the aligned pairs.
        checkAlignedPairs_kmers(testCase, alignedPairs, SsX->length, SsY->length);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        stList_destruct(alignedPairs);
    }
}


static void test_kmers_getAlignedPairsWithRaggedEnds(CuTest *testCase) {
    for (int64_t test = 0; test < 1000; test++) {
        //Make a pair of sequences
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
        //printf("Sequence X to align: %s END\n", SsX->repr);
        //printf("Sequence Y to align: %s END\n", SsY->repr);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine5_kmer_construct(fiveState);
        stList *alignedPairs = getAlignedPairs(sM, SsX, SsY, kmer, p, 1, 1);

        //printf("Before filtering alignedPairs Length: %lld\n", (int64_t) stList_length(alignedPairs));
        alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, SsX->repr, SsY->repr, 0.4);
        //printf("After filtering alignedPairs Length: %lld\n", (int64_t) stList_length(alignedPairs));

        //Check the aligned pairs.
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
        stList_destruct(alignedPairs);
    }
}


CuSuite* kmerTestSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_Kmers_cell);
    SUITE_ADD_TEST(suite, test_Kmers_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_Kmers_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_kmers_getAlignedPairs);
    SUITE_ADD_TEST(suite, test_kmers_getAlignedPairsWithRaggedEnds);
    return suite;
}
