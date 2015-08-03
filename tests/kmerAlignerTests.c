// Art Rand

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "shim.h"
#include "CuTest.h"
#include "../inc/shim.h"
#include "../../sonLib/lib/CuTest.h"
#include "../inc/stateMachine.h"
#include "../inc/pairwiseAligner.h"
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

    Sequence* xSeq = sequenceConstruct(4, testXseq, getKmer);
    Sequence* ySeq = sequenceConstruct(4, testYseq, getKmer);

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
    // TODO TODO TODO
    // Make some simpler sequences and walk thought the alignment process.
    // Right now I suspect that the problem is at the end of the alignment ie.
    // aligning the final kmers.
    //const char *sX = "AGTG";
    //const char *sY = "AGTG";

    // set lX and lY to the lengths of those sequences
    int64_t lX = strlen(sX)-1;
    int64_t lY = strlen(sY)-1;

    // construct a sequence struct from those sequences and assign the get function as get base
    Sequence* sX2 = sequenceConstruct(lX+1, sX, getKmer);
    Sequence* sY2 = sequenceConstruct(lY+1, sY, getKmer);

    // construct a 5-state state machine, the forward and reverse DP Matrices, the band, the band
    // iterators and the anchor pairs
    StateMachine *sM = stateMachine5_kmer_construct(fiveState);
    printf("just finished running stateMachine5 kmer construct, states: %lld\n", sM->stateNumber);

    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    printf("just finished running dpMatrix_construct, twice\n");

    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    printf("just finished running band_construct\n");
    BandIterator *bandIt = bandIterator_construct(band);
    printf("just finished running bandIterator_construct\n");

    //Initialise matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //initialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    printf("just finished initializing matrices\n");

    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);
    printf("just finished initializing values\n");

    //Forward algorithm
    printf("\n-->At forward algorithm\n");
    for (int64_t i = 1; i <= lX + lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(sM, i, dpMatrixForward, sX2, sY2);
    }
    //Backward algorithm
    printf("\n-->At Backward algorithm\n");
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, sX2, sY2);
    }

    //Calculate total probabilities
    // NOTE this is where the problem is
    double totalProbForward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    printf("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    //Check the forward and back probabilities are about equal

    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001);

    // Test calculating the posterior probabilities along the diagonals of the
    // matrix.
    printf("\n-->Calculating posterior probabilities\n");
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
        p->threshold = 0.4 // used to be 0.2 needs to be doubled for kmers;
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

CuSuite* kmerTestSuite() {
    CuSuite* suite = CuSuiteNew();

    //SUITE_ADD_TEST(suite, test_Kmers_cell);
    SUITE_ADD_TEST(suite, test_Kmers_diagonalDPCalculations);



    return suite;
}
