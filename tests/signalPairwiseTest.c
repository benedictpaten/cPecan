#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <discreteHmm.h>
#include <nanopore.h>
#include "randomSequences.h"
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "emissionMatrix.h"



static void test_signal_cell(CuTest *testCase) {
    StateMachineFunctions *sMfs = stateMachineFunctions_construct(emissions_signal_getKmerGapProb,
                                                                  emissions_signal_getEventGapProb,
                                                                  //emissions_signal_getLogPhiMatchProb);
                                                                  emissions_signal_getlogGaussPDFMatchProb);
    char *modelFile = stString_print("../../cPecan/models/template.eTable.model");
    StateMachine *sM = getSignalStateMachine3(modelFile, sMfs);
    double lowerF[sM->stateNumber], middleF[sM->stateNumber], upperF[sM->stateNumber], currentF[sM->stateNumber];
    double lowerB[sM->stateNumber], middleB[sM->stateNumber], upperB[sM->stateNumber], currentB[sM->stateNumber];
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
    int64_t testLength = 5;

    double fakeEventSeq[5] = {
                              60.032615, //ATGACA
                              60.332089, //TGACAC
                              61.618848, //GACACA
                              66.015805, //ACACAT
                              59.783408, //CACATT
                              };
    char *referenceSeq = "ATGACACATT";
    int64_t correctedLength = correctSeqLength(strlen(referenceSeq), event);
    CuAssertIntEquals(testCase, testLength, correctedLength);

    // make sequence objects
    Sequence *eventSeq = sequence_sequenceConstruct(testLength, fakeEventSeq, sequence_getEvent);
    Sequence *referSeq = sequence_sequenceConstruct(correctedLength, referenceSeq, sequence_getKmer);

    // test sequence_getEvent
    for (int64_t i = 0; i < testLength; i++) {
        CuAssertDblEquals(testCase, *(double *)eventSeq->get(eventSeq->elements, i), fakeEventSeq[i], 0.0);
    }

    // get one element from each sequence
    void *kX = referSeq->get(referSeq->elements, 4);
    void *eY = eventSeq->get(eventSeq->elements, 4);

    //Do forward
    cell_calculateForward(sM, lowerF, NULL, NULL, middleF, kX, eY, NULL);
    cell_calculateForward(sM, upperF, middleF, NULL, NULL, kX, eY, NULL);
    cell_calculateForward(sM, currentF, lowerF, middleF, upperF, kX, eY, NULL);
    //Do backward
    cell_calculateBackward(sM, currentB, lowerB, middleB, upperB, kX, eY, NULL);
    cell_calculateBackward(sM, upperB, middleB, NULL, NULL, kX, eY, NULL);
    cell_calculateBackward(sM, lowerB, NULL, NULL, middleB, kX, eY, NULL);
    double totalProbForward = cell_dotProduct2(currentF, sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, sM, sM->startStateProb);
    st_logInfo("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal
}

static void test_signal_diagonalDPCalculations(CuTest *testCase) {
    // make some DNA sequences
    //char *sX = "ATGACACATT";
    char *sX =   "ATGACATTCATT"; // has TT insert
    double sY[5] = {
            60.032615, //ATGACA
            60.332089, //TGACAC
            61.618848, //GACACA
            66.015805, //ACACAT
            59.783408, //CACATT
    };
    // make variables for the (corrected) length of the sequences
    int64_t lX = correctSeqLength(strlen(sX), event);
    int64_t lY = 5;
    // make Sequence objects
    Sequence *SsX = sequence_sequenceConstruct(lX, sX, sequence_getKmer);
    Sequence *SsY = sequence_sequenceConstruct(lY, sY, sequence_getEvent);

    // make stateMachine, forward and reverse DP matrices and banding stuff
    StateMachineFunctions *sMfs = stateMachineFunctions_construct(emissions_signal_getKmerGapProb,
                                                                  emissions_signal_getEventGapProb,
                                                                  emissions_signal_getlogGaussPDFMatchProb);
    char *modelFile = stString_print("../../cPecan/models/template.eTable.model");
    StateMachine *sM = getSignalStateMachine3(modelFile, sMfs);

    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, SsX->length, SsY->length, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    // Initialize Matrices
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
        diagonalCalculationForward(sM, i, dpMatrixForward, SsX, SsY);
    }
    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, SsX, SsY);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    // Test the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       SsX, SsY);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01);
    }

    // Now do the posterior probabilities, get aligned pairs with posterior match probs above threshold
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    // Perform alignment
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.2;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, SsX, SsY,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    // Make a list of the correct anchor points
    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                          (void (*)(void *)) stIntTuple_destruct);


    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 2));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(5, 3));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(6, 4));

    // make sure alignedPairs is correct
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, 5, (int) stList_length(alignedPairs));

    // clean up
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
}

static void test_scaleModel(CuTest *testCase) {
    char *modelFile = stString_print("../../cPecan/models/template.eTable.model");
    StateMachineFunctions *sMfs = stateMachineFunctions_construct(emissions_signal_getKmerGapProb,
                                                                  emissions_signal_getEventGapProb,
                                                                  emissions_signal_getlogGaussPDFMatchProb);
    StateMachine *sM = getSignalStateMachine3(modelFile, sMfs);

    char *npReadFile = stString_print("../../cPecan/tests/ZymoC_file1.npRead");
    NanoporeRead *npRead = loadNanoporeReadFromFile(npReadFile);

    emissions_signal_scaleModel(sM, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd);

    StateMachine *sM2 = getSignalStateMachine3(modelFile, sMfs); // unscaled model

    for (int64_t i = 0; i < sM->parameterSetSize * MODEL_PARAMS; i += 4) {
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i],
                          (sM2->EMISSION_MATCH_PROBS[i] * npRead->templateParams.scale + npRead->templateParams.shift),
                          0.0);
    }

}


static void test_signal_strandAlignmentNoBanding(CuTest *testCase) {
    /*
     * test strand event alignment to it's 2D read
     */
    char *ZymoReference = stString_print("../../cPecan/tests/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);

    // load NanoporeRead
    char *npReadFile = stString_print("../../cPecan/tests/ZymoC_file1.npRead");
    NanoporeRead *npRead = loadNanoporeReadFromFile(npReadFile);

    // make Sequence objects
    // 2D read sequence
    int64_t lX = correctSeqLength(npRead->readLength, event);
    Sequence *twoDreadSeq = sequence_sequenceConstruct(lX, npRead->twoDread, sequence_getKmer);

    // reference nucleotide sequence
    //int64_t lX = correctSeqLength(strlen(ZymoReferenceSeq), event);
    //Sequence *twoDreadSeq = sequence_sequenceConstruct(lX, ZymoReferenceSeq, sequence_getKmer);

    // event sequence
    int64_t lY = npRead->nbTemplateEvents;
    Sequence *templateEventSeq = sequence_sequenceConstruct(lY, npRead->templateEvents, sequence_getEvent);

    // load stateMachine and model
    char *modelFile = stString_print("../../cPecan/models/template.eTable.model");
    StateMachineFunctions *sMfs = stateMachineFunctions_construct(emissions_signal_getKmerGapProb,
                                                                  emissions_signal_getEventGapProb,
                                                                  emissions_signal_getlogGaussPDFMatchProb);
    StateMachine *sM = getSignalStateMachine3(modelFile, sMfs);
    emissions_signal_scaleModel(sM, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd);

    // matrix and bands
    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, twoDreadSeq->length, templateEventSeq->length, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    // Ini matrices
    // Initialize Matrices
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
        diagonalCalculationForward(sM, i, dpMatrixForward, twoDreadSeq, templateEventSeq);
    }
    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, twoDreadSeq, templateEventSeq);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);
    st_uglyf("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    // Test the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       twoDreadSeq, templateEventSeq);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.02);
    }

    // Now do the posterior probabilities, get aligned pairs with posterior match probs above threshold
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    // Perform alignment
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.4;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, twoDreadSeq, templateEventSeq,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }
    st_uglyf("No. aligned pairs: %lld\n", stList_length(alignedPairs));
    //for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
    //    stIntTuple *pair = stList_get(alignedPairs, i);
    //    int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
    //    st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
    //    st_uglyf("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
    //}
}


CuSuite *signalPairwiseTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_signal_cell);
    SUITE_ADD_TEST(suite, test_signal_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_scaleModel);
    SUITE_ADD_TEST(suite, test_signal_strandAlignmentNoBanding);
    return suite;
}