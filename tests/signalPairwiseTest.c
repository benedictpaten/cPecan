#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <nanopore.h>
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"


// brute force probability formulae
static double test_standardNormalPdf(double x) {
    double pi = 3.141592653589793;
    double inv_sqTwoPi = 1 / (sqrt(2*pi));
    double result = inv_sqTwoPi * exp(-(x * x) / 2);
    return result;
}

static double test_normalPdf(double x, double mu, double sigma) {
    double pi = 3.141592653589793;
    double inv_sqTwoPi = 1 / (sqrt(2*pi));
    double c = inv_sqTwoPi * (1/sigma);
    double a = (x - mu) / sigma;
    double result = c * exp(-0.5f * a * a);
    return result;
}

static double test_inverseGaussianPdf(double x, double mu, double lambda) {
    double pi = 3.141592653589793;
    double c = lambda / (2 * pi * pow(x, 3.0));
    double sqt_c = sqrt(c);
    double xmmu = x - mu;
    double result = sqt_c * exp((-lambda * xmmu * xmmu)/
                                (2 * mu * mu * x));
    return result;
}

static void test_getLogGaussPdfMatchProb(CuTest *testCase) {
    // standard normal distribution
    double eventModel[] = {0, 0, 1.0};
    double control = test_standardNormalPdf(0);
    char *kmer1 = "AAAAAA";
    double event1[] = {0};
    double test = emissions_signal_logGaussMatchProb(eventModel, kmer1, event1);
    double expTest = exp(test);
    CuAssertDblEquals(testCase, expTest, control, 0.001);
    CuAssertDblEquals(testCase, test, log(control), 0.001);

    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);
    double event2[] = {62.784241};
    double control2 = test_normalPdf(62.784241, sM->EMISSION_MATCH_PROBS[1], sM->EMISSION_MATCH_PROBS[2]);
    double test2 = emissions_signal_logGaussMatchProb(sM->EMISSION_MATCH_PROBS, kmer1, event2);
    CuAssertDblEquals(testCase, test2, log(control2), 0.001);
    stateMachine_destruct(sM);
}

static void test_bivariateGaussPdfMatchProb(CuTest *testCase) {
    // standard normal distribution
    double eventModel[] = {0, 0, 1.0, 0, 1.0};
    double control = test_standardNormalPdf(0);
    double controlSq = control * control;
    char *kmer1 = "AAAAAA";
    double event1[] = {0, 0}; // mean noise
    double test = emissions_signal_getBivariateGaussPdfMatchProb(eventModel, kmer1, event1);
    double eTest = exp(test);
    CuAssertDblEquals(testCase, controlSq, eTest, 0.001);

    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);
    double event2[] = {62.784241, 0.664989};
    double control2a = test_normalPdf(62.784241, sM->EMISSION_MATCH_PROBS[1], sM->EMISSION_MATCH_PROBS[2]);
    double control2b = test_normalPdf(0.664989, sM->EMISSION_MATCH_PROBS[3], sM->EMISSION_MATCH_PROBS[4]);
    double control2Sq = control2a * control2b;
    double test2 = emissions_signal_getBivariateGaussPdfMatchProb(sM->EMISSION_MATCH_PROBS, kmer1, event2);
    double eTest2 = exp(test2);
    CuAssertDblEquals(testCase, eTest2, control2Sq, 0.001);
    stateMachine_destruct(sM);
}

static void test_twoDistributionPdf(CuTest *testCase) {
    // load up a stateMachine
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);

    double event[] = {62.784241, 0.664989}; // level_mean and noise_mean for AAAAAA
    char *kmer1 = "AAAAAA";

    // get a sample match prob
    double test = emissions_signal_getEventMatchProbWithTwoDists(sM->EMISSION_MATCH_PROBS, kmer1, event);
    double control_level = test_normalPdf(62.784241, sM->EMISSION_MATCH_PROBS[1], sM->EMISSION_MATCH_PROBS[2]);
    double control_noise = test_inverseGaussianPdf(0.664989, sM->EMISSION_MATCH_PROBS[3], sM->EMISSION_MATCH_PROBS[5]);
    double control_prob = log(control_level) + log(control_noise);
    CuAssertDblEquals(testCase, test, control_prob, 0.001);
}


static void test_signal_cell(CuTest *testCase) {
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);
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

    double fakeEventSeq[15] = {
                              60.032615, 0.791316, 0.005, //ATGACA
                              60.332089, 0.620198, 0.012, //TGACAC
                              61.618848, 0.747567, 0.008, //GACACA
                              66.015805, 0.714290, 0.021, //ACACAT
                              59.783408, 1.128591, 0.002, //CACATT
                              };
    char *referenceSeq = "ATGACACATT";
    int64_t correctedLength = sequence_correctSeqLength(strlen(referenceSeq), event);
    CuAssertIntEquals(testCase, testLength, correctedLength);

    // make sequence objects
    Sequence *eventSeq = sequence_construct(testLength, fakeEventSeq, sequence_getEvent);
    Sequence *referSeq = sequence_construct(correctedLength, referenceSeq, sequence_getKmer2);

    // test sequence_getEvent
    for (int64_t i = 0; i < testLength; i++) {
        CuAssertDblEquals(testCase, *(double *)eventSeq->get(eventSeq->elements, i),
                          fakeEventSeq[i * NB_EVENT_PARAMS], 0.0);
    }

    // get one element from each sequence
    void *kX = referSeq->get(referSeq->elements, 1);
    void *eY = eventSeq->get(eventSeq->elements, 1);

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

static void test_echelon_cell(CuTest *testCase) {
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    st_uglyf("SENTINAL - about to construct stateMachine\n");
    StateMachine *sM = getEchelonStateMachine(modelFile);
    st_uglyf("SENTINAL - stateMachine constructed and poreModel loaded\n");
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

    double fakeEventSeq[15] = {
            60.032615, 0.791316, 0.005, //ATGACA
            60.332089, 0.620198, 0.012, //TGACAC
            61.618848, 0.747567, 0.008, //GACACA
            66.015805, 0.714290, 0.021, //ACACAT
            59.783408, 1.128591, 0.002, //CACATT
    };
    char *referenceSeq = "ATGACACATT";
    char *endPadding = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    int64_t correctedLength = sequence_correctSeqLength(strlen(referenceSeq), event);
    CuAssertIntEquals(testCase, testLength, correctedLength);
    referenceSeq = stString_print("%s%s", referenceSeq, endPadding);
    st_uglyf("SENTINAL - padded sequence %s\n", referenceSeq);

    // make sequence objects
    Sequence *referSeq = sequence_construct(correctedLength, referenceSeq, sequence_getKmer2);
    Sequence *eventSeq = sequence_construct(testLength, fakeEventSeq, sequence_getEvent);

    // test sequence_getEvent
    for (int64_t i = 0; i < testLength; i++) {
        st_uglyf("SENTINAL - testing getEvent %lld\n", i);
        CuAssertDblEquals(testCase, *(double *)eventSeq->get(eventSeq->elements, i),
                          fakeEventSeq[i * NB_EVENT_PARAMS], 0.0);
    }

    // get one element from each sequence
    void *kX = referSeq->get(referSeq->elements, 1);
    void *eY = eventSeq->get(eventSeq->elements, 1);

    st_uglyf("SENTINAL - got elements kX: %s, eY: %f\n", (char *) kX, *(double *) eY);


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
    st_uglyf("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal
}


static void test_signal_diagonalDPCalculations(CuTest *testCase) {
    // make some DNA sequences
    char *sX = "ATGACACATT";
    double sY[12] = {
            60.032615, 0.791316, 0.005, //ATGACA
            60.332089, 0.620198, 0.012, //TGACAC
            //61.618848, 0.747567, 0.008, //GACACA //skip
            66.015805, 0.714290, 0.021, //ACACAT
            59.783408, 1.128591, 0.002, //CACATT
    };
    // make variables for the (corrected) length of the sequences
    int64_t lX = sequence_correctSeqLength(strlen(sX), event);
    int64_t lY = 4;
    // make Sequence objects
    Sequence *SsX = sequence_construct(lX, sX, sequence_getKmer2);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent);

    // make stateMachine, forward and reverse DP matrices and banding stuff
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);

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
        p->threshold = 0.5;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, SsX, SsY,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    // Make a list of the correct anchor points
    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                          (void (*)(void *)) stIntTuple_destruct);


    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(3, 2));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(4, 3));


    // make sure alignedPairs is correct
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, 4, (int) stList_length(alignedPairs));

    // clean up
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
}

static void test_scaleModel(CuTest *testCase) {
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);

    char *npReadFile = stString_print("../../cPecan/tests/ZymoC_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    emissions_signal_scaleModel(sM, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd);

    // unpack npRead to make things easier
    double scale = npRead->templateParams.scale, shift = npRead->templateParams.shift,
        var = npRead->templateParams.var, scale_sd = npRead->templateParams.scale_sd,
        var_sd = npRead->templateParams.var_sd;


    StateMachine *sM2 = getSignalStateMachine3(modelFile); // unscaled model
    // TODO build this out to check everything
    for (int64_t i = 1; i < 1 + (sM->parameterSetSize * MODEL_PARAMS); i += 5) {
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i],
                          (sM2->EMISSION_MATCH_PROBS[i] * scale + shift), 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i+1],
                          sM2->EMISSION_MATCH_PROBS[i+1] * var, 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i+2],
                          sM2->EMISSION_MATCH_PROBS[i+2] * scale_sd, 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i+4],
                          sM2->EMISSION_MATCH_PROBS[i+4] * var_sd, 0.0);
        CuAssertDblEquals(testCase, sM->EMISSION_MATCH_PROBS[i+3],
                          sqrt(pow(sM->EMISSION_MATCH_PROBS[i+2], 3.0) / sM->EMISSION_MATCH_PROBS[i+4]),
                          0.0);
    }
}

static stList *compareAlignedPairs(stList *pairs1, stList *pairs2) {
    /*
     * doesn't work, need to remove probs
     */
    stSortedSet *sortedSet1 = stList_getSortedSet(pairs1, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    stSortedSet *sortedSet2 = stList_getSortedSet(pairs2, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    stSortedSet *intersection = stSortedSet_getIntersection(sortedSet1, sortedSet2);
    stList *list = stSortedSet_getList(intersection);
    return list;
}

static void checkAlignedPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    //st_uglyf("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
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
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

static void test_signal_strandAlignmentNoBanding2(CuTest *testCase) {
    // load reference
    char *ZymoReference = stString_print("../../cPecan/tests/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);

    // load NanoporeRead
    char *npReadFile = stString_print("../../cPecan/tests/ZymoC_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    // load stateMachine and model
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sM = getSignalStateMachine3(modelFile);
    emissions_signal_scaleModel(sM, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd); // clunky

    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.5; //0.4;

    stList *alignedPairs = getAlignedPairsWithoutBanding(sM, ZymoReferenceSeq, npRead->templateEvents, lX, lY, p,
                                                         sequence_getKmer2, sequence_getEvent, 0, 0);
    st_uglyf("there are %lld aligned pairs without banding, first time\n", stList_length(alignedPairs));
    checkAlignedPairs(testCase, alignedPairs, lX, lY);
}

static stList *scoreAnchorPairs(stList *anchorPairs, stList *alignedPairs) {
    /*
     * Selects the aligned pairs contained in anchor pairs.
     */
    stSortedSet *anchorPairsSet = stList_getSortedSet(anchorPairs, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    assert(stList_length(anchorPairs) == stSortedSet_size(anchorPairsSet));
    stList *scoredAnchorPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        stIntTuple *j = stIntTuple_construct2(stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2));
        if(stSortedSet_search(anchorPairsSet, j) != NULL) {
            stList_append(scoredAnchorPairs, stIntTuple_construct3(stIntTuple_get(aPair, 0), stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2)));
            stSortedSet_remove(anchorPairsSet, j);
        }
        stIntTuple_destruct(j);
    }

    //The following should not really be needed, and may be masking a bug/numerical precision issues
    stSortedSetIterator *it = stSortedSet_getIterator(anchorPairsSet);
    stIntTuple *pair;
    while((pair = stSortedSet_getNext(it))) {
        stList_append(scoredAnchorPairs, stIntTuple_construct3(0, stIntTuple_get(pair, 0), stIntTuple_get(pair, 1)));
    }
    stSortedSet_destructIterator(it);

    stSortedSet_destruct(anchorPairsSet);
    assert(stList_length(anchorPairs) == stList_length(scoredAnchorPairs));

    return scoredAnchorPairs;
}

static void test_signal_getAlignedPairsWithBanding(CuTest *testCase) {
    // load up test stuff and make stateMachine
    char *ZymoReference = stString_print("../../cPecan/tests/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *npReadFile = stString_print("../../cPecan/tests/ZymoC_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    char *templateModelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sMt = getSignalStateMachine3(templateModelFile);
    emissions_signal_scaleModel(sMt, npRead->templateParams.scale, npRead->templateParams.shift,
                                npRead->templateParams.var, npRead->templateParams.scale_sd,
                                npRead->templateParams.var_sd); // clunky

    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.5;

    // get anchors
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);
    // remap and filter
    nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(anchorPairs);

    // make Sequences
    Sequence *refSeq = sequence_construct(lX, ZymoReferenceSeq, sequence_getKmer2);
    Sequence *templateSeq = sequence_construct(lY, npRead->templateEvents, sequence_getEvent);

    // do alignment
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p, 0, 0);
    checkAlignedPairs(testCase, alignedPairs, lX, lY);
    st_uglyf("there are %lld aligned pairs using anchors\n", stList_length(alignedPairs));
    // do alignment without banding
    stList *alignedPairs2 = getAlignedPairsWithoutBanding(sMt, ZymoReferenceSeq, npRead->templateEvents, lX, lY, p,
                                                          sequence_getKmer2, sequence_getEvent, 0, 0);
    st_uglyf("there are %lld aligned pairs without banding\n", stList_length(alignedPairs2));
    checkAlignedPairs(testCase, alignedPairs2, lX, lY);
}

CuSuite *signalPairwiseTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getLogGaussPdfMatchProb);
    SUITE_ADD_TEST(suite, test_bivariateGaussPdfMatchProb);
    SUITE_ADD_TEST(suite, test_twoDistributionPdf);
    //SUITE_ADD_TEST(suite, test_signal_cell);
    SUITE_ADD_TEST(suite, test_echelon_cell);
    //SUITE_ADD_TEST(suite, test_signal_diagonalDPCalculations);
    //SUITE_ADD_TEST(suite, test_scaleModel);
    //SUITE_ADD_TEST(suite, test_signal_strandAlignmentNoBanding2);
    //SUITE_ADD_TEST(suite, test_signal_getAlignedPairsWithBanding);
    return suite;
}