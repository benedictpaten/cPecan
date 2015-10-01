#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "randomSequences.h"
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "emissionMatrix.h"
#include "../../sonLib/lib/CuTest.h"
#include "../inc/pairwiseAligner.h"
#include "../../sonLib/lib/sonLibExcept.h"
#include "../../sonLib/lib/sonLibCommon.h"
#include "../../sonLib/lib/sonLibList.h"
#include "../../sonLib/lib/sonLibTuples.h"
#include "../../sonLib/lib/sonLibRandom.h"
#include "../inc/emissionMatrix.h"
#include "../../sonLib/lib/sonLibString.h"
#include "../../sonLib/lib/sonLibFile.h"
#include "../../sonLib/lib/sonLibSortedSet.h"
#include "../inc/multipleAligner.h"
#include "../inc/discreteHmm.h"


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

    // make variables for the (corrected) length of the sequences

    // make Sequence objects

    // make stateMachine, forward and reverse DP matrices and banding stuff


}

CuSuite *signalPairwiseTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_signal_cell);
    SUITE_ADD_TEST(suite, test_signal_diagonalDPCalculations);
    return suite;
}