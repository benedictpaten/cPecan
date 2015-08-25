// Art Rand


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
//#include "shim.h"
#include "CuTest.h"
#include "randomSequences.h"
//#include "../inc/shim.h"
#include "../../sonLib/lib/CuTest.h"
#include "../inc/stateMachine.h"
#include "../inc/pairwiseAligner.h"
#include "../inc/emissionMatrix.h"

static void test_sequenceConstruct(CuTest* testCase) {
    char *tS = getRandomSequence(100);
    Sequence* testSequence = sequenceConstruct(100, tS, nucleotide);
    CuAssertIntEquals(testCase, 100, testSequence->length);
    CuAssertStrEquals(testCase, tS, testSequence->repr);
    free(tS);
    sequenceDestroy(testSequence);
}


// test chars (nucleotides)
static void test_chars(CuTest* testCase) {
    char* testSequence = "GATACA";
    int dnaLen = (int) strlen(testSequence);
    Sequence* charSeq = sequenceConstruct(dnaLen, testSequence, nucleotide);
    for (int c = 0; c < dnaLen; c++) {
        char* b = charSeq->get(charSeq->elements, c);
        CuAssertStrEquals(testCase, &testSequence[c], b);
    }
    sequenceDestroy(testSequence);
}

static void test_kmers(CuTest* testCase) {
    char* testSequence = "GATACAGATACAGATACA";
    int len = (int) strlen(testSequence);
    char testKmers[14][6] = {"GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                             "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                             "GATAC", "ATACA"};
    Sequence* kmerSeq = sequenceConstruct(len, testSequence, kmer);
    for (int x = 0; x < 5; x++) {
        char* kmer = kmerSeq->get(kmerSeq->elements, x);
        CuAssertStrEquals(testCase, &testKmers[x], kmer);
    }
}
/*
static void test_events(CuTest* testCase) {
    char testKmers[12][6] = { "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                              "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA" };
    double testMeans[] = { 12.12, 23.23, 34.34, 45.45, 56.56, 67.67,
                           12.12, 23.23, 34.34, 45.45, 56.56, 67.67 };

    void* eventSequence = eventSequenceConstruct(12, testMeans, *testKmers);
    Sequence* eventSeq = sequenceConstruct(12, eventSequence, event);

    for (int x = 0; x < 12; x++) {
        Event *ev = eventSeq->get(eventSeq->elements, x);
        // test for kmer
        CuAssertStrEquals(testCase, &testKmers[x], ev->kmer);
        CuAssertDblEquals(testCase, testMeans[x], ev->mean, 0.0);
    }
}
*/

static void test_getKmerIndex(CuTest* testCase) {
    char kmer2Set[25][3] = {"AA", "AC", "AG", "AT", "AN", "CA", "CC", "CG", "CT", "CN",
                           "GA", "GC", "GG", "GT", "GN", "TA", "TC", "TG", "TT", "TN",
                           "NA", "NC", "NG", "NT", "NN"};

    for (int64_t i = 0; i < NUM_OF_KMERS; i++) {
        int64_t kmer_index_result = getKmerIndex(kmer2Set[i]);
        CuAssertIntEquals(testCase, kmer_index_result, i);
    }
}

CuSuite* shimTestSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_sequenceConstruct);
    SUITE_ADD_TEST(suite, test_getKmerIndex);
//    SUITE_ADD_TEST(suite, test_chars);
//    SUITE_ADD_TEST(suite, test_kmers);
//    SUITE_ADD_TEST(suite, test_events);

    return suite;
}


/*
 * static void test_getKmerIndex(CuTest* testCase) {
    char* testA = "A";
    int64_t iA = getKmerIndex(testA);

    char* testC = "C";
    int64_t iC = getKmerIndex(testC);

    char* testG = "G";
    int64_t iG = getKmerIndex(testG);

    char* testT = "T";
    int64_t iT = getKmerIndex(testT);

    char* testN = "N";
    int64_t iN = getKmerIndex(testN);

    printf("%s, %lld\n", testA, iA);
    printf("%s, %lld\n", testC, iC);
    printf("%s, %lld\n", testG, iG);
    printf("%s, %lld\n", testT, iT);
    printf("%s, %lld\n", testN, iN);
}
 */