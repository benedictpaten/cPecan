// Art Rand


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "shim.h"
#include "CuTest.h"
#include "../inc/shim.h"

// test chars (nucleotides)
static void test_chars(CuTest* testCase) {
    char* testSequence = "GATACA";
    int dnaLen = strlen(testSequence);
    Sequence* charSeq = sequenceConstruct(dnaLen, testSequence, getBase);
    for (int c = 0; c < dnaLen; c++) {
        char* b = charSeq->get(charSeq->elements, c);
        CuAssertStrEquals(testCase, &testSequence[c], b);
    }
}

static void test_kmers(CuTest* testCase) {
    char* testSequence = "GATACAGATACAGATACA";
    int len = strlen(testSequence);
    char testKmers[14][6] = {"GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                             "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                             "GATAC", "ATACA"};
    Sequence* kmerSeq = sequenceConstruct(len, testSequence, getKmer);
    for (int x = 0; x < 5; x++) {
        char* kmer = kmerSeq->get(kmerSeq->elements, x);
        CuAssertStrEquals(testCase, &testKmers[x], kmer);
    }
}

static void test_events(CuTest* testCase) {
    char testKmers[12][6] = { "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA",
                              "GATAC", "ATACA", "TACAG", "ACAGA", "CAGAT", "AGATA" };
    double testMeans[] = { 12.12, 23.23, 34.34, 45.45, 56.56, 67.67,
                           12.12, 23.23, 34.34, 45.45, 56.56, 67.67 };
    
    void* eventSequence = eventSequenceConstruct(12, testMeans, *testKmers);
    Sequence* eventSeq = sequenceConstruct(12, eventSequence, getEvent);
    
    for (int x = 0; x < 12; x++) {
        Event *ev = eventSeq->get(eventSeq->elements, x);
        // test for kmer
        CuAssertStrEquals(testCase, &testKmers[x], ev->kmer);
        CuAssertDblEquals(testCase, testMeans[x], ev->mean, 0.0);
    }
}

CuSuite* shimTestSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_chars);
    SUITE_ADD_TEST(suite, test_kmers);
    SUITE_ADD_TEST(suite, test_events);

    return suite;
}

/*
{
    printf("Quick Tests\n\n");
    char testKmers[3][6] = {"ATGAC", "TGACA", "GACAT"};
    double testMeans[] = {45.45, 65.65, 101.101};
    char testDNA[] = "GATACAGATACA";
    
    // Test chars
    int dnaLen = strlen(testDNA);
    Sequence* charSeq = sequenceConstruct(dnaLen, testDNA, getBase);
    for (int c = 0; c < dnaLen; c++) {
        char* b = charSeq->get(charSeq->elements, c);
        printf("%c is the base at position %d\n", *b, c);
    }
    printf("\n");
    // Test kmers
    Sequence* kmerSeq = sequenceConstruct(dnaLen, testDNA, getKmer);
    
    for (int k = 0; k < (dnaLen - 4); k++) {
        char* kmer = kmerSeq->get(kmerSeq->elements, k);
        printf("%s is the kmer at position %d\n", kmer, k);
        
    }
    printf("\n");
    
    
    // Test Events
    void* eventSequence = eventSequenceConstruct(3, testMeans, *testKmers);
    Sequence* eventSeq = sequenceConstruct(3, eventSequence, getEvent);
    for (int x = 0; x < 3; x++) {
//        Event *ev;
        Event *ev = eventSeq->get(eventSeq->elements, x);
        printf("Event %d retrieved: kmer: %s mean: %f base: %c \n", x, 
                                                                    ev->kmer, 
                                                                    ev->mean, 
                                                                    ev->base);
        
    }
        
    return 0;
}
*/