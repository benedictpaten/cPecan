/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "realigner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randomSequences.h"

static void test_poa_getReferenceGraph(CuTest *testCase) {

}

static void test_poa_augment(CuTest *testCase) {

}

static void test_poa_realign(CuTest *testCase) {

}

CuSuite* realignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
    SUITE_ADD_TEST(suite, test_poa_augment);
    SUITE_ADD_TEST(suite, test_poa_realign);

    return suite;
}
