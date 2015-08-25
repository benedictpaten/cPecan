/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sonLib.h"
#include "../../sonLib/lib/CuTest.h"
#include "../../sonLib/lib/sonLibCommon.h"

CuSuite* pairwiseAlignmentTestSuite(void);
CuSuite* multipleAlignerTestSuite(void);
CuSuite* pairwiseAlignmentLongTestSuite(void);
//CuSuite* shimTestSuite(void);
//CuSuite* kmerTestSuite(void);

int stBaseAlignerRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, pairwiseAlignmentTestSuite());
    //CuSuiteAddSuite(suite, multipleAlignerTestSuite());
    //CuSuiteAddSuite(suite, pairwiseAlignmentLongTestSuite());
    //CuSuiteAddSuite(suite, shimTestSuite());
    //CuSuiteAddSuite(suite, kmerTestSuite());
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
  return stBaseAlignerRunAllTests();
}
