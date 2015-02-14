#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cPecan.cPecanLibTest import TestCase as cPecanLibTest
from cPecan.cactus_realignTest import TestCase as realignTest
from cPecan.cactus_expectationMaximisationTest import TestCase as expectationMaximisationTest

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(cPecanLibTest, 'test'),
                                   unittest.makeSuite(realignTest, 'test'),
                                   unittest.makeSuite(expectationMaximisationTest, 'test'),                          
                                   progressiveSuite()))
    return allTests
        
def main():
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
                
