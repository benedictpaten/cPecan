#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cPecan.cPecanLibTest import TestCase as cPecanLibTest
from cPecan.cPecanRealignTest import TestCase as realignTest
from cPecan.cPecanEmTest import TestCase as eMTest

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(cPecanLibTest, 'test'),
                                   unittest.makeSuite(realignTest, 'test'),
                                   unittest.makeSuite(eMTest, 'test')))
    return allTests
        
def main():
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
                
