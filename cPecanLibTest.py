#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
from sonLib.bioio import parseSuiteTestOptions, getLogLevelString, logger, system

"""Tests cPecanLib. 
"""

class TestCase(unittest.TestCase):
    def testCPecanLib(self):
        """Run all the cPecanLib CuTests, fail if any of them fail.
        """
        system("cPecanLibTests %s" % getLogLevelString())

def main():
    parseSuiteTestOptions()
    unittest.main()
        
if __name__ == '__main__':
    main()
