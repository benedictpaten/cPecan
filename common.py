#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs
"""

import os
import random
import sys

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getLogLevelString

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString

def runCPecanEm(sequenceFiles, alignmentsFile, outputModelFile,
                 inputModelFile=None, 
                 modelType=None,
                 jobTreeDir=None,
                 iterations=None, randomStart=None, 
                 trials=None,
                 optionsToRealign=None,
                 logLevel=None, 
                 updateTheBand=None,
                 maxAlignmentLengthPerJob=None,
                 maxAlignmentLengthToSample=None,
                 useDefaultModelAsStart=None, 
                 setJukesCantorStartingEmissions=None,
                 trainEmissions=None,
                 tieEmissions=None,
                 outputTrialHmms = None,
                 outputXMLModelFile = None,
                 blastScoringMatrixFile=None):
    logLevel = getLogLevelString2(logLevel)
    jobTreeDir= nameValue("jobTree", jobTreeDir, str)
    inputModelFile= nameValue("inputModel", inputModelFile, str)
    modelType = nameValue("modelType", modelType, str)
    iterations = nameValue("iterations", iterations, int)
    trials = nameValue("trials", trials, int)
    randomStart = nameValue("randomStart", randomStart, bool)
    updateTheBand = nameValue("updateTheBand", updateTheBand, bool)
    maxAlignmentLengthPerJob = nameValue("maxAlignmentLengthPerJob", maxAlignmentLengthPerJob, int)
    maxAlignmentLengthToSample = nameValue("maxAlignmentLengthToSample", maxAlignmentLengthToSample, int)
    optionsToRealign = nameValue("optionsToRealign", optionsToRealign, quotes=True)
    useDefaultModelAsStart = nameValue("useDefaultModelAsStart", useDefaultModelAsStart, bool) 
    trainEmissions = nameValue("trainEmissions", trainEmissions, bool)
    tieEmissions = nameValue("tieEmissions", tieEmissions, bool)
    setJukesCantorStartingEmissions = nameValue("setJukesCantorStartingEmissions", setJukesCantorStartingEmissions, float)
    outputTrialHmms = nameValue("outputTrialHmms", outputTrialHmms, bool)
    outputXMLModelFile = nameValue("outputXMLModelFile", outputXMLModelFile, str)
    blastScoringMatrixFile = nameValue("blastScoringMatrixFile", blastScoringMatrixFile, str)
    
    system("cPecanEm --sequences '%s' --alignments %s --outputModel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
           (" ".join(sequenceFiles), alignmentsFile, outputModelFile, iterations, trials, randomStart, 
            jobTreeDir, inputModelFile, optionsToRealign, modelType,
            maxAlignmentLengthPerJob, maxAlignmentLengthToSample, updateTheBand, useDefaultModelAsStart, 
            trainEmissions, tieEmissions, setJukesCantorStartingEmissions, outputTrialHmms, 
            outputXMLModelFile, blastScoringMatrixFile))
