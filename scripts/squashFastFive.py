#!/usr/bin/env python
"""Keeping around for comparison of extracted nanopore models
"""

from __future__ import print_function
from numpy import log2
import sys
import h5py

# eventually I want this program to extract the template model, the complement
# model and the string of level measurements. Probably a good idea to start
# thinking about what format to use...

if len(sys.argv) < 2:
    print("USAGE: python squashFastFive.py /path/to/fast5")

if len(sys.argv) == 2:
    # get input fast5
    path_to_fastFive = sys.argv[1]
    # make paths to models
    template_file_path = "../models/template.eTable.model"
    complement_file_path = "../models/complement.eTable.model"
    file_paths = [template_file_path, complement_file_path]
    #template_file = open(template_file_path, 'w')
    #complement_file = open(complement_file_path, 'w')
    # open fast5
    fastFive = h5py.File(path_to_fastFive, 'r')
    template_model = fastFive['/Analyses/Basecall_2D_000/BaseCalled_template/Model']
    complement_model = fastFive['/Analyses/Basecall_2D_000/BaseCalled_complement/Model']
    stay_prob = fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["stay_prob"]
    stay_prob = log2(stay_prob)
    skip_prob = fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["skip_prob"]
    skip_prob = log2(skip_prob)

    for i, model in enumerate([template_model, complement_model]):
        out_file = open(file_paths[i], 'w')
        nb_params = 0
        print("0", end=' ', file=out_file) # placeholder for covariance parameter
        for kmer, level_mean, level_stdev, sd_mean, sd_stdev, weight in model:
            print(level_mean, level_stdev, sd_mean, sd_stdev, end=' ', file=out_file)
            nb_params += 1
        print("", end="\n", file=out_file)
        for _ in xrange(nb_params):
            print(skip_prob, end=' ', file=out_file)
        print("", end="\n", file=out_file)
        for _ in xrange(nb_params):
            print(stay_prob, end=' ', file=out_file)
        print("", end="\n", file=out_file)