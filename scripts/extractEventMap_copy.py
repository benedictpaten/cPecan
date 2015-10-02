#!/usr/bin/env python
"""I will later put this all together into something useful
"""

from __future__ import print_function
from itertools import islice
import sys
import h5py

# eventually I want this program to extract the template model, the complement
# model and the string of level measurements. Probably a good idea to start
# thinking about what format to use...

if len(sys.argv) < 2:
    print("USAGE: python thisScript.py /path/to/fast5")

if len(sys.argv) == 2:
    # get input fast5
    path_to_fastFive = sys.argv[1]
    # make paths to models
    template_file_path = "../models/template.eTable.model"
    complement_file_path = "../models/complement.eTable.model"
    file_paths = [template_file_path, complement_file_path]

    # open fast5
    fastFive = h5py.File(path_to_fastFive, 'r')

    template_events = fastFive['/Analyses/Basecall_2D_000/BaseCalled_template/Events']
    complement_events = fastFive['/Analyses/Basecall_2D_000/BaseCalled_complement/Events']

    template_map = [0]
    previous_prob = 0
    for i, line in islice(enumerate(template_events), 1, None):
        index = i
        move = line[6]
        this_prob = line[8]
        if move == 1:
            template_map.append(i)
        if move > 1:
            for skip in xrange(move - 1):
                template_map.append(i - 1)
            template_map.append(i)
        if move == 0:
            if this_prob > previous_prob:
                template_map[-1] = i
        previous_prob = this_prob
    final_event_index = [template_map[-1]]
    padding = final_event_index * 5 # make this a kmer-measured thing
    template_map = template_map + padding


print(template_map, len(template_map))