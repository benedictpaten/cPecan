#!/usr/bin/env python
"""Small library for working with MinION data
"""

from __future__ import print_function
from numpy import log2
from itertools import islice
from serviceCourse.parsers import read_fastq
import h5py


class NanoporeRead(object):
    def __init__(self, fast_five_file):
        self.fastFive = h5py.File(fast_five_file, 'r')
        self.template_event_table = self.fastFive['/Analyses/Basecall_2D_000/BaseCalled_template/Events']
        self.complement_event_table = self.fastFive['/Analyses/Basecall_2D_000/BaseCalled_complement/Events']
        self.template_scale = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["scale"]
        self.template_shift = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["shift"]
        self.template_drift = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["drift"]
        self.template_var = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["var"]
        self.template_scale_sd = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["scale_sd"]
        self.template_var_sd = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["var_sd"]
        self.complement_scale = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["scale"]
        self.complement_shift = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["shift"]
        self.complement_drift = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["drift"]
        self.complement_var = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["var"]
        self.complement_scale_sd = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["scale_sd"]
        self.complement_var_sd = self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["var_sd"]
        self.template_map = []
        self.complement_map = []
        self.template_read_sequence = ''
        self.complement_read_sequence = ''
        self.template_read_length = 0
        self.complement_read_length = 0
        self.get_event_map()
        self.get_read_sequence()

    def get_event_map(self):
        def make_map(events):
            event_map = [0]
            previous_prob = 0
            for i, line in islice(enumerate(events), 1, None):
                move = line[6]
                this_prob = line[8]
                if move == 1:
                    event_map.append(i)
                if move > 1:
                    for skip in xrange(move - 1):
                        event_map.append(i - 1)
                    event_map.append(i)
                if move == 0:
                    if this_prob > previous_prob:
                        event_map[-1] = i
                previous_prob = this_prob
            final_event_index = [event_map[-1]]
            padding = final_event_index * 5 # make this a kmer-measured thing
            event_map = event_map + padding
            return event_map
        self.template_map = make_map(self.template_event_table)
        self.complement_map = make_map(self.complement_event_table)
        return

    def event_generator(self, event_table):
        for line in event_table:
            yield line[0], line[2], line[3] # mean, stdev, length

    def get_read_sequence(self):
        template_path_to_fastq = "/Analyses/Basecall_2D_000/BaseCalled_template/Fastq"
        complement_path_to_fastq = "/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq"
        self.template_read_sequence = self.fastFive[template_path_to_fastq][()].split()[2]
        self.template_read_length = len(self.template_read_sequence)
        self.complement_read_sequence = self.fastFive[complement_path_to_fastq][()].split()[2]
        self.complement_read_length = len(self.complement_read_sequence)
        return




