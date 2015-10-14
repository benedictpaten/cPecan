#!/usr/bin/env python
"""Investigating the duration of events
"""

from __future__ import print_function
from nanoporeLib import NanoporeRead
from glob import glob
from math import floor
from random import sample
import sys


def write_event_durations_to_file(directories, destination):
    # make the files to hold durations
    template_file = open(destination + "template_durations.csv", 'w')
    complement_file = open(destination + "complement_durations.csv", 'w')
    all_together_file = open(destination + "all_durations.csv", 'w')

    # go through each directory and gather up durations
    for directory in directories:
        directory_files = directory + "*.fast5"
        nb_files = len(glob(directory_files))
        sample_set = int(floor(0.25 * nb_files))
        print("Found {nb_files} in folder {folder} sampling {nb_sample_set}".format(nb_files=nb_files,
                                                                                    folder=directory,
                                                                                    nb_sample_set=sample_set),
              file=sys.stderr)
        # make this random.sample, sample 10% of the files in the directory
        for f in sample(glob(directory_files), sample_set):
            npRead = NanoporeRead(f)
            if hasattr(npRead, 'template_event_table') and hasattr(npRead, 'complement_event_table'):
                file_name = f.split('/')[-1]
                print(file_name, file=sys.stderr)
                for t, c in map(None, npRead.template_event_table, npRead.complement_event_table):
                    if t is not None:
                        print(t[3], t[6], 't', sep=',', end='\n', file=template_file)
                        print(t[3], t[6], 't',sep=',', end='\n', file=all_together_file)
                    if c is not None:
                        print(c[3], c[6], 'c', sep=',', end='\n', file=complement_file)
                        print(c[3], c[6], 'c', sep=',', end='\n', file=all_together_file)


def main():
    #print(sys.argv[1:-1])  # gets all arguments to the script except the last one
    #print(sys.argv[-1])  # gets the last argument
    files = sys.argv[1:-1]
    write_event_durations_to_file(files, sys.argv[-1])


if __name__ == '__main__':
    main()
