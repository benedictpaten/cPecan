#!/usr/bin/env python
"""I will later put this all together into something useful
"""

from __future__ import print_function
from nanoporeLib import NanoporeRead
import sys
import h5py


def main():
    if len(sys.argv) < 3:
        print("USAGE: python thisScript.py /path/to/fast5")

    if len(sys.argv) == 3:
        out_file = open(sys.argv[2], 'w')
        npRead = NanoporeRead(sys.argv[1])
        print(npRead.template_read_length, end=' ', file=out_file)
        print(len(npRead.template_event_table), end=' ', file=out_file)
        print(npRead.template_scale, end=' ', file=out_file)
        print(npRead.template_shift, end='\n', file=out_file)
        print(npRead.template_read_sequence, end='\n', file=out_file)
        for _ in npRead.template_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)
        for mean, stdev, length in npRead.event_generator(npRead.template_event_table):
            print(mean, stdev, length, end=' ', file=out_file)


if __name__ == '__main__':
    main()
