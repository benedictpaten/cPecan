#!/usr/bin/env python
"""Convert a MinION fast5 to a npRead for cPecan

Format:
line 1 [2D read length] [# of template events] [# of complement events]
       [template scale] [template shift] [template var] [template scale_sd] [template var_sd]
       [complement scale] [complement shift] [complement var] [complement scale_sd] [complement var_sd] \n
line 2 [2D read sequence] \n
line 3 [template event map] \n
line 4 [template events (mean, stddev, length)] \n
line 5 [complement event map] \n
line 6 [complement events (mean, stddev, length)] \n
"""

from __future__ import print_function
from nanoporeLib import NanoporeRead
import sys


def main():
    if len(sys.argv) < 3:
        print("USAGE: python thisScript.py /path/to/fast5")

    if len(sys.argv) == 3:
        # setup
        out_file = open(sys.argv[2], 'w')
        # load and transform
        npRead = NanoporeRead(sys.argv[1])
        npRead.transform_events(npRead.template_events, npRead.template_drift)
        npRead.transform_events(npRead.complement_events, npRead.complement_drift)
        # output
        # line 1
        print(len(npRead.twoD_read_sequence), end=' ', file=out_file) # 2D read length
        print(len(npRead.template_events), end=' ', file=out_file)    # nb of template events
        print(len(npRead.complement_events), end=' ', file=out_file)  # nb of complement events
        print(npRead.template_scale, end=' ', file=out_file)          # template scale
        print(npRead.template_shift, end=' ', file=out_file)          # template shift
        print(npRead.template_var, end=' ', file=out_file)            # template var
        print(npRead.template_scale_sd, end=' ', file=out_file)       # template scale_sd
        print(npRead.template_var_sd, end=' ', file=out_file)         # template var_sd
        print(npRead.complement_scale, end=' ', file=out_file)        # complement scale
        print(npRead.complement_shift, end=' ', file=out_file)        # complement shift
        print(npRead.complement_var, end=' ', file=out_file)         # complement var
        print(npRead.complement_scale_sd, end=' ', file=out_file)     # complement scale_sd
        print(npRead.complement_var_sd, end='\n', file=out_file)      # complement var_sd

        # line 2
        print(npRead.twoD_read_sequence, end='\n', file=out_file)

        # line 3
        for _ in npRead.template_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 4
        for mean, start, stdev, length in npRead.template_events:
            print(mean, stdev, length, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 5
        for _ in npRead.complement_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 6
        for mean, start, stdev, length in npRead.complement_events:
            print(mean, stdev, length, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)


if __name__ == '__main__':
    main()
