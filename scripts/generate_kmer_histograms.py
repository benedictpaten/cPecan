#!/usr/bin/env python

import sys
from alignmentAnalysisLib import Kmer_histogram
from itertools import product
from argparse import ArgumentParser
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--alignments', '-a', action='store',
                        dest='alns', required=False, type=str, default=None,
                        help="alignment files")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    parser.add_argument('--kmer', '-k', action='append', default=None,
                        dest='kmers', required=False, type=str,
                        help="which kmers to plot")
    parser.add_argument('--strand', '-s', action='store',
                        dest='strand', required=False, type=str, default=None,
                        help="strand")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.25, dest='threshold')
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()

def histogram_runner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            k = Kmer_histogram(**f)
            k.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def main(args):
    # get the kmers we want
    args = parse_args()
    kmers_of_interest = []
    if args.kmers is None:
        for kmer in product("ACTG", repeat=6):
            kmer = ''.join(kmer)
            if "C" in kmer:
                kmers_of_interest.append(kmer)
            else:
                continue
    else:
        kmers_of_interest = args.kemrs

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    for kmer in kmers_of_interest:
        hist_args = {
            "path_to_alignments": args.alns,
            "kmer": kmer,
            "strand": args.strand,
            "threshold": args.threshold,
            "max_assignments": 1000,
            "out_dir": args.out,
        }
        work_queue.put(hist_args)

    for w in xrange(workers):
        p = Process(target=histogram_runner, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')

if __name__ == "__main__":
    sys.exit(main(sys.argv))
