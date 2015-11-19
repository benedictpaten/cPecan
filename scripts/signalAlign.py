#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import sys
sys.path.append("../")
from nanoporeLib import *
from multiprocessing import Process, Queue, current_process, Manager
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with fast5s for alignment")
    parser.add_argument('-nb_files', '-nb', action='store', dest='nb_files', required=False,
                        default=50, type=int, help="maximum number of reads to align")

    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")

    # reference
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str, help="reference sequence to align to, in FASTA")
    parser.add_argument('--index', '-i', action='store', dest='bwa_index',
                        required=False, type=str, default=None,
                        help="path to bwa index files, if you've already made them")
    # input HMMs
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")

    # output
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        required=False, default=None, help="decide which model to use, vanilla by default")

    args = parser.parse_args()
    return args


def aligner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.do_alignment()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def main(args):
    # parse args
    args = parse_args()

    start_message = """
    Starting Signal Align
    Aligning files from: {fileDir}
    Aligning to reference: {reference}
    Aligning {nbFiles}
    Input template HMM: {inThmm}
    Input complement HMM: {inChmm}
    """.format(fileDir=args.files_dir, reference=args.ref, nbFiles=args.nb_files,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm)

    print(start_message, file=sys.stderr)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "tempFiles_alignment")

    # index the reference for bwa, if needed
    if args.bwa_index is None:
        print("signalAlign - indexing reference", file=sys.stderr)
        bwa_ref_index = get_bwa_index(args.ref, temp_dir_path)
        print("signalAlign - indexing reference, done", file=sys.stderr)
    else:
        bwa_ref_index = args.bwa_index

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    fast5s = [x for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    nb_files = args.nb_files
    if nb_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:nb_files]

    for fast5 in fast5s:
        alignment_args = {
            "reference": args.ref,
            "destination": temp_dir_path,
            "stateMachineType": args.stateMachineType,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
            "in_fast5": args.files_dir + fast5
        }
        work_queue.put(alignment_args)

    for w in xrange(workers):
        p = Process(target=aligner, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')
    print("signalAlign - finished alignments\n", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
