#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""
from __future__ import print_function, division
import subprocess
import os
import h5py as h5
import numpy as np
import sys
sys.path.append("../")
from multiprocessing import Process, Queue, current_process, Manager
from nanoporeLib import *
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser (description=__doc__)

    parser.add_argument('--file_directory', '-d', action='append',
                        dest='files_dir', required=True, type=str,
                        help="directories with fast5 files to train on")

    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="location of refrerence sequence in FASTA")

    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the trained model, and use for working directory.")

    parser.add_argument('--iterations', '-t', action='store', dest='iter',
                        required=True, type=int)

    parser.add_argument('--train_amount', '-m', action='store', dest='amount',
                        required=True, type=int,
                        help="limit the total length of sequence to use in training.")

    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")

    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")

    parser.add_argument('--banded', '-b', action='store_true', dest='banded',
                        default=False, help='flag, use banded alignment heuristic')

    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=True,
                        type=int, help="number of jobs to run concurrently")

    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        required=True, help="decide which model to use, vanilla by default")

    args = parser.parse_args()
    return args


def get_2d_length(fast5):
    read = h5.File(fast5, 'r')
    read_length = 0
    twoD_read_sequence_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
    if not (twoD_read_sequence_address in read):
        print("This read didn't have a 2D read", fast5, end='\n', file=sys.stderr)
        read.close()
        return 0
    else:
        read_length = len(read[twoD_read_sequence_address][()].split()[2])
        read.close()
        return read_length


def cull_training_files(directories, training_amount):
    print("trainModels - culling training files.\n", end="", file=sys.stderr)

    training_files = []

    for directory in directories:
        fast5s = [x for x in os.listdir(directory) if x.endswith(".fast5")]
        shuffle(fast5s)

        total_amount = 0
        n = 0
        for i in xrange(len(fast5s)):
            training_files.append(directory + fast5s[i])
            n += 1
            total_amount += get_2d_length(directory + fast5s[i])
            if total_amount >= training_amount:
                break
        print("Culled {nb_files} training files from {dir}.".format(nb_files=n,
                                                                    dir=directory),
              end="\n", file=sys.stderr)

    return training_files


def get_expectations(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run(get_expectations=True)
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def get_model(type, symbol_set_size):
    if type == "threeState":
        return ContinuousPairHmm(model_type=type, symbol_set_size=symbol_set_size)
    if type == "vanilla":
        return ConditionalSignalHmm(model_type=type, symbol_set_size=symbol_set_size)


def add_and_norm_expectations(path, files, model, hmm_file):
    model.likelihood = 0
    for f in files:
        model.add_expectations_file(path + f)
        #os.remove(path + f)  # try reactivating this
    model.normalize()
    model.write(hmm_file)
    model.running_likelihoods.append(model.likelihood)


def main(argv):
    # parse command line arguments
    args = parse_args()

    start_message = """\n
    Starting Baum-Welch training.
    Directories with training files: {files_dir}
    Training on {amount} bases.
    Using reference sequence: {ref}
    Input template/complement models: {inTHmm}/{inCHmm}
    Writing trained models to: {outLoc}
    Performing {iterations} iterations.
    Using model: {model}
    \n
    """.format(files_dir=args.files_dir, amount=args.amount, ref=args.ref,
               inTHmm=args.in_T_Hmm, inCHmm=args.in_C_Hmm, outLoc=args.out,
               iterations=args.iter, model=args.stateMachineType)
    print(start_message, file=sys.stdout)

    if not os.path.isfile(args.ref):  # TODO make this is_fasta(args.ref)
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make directory to put the files we're using files
    working_folder = FolderHandler()
    working_directory_path = working_folder.open_folder(args.out + "tempFiles_expectations")

    # index the reference for bwa
    print("signalAlign - indexing reference", file=sys.stderr)
    bwa_ref_index = get_bwa_index(args.ref, working_directory_path)
    print("signalAlign - indexing reference, done", file=sys.stderr)

    # make model objects, these handle normalizing, loading, and writing
    template_model = get_model(type=args.stateMachineType, symbol_set_size=4096)
    complement_model = get_model(type=args.stateMachineType, symbol_set_size=4096)

    # make some paths to files to hold the HMMs
    template_hmm = working_folder.add_file_path("template_trained.hmm")
    complement_hmm = working_folder.add_file_path("complement_trained.hmm")

    print("Starting {iterations} iterations.\nRunning likelihoods\nTempalte\tComplement".format(
        iterations=args.iter), file=sys.stdout)

    for i in xrange(args.iter):
        # if we're starting there are no HMMs
        if i == 0:
            in_template_hmm = None
            in_complement_hmm = None
        else:
            in_template_hmm = template_hmm
            in_complement_hmm = complement_hmm

        # first cull a set of files to get expectations on
        training_files = cull_training_files(args.files_dir, args.amount)

        # setup
        workers = args.nb_jobs
        work_queue = Manager().Queue()
        done_queue = Manager().Queue()
        jobs = []

        # get expectations for all the files in the queue
        for fast5 in training_files:
            alignment_args = {
                "reference": args.ref,
                "destination": working_directory_path,
                "stateMachineType": args.stateMachineType,
                "bwa_index": bwa_ref_index,
                "in_templateHmm": in_template_hmm,
                "in_complementHmm": in_complement_hmm,
                "banded": args.banded,
                "in_fast5": fast5
            }
            work_queue.put(alignment_args)

        for w in xrange(workers):
            p = Process(target=get_expectations, args=(work_queue, done_queue))
            p.start()
            jobs.append(p)
            work_queue.put('STOP')

        for p in jobs:
            p.join()

        done_queue.put('STOP')

        # load then normalize the expectations
        template_expectations_files = [x for x in os.listdir(working_directory_path)
                                       if x.endswith(".template.expectations")]

        complement_expectations_files = [x for x in os.listdir(working_directory_path)
                                         if x.endswith(".complement.expectations")]

        if len(template_expectations_files) > 0:
            add_and_norm_expectations(path=working_directory_path,
                                      files=template_expectations_files,
                                      model=template_model,
                                      hmm_file=template_hmm)

        if len(complement_expectations_files) > 0:
            add_and_norm_expectations(path=working_directory_path,
                                      files=complement_expectations_files,
                                      model=complement_model,
                                      hmm_file=complement_hmm)

        if len(template_model.running_likelihoods) > 0 and len(complement_model.running_likelihoods) > 0:
            print("{t_likelihood}\t{c_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1],
                                                          c_likelihood=complement_model.running_likelihoods[-1]))

    print("signalAlign - finished training routine", file=sys.stdout)
    print("signalAlign - finished training routine", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

