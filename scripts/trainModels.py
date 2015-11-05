#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""
from __future__ import print_function
import h5py as h5
import subprocess
import os
import sys
sys.path.append("../")
from serviceCourse.parsers import read_fasta
from serviceCourse.sequenceTools import reverse_complement
from nanoporeLib import *
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser (description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str,
                        help="directory with fast5 files to train on, should have trailing /")

    parser.add_argument('--train_amount', '-m', action='store', dest='amount',
                        required=False, type=int, default='1000000')

    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str)

    parser.add_argument('--index', '-i', action='store', dest='bwa_index',
                        required=False, type=str, default=None)

    parser.add_argument('--inputHmm', '-y', action='store', dest='inHmm',
                        required=False, type=str, default=None)

    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the trained model, should have trailing slash")

    parser.add_argument('--strand', '-s', action='store', dest='strand',
                        required=True, type=str)

    parser.add_argument('--iterations', '-t', action='store', dest='iter',
                        required=False, type=str, default='50')

    parser.add_argument('--strawMan', '-sm', action='store_true', dest='strawMan',
                        required=False, default=False)

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


def cull_training_files(directory, training_amount):
    print("trainModels - culling training files.  ", end="", file=sys.stderr)
    fast5s = [x for x in os.listdir(directory) if x.endswith(".fast5")]
    shuffle(fast5s)
    training_files = []
    total_amount = 0

    for i in xrange(len(fast5s)):
        training_files.append(directory+fast5s[i])
        total_amount += get_2d_length(directory+fast5s[i])
        if total_amount >= training_amount:
            break
    print("Culled {nb_files} training files.".format(nb_files=len(training_files)),
          end="\n", file=sys.stderr)

    return training_files


def doEM(in_fast5, reference, destination, strawMan_flag, strand="Template", bwa_index=None, in_hmm=None, iterations=50):
    """
    :param in_fast5: path to the fast5 file
    :param reference: path to the reference sequence (fasta)
    :param strand, string, "Template" or "Complement"
    :param strawMan_flag, set to True to use straw man model
    :param bwa_index: path to index files (no suffix)
    :param in_hmm: input hmm to continue training
    :param destination: place to put the trained hmm
    :param iterations: number of iterations to do
    :return: void, makes a trained hmm (or updates the existing one)
    """
    # Preamble: first we align the 2D read to the reference and get the orientation, ie does
    # it align to the reference as it is in the fastA or to the reverse complement

    # containers and defaults
    temp_dir = destination + "tempFiles_{strand}/".format(strand=strand)
    # make the temp directory, if needed
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))
    temp_np_read = temp_dir + "temp_nanoporeRead.npRead"
    temp_2d_read = temp_dir + "temp_2d_read.fa"
    make_npRead_and_2d_seq(in_fast5, temp_np_read, temp_2d_read)

    # if there is no bwa index given, make one
    if bwa_index is None:
        bwa = Bwa(reference)
        bwa.build_index(temp_dir)
        bwa_ref_index = temp_dir + "temp_bwaIndex"
    else:
        bwa_ref_index = bwa_index

    # align with bwa
    #bwa_dir = "/Users/Rand/projects/BGCs/submodules/bwa/"  # todo require bwa in path remove this
    #command = "{bwaDir}bwa mem -x ont2d {index} {query}".format(bwaDir=bwa_dir, index=bwa_ref_index,
    #                                                            query=temp_2d_read)
    # this is a small SAM file that comes from bwa
    #print("trainModels - mapping read from {inFile}".format(inFile=in_fast5), file=sys.stderr)
    #aln = subprocess.check_output(command.split())
    #aln = aln.split("\t") # split

    # migrated function outside
    orientation = orient_read_with_bwa(bwa_index=bwa_ref_index, query=temp_2d_read)

    # forward strand
    if orientation == 0:
        forward = True

    # backward strand
    if orientation == 16:
        forward = False

    # didn't map
    elif (orientation != 0) and (orientation != 16):
        print("\n\ntrainModels - read didn't map", file=sys.stderr)
        return  # todo double check if this works correctly

    # EM training routine: now we can run the training, we run training on either the template or
    # complement, so that we can run this program in parallel and really do both at once

    # containers and defaults
    temp_ref_seq = temp_dir + "temp_ref_seq.txt"
    path_to_vanillaAlign = "./vanillaAlign"
    strand_flags = ["-y", "-t"]
    training_hmm = destination + "{strand}_trained.hmm".format(strand=strand)

    # make the temp sequence file
    make_temp_sequence(reference, forward, temp_ref_seq)

    # switch flags if we're doing the complement
    if strand == "Complement":
        strand_flags = ["-z", "-c"]

    use_strawMan_model = ""
    if strawMan_flag is True:
        use_strawMan_model = "--s "

    print("trainModels - starting B-W on file: {inFile}".format(inFile=in_fast5), end="\n", file=sys.stderr)
    # training commands
    em_command_start = "{vanillaAlign} {straw}-r {ref} -q {npRead} {outHmmFlag} {outHmm} -i {iter}".format(
        vanillaAlign=path_to_vanillaAlign, straw=use_strawMan_model, ref=temp_ref_seq, npRead=temp_np_read,
        outHmmFlag=strand_flags[1], outHmm=training_hmm, iter=iterations)
    em_command_continue = \
        "{vanillaAlign} {straw}-r {ref} -q {npRead} {inHmmFlag} {inHmm} {outHmmFlag} {outHmm} -i {iter}"\
        .format(vanillaAlign=path_to_vanillaAlign, straw=use_strawMan_model, inHmmFlag=strand_flags[0], inHmm=in_hmm,
                ref=temp_ref_seq, npRead=temp_np_read, outHmmFlag=strand_flags[1],
                outHmm=training_hmm, iter=iterations)

    # if we are given an input HMM
    if in_hmm is None:
        print("trainModels - running command", em_command_start, end="\n", file=sys.stderr)
        os.system(em_command_start)

    # if we're starting from scratch
    else:
        print("trainModels - running command", em_command_continue, end="\n", file=sys.stderr)
        os.system(em_command_continue)


def main(args):
    # parse command line arguments
    args = parse_args()

    start_message = """\n
    Starting Baum-Welch training.
    Directory with training files: {files_dir}
    Training on {amount} bases.
    Using reference sequence: {ref}
    Input hmm: {inHmm}
    Writing trained hmm to: {outLoc}
    Training on strand: {strand}
    Performing {iterations} iterations.
    Using strawMan model: {straw}
    \n
    """.format(files_dir=args.files_dir, amount=args.amount, ref=args.ref, inHmm=args.inHmm,
               outLoc=args.out, strand=args.strand, iterations=args.iter, straw=args.strawMan)
    print(start_message, file=sys.stderr)

    # make directory to put temporary files
    temp_dir = args.out + "tempFiles_{strand}/".format(strand=args.strand)
    # make the temp directory, if needed
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))

    # index the reference for bwa, if needed
    bwa_ref_index = ''
    if args.bwa_index is None:
        print("trainModels - indexing reference", file=sys.stderr)
        bwa_ref_index = get_bwa_index(args.ref, temp_dir)
        print("trainModels - indexing reference, done\n", file=sys.stderr)
    else:
        bwa_ref_index = args.bwa_index

    # get a random list of files containing the number of bases we want to train
    training_file_list = cull_training_files(args.files_dir, args.amount)

    training_hmm = args.out + "{strand}_trained.hmm".format(strand=args.strand)

    # get started. if there is no input training hmm, then we start from scratch, and train on one file
    if args.inHmm is None:
        get_started_file = training_file_list.pop()
        doEM(in_fast5=get_started_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=None, destination=args.out, iterations=args.iter, strawMan_flag=args.strawMan)

    # otherwise train with the input hmm on one file
    if args.inHmm is not None:
        get_started_file = training_file_list.pop()
        doEM(in_fast5=get_started_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=args.inHmm, destination=args.out, iterations=args.iter, strawMan_flag=args.strawMan)

    # train on the rest of the files
    for training_file in training_file_list:
        doEM(in_fast5=training_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=training_hmm, destination=args.out, iterations=args.iter, strawMan_flag=args.strawMan)

    print("\nFinished Training routine, exiting.\n", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

