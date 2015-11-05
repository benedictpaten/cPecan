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
from nanoporeLib import NanoporeRead, NanoporeModel
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

    args = parser.parse_args()
    return args


def make_temp_npRead(fast5, npRead_dest, twod_read_dest):
    """process a MinION .fast5 file into a npRead file for use with signalAlign also extracts
    the 2D read into fasta format
    """
    # setup
    out_file = open(npRead_dest, 'w')
    temp_fasta = open(twod_read_dest, "w")

    # load and transform
    npRead = NanoporeRead(fast5)
    npRead.get_2D_event_map()
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
    print(npRead.complement_var, end=' ', file=out_file)          # complement var
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

    # make temp read
    npRead.extract_2d_read(temp_fasta)
    npRead.close()
    return


def make_temp_sequence(fasta, forward, destination):
    """extract the sequence from a fasta and put into a simple file that is used by signalAlign
    """
    out_file = open(destination, "w")
    for header, comment, sequence in read_fasta(fasta):
        if forward is False:
            sequence = reverse_complement(sequence)
        print(sequence, end='\n', file=out_file)


class Bwa(object):
    """run BWA easily
    """
    def __init__(self, target):
        self.target = target
        self.bwa_dir = "/Users/Rand/projects/BGCs/submodules/bwa/"
        self.db_handle = ''

    def build_index(self, destination):
        # make a place to put the database
        path_to_bwa_index = destination

        # build database
        self.db_handle = path_to_bwa_index + '/temp_bwaIndex'
        os.system("{0}bwa index -p {1} {2}".format(self.bwa_dir, self.db_handle, self.target))

    def run(self, query):
        # run alignment
        os.system("{0}bwa mem -x ont2d {1} {2}".format(self.bwa_dir, self.db_handle, query))


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
    fast5s = [x for x in os.listdir(directory) if x.endswith(".fast5")]
    shuffle(fast5s)
    training_files = []
    total_amount = 0

    for i in xrange(len(fast5s)):
        training_files.append(directory+fast5s[i])
        total_amount += get_2d_length(directory+fast5s[i])
        if total_amount >= training_amount:
            break

    return training_files


def doEM(in_fast5, reference, destination, strand="Template", bwa_index=None, in_hmm=None, iterations=50):
    """
    :param in_fast5: path to the fast5 file
    :param reference: path to the reference sequence (fasta)
    :param strand, string, "Template" or "Complement"
    :param bwa_index: path to index files (no suffix)
    :param in_hmm: input hmm to continue training
    :param destination: place to put the trained hmm
    :param iterations: number of iterations to do
    :return: void, makes a trained hmm (or updates the existing one)
    """
    # Preamble: first we align the 2D read to the reference and get the orientation, ie does
    # it align to the reference as it is in the fastA or to the reverse complement

    # containers and defaults
    temp_dir = destination + "tempFiles/"
    # make the temp directory, if needed
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))
    temp_np_read = temp_dir + "temp_nanoporeRead.npRead"
    temp_2d_read = temp_dir + "temp_2d_read.fa"
    make_temp_npRead(in_fast5, temp_np_read, temp_2d_read)

    # if there is no bwa index given, make one
    if bwa_index is None:
        bwa = Bwa(reference)
        bwa.build_index(temp_dir)
        bwa_ref_index = temp_dir + "temp_bwaIndex"
    else:
        bwa_ref_index = bwa_index

    # align with bwa
    bwa_dir = "/Users/Rand/projects/BGCs/submodules/bwa/"  # todo require bwa in path remove this
    command = "{bwaDir}bwa mem -x ont2d {index} {query}".format(bwaDir=bwa_dir, index=bwa_ref_index,
                                                                query=temp_2d_read)
    # this is a small SAM file that comes from bwa
    aln = subprocess.check_output(command.split())
    aln = aln.split("\t") # split

    # forward strand
    if aln[7] == "0":
        forward = True

    # backward strand
    if aln[7] == "16":
        forward = False

    # didn't map
    elif (aln[7] != "0") and (aln[7] != "16"):
        print("trainModels - read didn't map", file=sys.stderr)
        return

    # EM training routine: now we can run the training, we run training on either the template or
    # complement, so that we can run this program in parallel and really do both at once

    # containers and defaults
    temp_ref_seq = temp_dir + "temp_ref_seq.txt"
    path_to_vanillaAlign = "./vanillaAlign"  # todo could also require sonlib/bin in path
    strand_flags = ["-y", "-t"]
    training_hmm = destination + "{strand}_trained.hmm".format(strand=strand)

    # make the temp sequence file
    make_temp_sequence(reference, forward, temp_ref_seq)

    # switch flags if we're doing the complement
    if strand == "Complement":
        strand_flags = ["-z", "-c"]

    print("\ntrainModels - starting B-W on file: {inFile}".format(inFile=in_fast5), end="\n", file=sys.stderr)
    # training commands
    em_command_start = "{vanillaAlign} -r {ref} -q {npRead} {outHmmFlag} {outHmm} -i {iter}".format(
        vanillaAlign=path_to_vanillaAlign, ref=temp_ref_seq, npRead=temp_np_read,
        outHmmFlag=strand_flags[1], outHmm=training_hmm, iter=iterations)
    em_command_continue = "{vanillaAlign} -r {ref} -q {npRead} {inHmmFlag} {inHmm} {outHmmFlag} {outHmm} -i {iter}"\
        .format(vanillaAlign=path_to_vanillaAlign, inHmmFlag=strand_flags[0], inHmm=in_hmm,
                ref=temp_ref_seq, npRead=temp_np_read, outHmmFlag=strand_flags[1],
                outHmm=training_hmm, iter=iterations)

    # if we are given an input HMM
    if in_hmm is None:
        os.system(em_command_start)

    # if we're starting from scratch
    else:
        os.system(em_command_continue)


def main(args):
    # parse command line arguments
    args = parse_args()

    start_message = """\n
    Starting Baum-Welch training.
    Directory with training files: {files_dir}
    Training on {amount} bases.
    Using reference sequence: {ref}
    Was given hmm {inHmm} to start with.
    Writing trained hmm to: {outLoc}
    Training on strand: {strand}
    Performing {iterations} iterations.
    \n
    """.format(files_dir=args.files_dir, amount=args.amount, ref=args.ref, inHmm=args.inHmm,
               outLoc=args.out, strand=args.strand, iterations=args.iter)
    print(start_message, file=sys.stderr)

    # make directory to put temporary files
    temp_dir = args.out + "tempFiles_{strand}/".format(srand=args.strand)
    # make the temp directory, if needed
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))

    # index the reference for bwa, if needed
    bwa_ref_index = ''
    if args.bwa_index is None:
        bwa = Bwa(args.ref)
        bwa.build_index(temp_dir)
        bwa_ref_index = temp_dir + "temp_bwaIndex"
    else:
        bwa_ref_index = args.bwa_index

    # get a random list of files containing the number of bases we want to train
    training_file_list = cull_training_files(args.files_dir, args.amount)

    training_hmm = args.out + "{strand}_trained.hmm".format(strand=args.strand)

    # get started. if there is no input training hmm, then we start from scratch, and train on one file
    if args.inHmm is None:
        get_started_file = training_file_list.pop()
        doEM(in_fast5=get_started_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=None, destination=args.out, iterations=args.iter)

    # otherwise train with the input hmm on one file
    if args.inHmm is not None:
        get_started_file = training_file_list.pop()
        doEM(in_fast5=get_started_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=args.inHmm, destination=args.out, iterations=args.iter)

    # train on the rest of the files
    for training_file in training_file_list:
        doEM(in_fast5=training_file, reference=args.ref, strand=args.strand, bwa_index=bwa_ref_index,
             in_hmm=training_hmm, destination=args.out, iterations=args.iter)

    print("\nFinished Training routine, exiting.\n", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

