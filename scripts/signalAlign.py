#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import sys
sys.path.append("../")
from nanoporeLib import *
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser (description=__doc__)

    # query files
    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with fast5 files to train on, should have trailing /")
    parser.add_argument('--in_file', '-f', action='store', dest='single_file',
                        required=False, type=str)
    parser.add_argument('-nb_files', '-nb', action='store', dest='nb_files', required=False,
                        default=50, type=int)

    # reference
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str)
    parser.add_argument('--index', '-i', action='store', dest='bwa_index',
                        required=False, type=str, default=None)
    # input HMMs
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None)
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None)

    # TODO make functionality to force use of new match model

    # output
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the trained model, should have trailing slash")
    parser.add_argument('--strawMan', '-sm', action='store_true', dest='strawMan',
                        required=False, default=False)

    args = parser.parse_args()
    return args


def do_alignment(in_fast5, reference, destination, strawMan_flag, posteriors_file,
                 bwa_index, in_Template_hmm=None, in_Complement_hmm=None,
                 in_Template_model=None, in_Complement_model=None):
    # Preamble set up before doing the alignment

    # containers and defaults
    read_label = in_fast5.split("/")[-1]
    temp_dir = destination + "tempFiles_{readLabel}/".format(readLabel=read_label)  # where the temp files go
    # make directory if it doesn't exist
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))
    temp_np_read = temp_dir + "temp_nanoporeRead.npRead"  # where the npRead goes
    temp_2d_read = temp_dir + "temp_2d_read.fa"  # where the fasta for the read goes
    # make the npRead and fasta
    temp_file_success = make_npRead_and_2d_seq(in_fast5,
                                               npRead_dest=temp_np_read,
                                               twod_read_dest=temp_2d_read)
    if temp_file_success is False:
        return False

    # check if we have a bwa index, and make one if needed
    #if bwa_index is None:
    #    bwa = Bwa(reference)
    #    bwa.build_index(temp_dir)
    #    bwa_ref_index = temp_dir + "temp_bwaIndex"
    #else:
    #    bwa_ref_index = bwa_index
    bwa_ref_index = bwa_index

    # get orientation from BWA
    orientation = orient_read_with_bwa(bwa_index=bwa_ref_index, query=temp_2d_read)

    # posteriors file handling
    posteriors_dest = ''
    model_label = ''

    # add an indicator for the model being used
    if strawMan_flag is True:
        model_label = ".sm"
    else:
        model_label = ".vl"

    # forward strand
    # this gives the format: /directory/for/files/file.model.orientation.tsv
    if orientation == 0:
        forward = True
        posteriors_dest = destination + posteriors_file + model_label + ".forward.tsv"

    # backward strand
    if orientation == 16:
        forward = False
        posteriors_dest = destination + posteriors_file + model_label + ".backward.tsv"

    # didn't map
    elif (orientation != 0) and (orientation != 16):
        print("\n\ntrainModels - read didn't map", file=sys.stderr)
        return False

    # Alignment routine

    # containers and defaults
    temp_ref_seq = temp_dir + "temp_ref_seq.txt" # temp file for vanillaAlign
    path_to_vanillaAlign = "./vanillaAlign"

    # make sequence for vanillaAlign, we orient the sequence so that the template events align to the
    # reference and the complement events align to the reverse complement of the reference
    make_temp_sequence(reference, forward, temp_ref_seq)

    # alignment flags

    # straw man model
    if strawMan_flag is True:
        use_strawMan_flag = "--s "
    else:
        use_strawMan_flag = ""

    # input (match) models
    if in_Template_model is not None:
        template_model_flag = "-T {hmm_loc} ".format(hmm_loc=in_Template_hmm)
    else:
        template_model_flag = ""
    if in_Complement_model is not None:
        complement_model_flag = "-C {hmm_loc} ".format(hmm_loc=in_Complement_hmm)
    else:
        complement_model_flag = ""

    # input HMMs
    if in_Template_hmm is not None:
        template_hmm_flag = "-y {hmm_loc} ".format(hmm_loc=in_Template_hmm)
    else:
        template_hmm_flag = ""
    if in_Complement_hmm is not None:
        complement_hmm_flag = "-z {hmm_loc} ".format(hmm_loc=in_Complement_hmm)
    else:
        complement_hmm_flag = ""

    # alignment commands
    alignment_command = \
        "{vA} {straw}-r {ref} -q {npRead} {t_model_flag}{c_model_flag}{t_hmm}{c_hmm} -u {posteriors} -L {readLabel}"\
        .format(vA=path_to_vanillaAlign, straw=use_strawMan_flag, ref=temp_ref_seq, readLabel=read_label,
                npRead=temp_np_read, t_model_flag=template_model_flag, c_model_flag=complement_model_flag,
                t_hmm=template_hmm_flag, c_hmm=complement_hmm_flag, posteriors=posteriors_dest)

    # run
    print("signalAlign - running command", alignment_command, end="\n", file=sys.stderr)
    os.system(alignment_command)
    return True


def main(args):
    # parse args
    args = parse_args()

    start_message = """
    Starting Signal Align
    """
    print(start_message, file=sys.stderr)

    # make directory to put temporary files
    temp_dir = args.out + "tempFiles_alignment/"
    # make the temp directory, if needed
    if not os.path.isdir(temp_dir):
        os.system("mkdir {dir}".format(dir=temp_dir))

    # index the reference for bwa, if needed
    bwa_ref_index = ''
    if args.bwa_index is None:
        print("signalAlign - indexing reference", file=sys.stderr)
        bwa_ref_index = get_bwa_index(args.ref, temp_dir)
        print("signalAlign - indexing reference, done", file=sys.stderr)
    else:
        bwa_ref_index = args.bwa_index

    # if only one file is specified for alignment
    if args.files_dir is None:
        # file handling
        in_file_name = args.single_file.split("/")[-1]

        # run the alignment
        aln_ = do_alignment(in_fast5=args.single_file, reference=args.ref, destination=args.out,
                            strawMan_flag=args.strawMan, posteriors_file=in_file_name, bwa_index=bwa_ref_index,
                            in_Template_hmm=args.in_T_Hmm, in_Complement_hmm=args.in_C_Hmm)
        if aln_ is True:
            print("signalAlign - aligned read {f} successfully".format(f=args.single_file), file=sys.stderr)
        if aln_ is False:
            print("signalALign - error while aligning {f}".format(f=args.single_file), file=sys.stderr)
    else:
        # get all the fast5s in the directory
        fast5s = [x for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

        # if we're picking a selection of them, pick randomly
        nb_files = args.nb_files
        if nb_files < len(fast5s):
            shuffle(fast5s)

        # if there are fewer files in the directory than our sampling amount, get all of them
        if nb_files > len(fast5s):
            nb_files = len(fast5s)

        # go through all the files and align them
        for i in xrange(nb_files):
            print("On file:{}\n".format(i), file=sys.stderr)
            f = fast5s[i]
            f_name = f.split("/")[-1][:-6]
            f = args.files_dir + f
            # run the alignment
            aln_ = do_alignment(in_fast5=f, reference=args.ref, destination=args.out, strawMan_flag=args.strawMan,
                                posteriors_file=f_name, bwa_index=bwa_ref_index, in_Template_hmm=args.in_T_Hmm,
                                in_Complement_hmm=args.in_C_Hmm)
            if aln_ is True:
                continue
                #print("signalAlign - aligned read {f} successfully".format(f=f), file=sys.stderr)
            if aln_ is False:
                print("signalALign - error while aligning {f}".format(f=f), file=sys.stderr)
                continue


if __name__ == "__main__":
    sys.exit(main(sys.argv))
