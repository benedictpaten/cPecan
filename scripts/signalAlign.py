#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import sys
sys.path.append("../")
from nanoporeLib import *
from threading import Thread
from multiprocessing import Pool
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
                        required=False, type=str, help="for single read alignment")
    parser.add_argument('-nb_files', '-nb', action='store', dest='nb_files', required=False,
                        default=50, type=int, help="maximum number of reads to align")

    # multithreading
    parser.add_argument('--threads', '-t', action='store', dest='threads', required=False,
                        type=int, default=2, help="number of threads to use, default=2")

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
    parser.add_argument('--strawMan', '-sm', action='store_true', dest='strawMan',
                        required=False, default=False, help="use strawman pair-HMM")

    args = parser.parse_args()
    return args


class SignalAlignment(object):
    def __init__(self, reference, destination, strawman, bwa_index, in_templateHmm, in_complementHmm):
        self.reference = reference
        self.destination = destination
        self.strawman = strawman
        self.bwa_index = bwa_index
        self.in_templateHmm = in_templateHmm
        self.in_complementHmm = in_complementHmm
        self.in_templateModel = None
        self.in_complementModel = None

    def do_alignment(self, in_fast5):
        # Preamble set up before doing the alignment

        # containers and defaults
        read_label = in_fast5.split("/")[-1]  # used in the posteriors file
        read_name = in_fast5.split("/")[-1][:-6]  # get the name without the '.fast5'
        temp_folder = FolderHandler()
        temp_dir_path = temp_folder.open_folder(self.destination + "tempFiles_{readLabel}".format(readLabel=read_label))

        temp_np_read = temp_folder.add_file_path("temp_{read}.npRead".format(read=read_label))
        temp_2d_read = temp_folder.add_file_path("temp_2Dseq_{read}.fa".format(read=read_label))
        temp_t_model = temp_folder.add_file_path("template_model.model")
        temp_c_model = temp_folder.add_file_path("complement_model.model")


        # make the npRead and fasta todo make this assert
        success, temp_t_model, temp_c_model = get_npRead_2dseq_and_models(fast5=in_fast5,
                                                                          npRead_path=temp_np_read,
                                                                          twod_read_path=temp_2d_read,
                                                                          template_model_path=temp_t_model,
                                                                          complement_model_path=temp_c_model)

        print("temp template model", temp_t_model)
        print("temp complement model", temp_c_model)

        if success is False:
            return False

        # add an indicator for the model being used
        if self.strawman is True:
            model_label = ".sm"
            use_strawMan_flag = "--s "
        else:
            model_label = ".vl"
            use_strawMan_flag = ""

        # this gives the format: /directory/for/files/file.model.orientation.tsv
        posteriors_file_path = ''

        # get orientation from BWA
        orientation = orient_read_with_bwa(bwa_index=self.bwa_index, query=temp_2d_read)

        # forward strand
        if orientation == 0:
            forward = True
            posteriors_file_path = self.destination + read_name + model_label + ".forward.tsv"

        # backward strand
        if orientation == 16:
            forward = False
            posteriors_file_path = self.destination + read_name + model_label + ".backward.tsv"

        # didn't map
        elif (orientation != 0) and (orientation != 16):
            print("\n\ntrainModels - read didn't map", file=sys.stderr)
            return False

        # Alignment routine

        # containers and defaults
        temp_ref_seq = temp_folder.add_file_path("temp_ref_seq.txt")

        path_to_vanillaAlign = "./vanillaAlign"  # todo could require this in path

        # make sequence for vanillaAlign, we orient the sequence so that the template events align to the
        # reference and the complement events align to the reverse complement of the reference
        make_temp_sequence(self.reference, forward, temp_ref_seq)

        # alignment flags

        # input (match) models
        if self.in_templateModel is not None:
            template_model_flag = "-T {model_loc} ".format(model_loc=self.in_templateModel)
        if temp_t_model is not None:
            template_model_flag = "-T {t_model} ".format(t_model=temp_t_model)
        else:
            template_model_flag = ""
        if self.in_complementModel is not None:
            complement_model_flag = "-C {model_loc} ".format(model_loc=self.in_complementModel)
        if temp_c_model is not None:
            complement_model_flag = "-C {c_model} ".format(c_model=temp_c_model)
        else:
            complement_model_flag = ""

        # input HMMs
        if self.in_templateHmm is not None:
            template_hmm_flag = "-y {hmm_loc} ".format(hmm_loc=self.in_templateHmm)
        else:
            template_hmm_flag = ""
        if self.in_complementHmm is not None:
            complement_hmm_flag = "-z {hmm_loc} ".format(hmm_loc=self.in_complementHmm)
        else:
            complement_hmm_flag = ""

        # alignment commands
        alignment_command = \
            "{vA} {straw}-r {ref} -q {npRead} {t_model}{c_model}{t_hmm}{c_hmm} -u {posteriors} -L {readLabel}"\
            .format(vA=path_to_vanillaAlign, straw=use_strawMan_flag, ref=temp_ref_seq, readLabel=read_label,
                    npRead=temp_np_read, t_model=template_model_flag, c_model=complement_model_flag,
                    t_hmm=template_hmm_flag, c_hmm=complement_hmm_flag, posteriors=posteriors_file_path)

        # run
        print("signalAlign - running command", alignment_command, end="\n", file=sys.stderr)
        os.system(alignment_command)
        #temp_folder.remove_folder()
        return True


def main(args):
    # parse args
    args = parse_args()

    start_message = """
    Starting Signal Align
    """
    print(start_message, file=sys.stderr)

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

    # if only one file is specified for alignment

    alignment_args = {
            "reference": args.ref,
            "destination": temp_dir_path, #+ args.out,
            "strawman": args.strawMan,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
    }

    alignment = SignalAlignment(**alignment_args)

    if args.files_dir is None:
        # run the alignment
        aln_ = alignment.do_alignment(in_fast5=args.single_file)

        if aln_ is True:
            print("signalAlign - aligned read {f} successfully".format(f=args.single_file), file=sys.stderr)
        if aln_ is False:
            print("signalALign - error while aligning {f}".format(f=args.single_file), file=sys.stderr)
        #temp_folder.remove_folder()

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
            # log which file we're on
            print("On file:{}".format(i), file=sys.stderr)

            f = fast5s[i]  # get the file from the list
            f = args.files_dir + f  # add the path to it

            aln_ = alignment.do_alignment(in_fast5=f)

            if aln_ is True:
                continue  # successful alignment log is within vanillaAlign
            if aln_ is False:
                print("signalAlign - error while aligning {f}".format(f=f), file=sys.stderr)
                continue


if __name__ == "__main__":
    sys.exit(main(sys.argv))
