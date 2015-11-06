#!/usr/bin/env python
"""Turn alignedPairs from vanillaAlign into EventAlign type files
"""
from __future__ import division, print_function
import sys
import os
import numpy as np
from argparse import ArgumentParser
from nanoporeLib import *

def parse_args():
    parser = ArgumentParser (description=__doc__)

    # query files
    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with tsv files from signalAlign")
    parser.add_argument('--fast5s', '-f', action='store', dest='fast5s',
                        required=False, type=str,
                        help="directory with fast5 files")
    # reference
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=False, type=str)

    args = parser.parse_args()
    return args

def calculate_lambda(noise_mean, noise_stdev):
    return (power(noise_mean, 3)) / (power(noise_stdev, 2))


def parse_signalAlign_output(tsv_file):
    for line in open(tsv_file, 'r'):
        line = line.strip()
        line = line.split("\t")
        strand = line[0]
        x = int(line[1])
        y = int(line[2])
        p = float(line[3])
        yield (x, y, p, strand)


def main(args):
    # parse args
    args = parse_args()

    file_header = "#F/B\tRefPos\trefKmer\tRead\tstrand\teMean\teNoise\teDuration\tposterior\tE(Mean)\tE(noise)\n"
    print(file_header, file=sys.stdout)

    # get the files
    files = [x for x in os.listdir(args.files_dir) if x.endswith(".tsv")]

    # get the reference sequence
    ref_seq = ''
    for header, comment, sequence in read_fasta(args.ref):
        ref_seq = sequence
        break

    # todo check

    for f in files:
        tsv = args.files_dir + f
        f = f.split(".")
        name = f[0] + "." + f[1] + ".fast5"  # there is a better way to deal with this
        model = f[2]
        orientation = f[3]

        # if the file isn't in the directory, start the loop over
        if not os.path.isfile(args.fast5s+name):
            continue

        # load the nanopore read
        npRead = NanoporeRead(args.fast5s+name)

        # make the models
        template_model = TemplateModel(args.fast5s+name).get_model_dict()
        complement_model = ComplementModel(args.fast5s+name).get_model_dict()

        for x, y, p, strand in parse_signalAlign_output(tsv):
            # get the reference kmer
            k_x = ref_seq[x:x+6]

            # look up the observations for the event
            if strand == "t":
                e_u = npRead.template_event_table[y][0]     # event mean
                e_o = npRead.template_event_table[y][2]     # event noise
                e_t = npRead.template_event_table[y][3]     # event duration

                # adjust the model
                Ek_x = template_model[k_x]                  # unscaled model

                # level mean
                level_E_u = Ek_x[0] * npRead.template_scale + npRead.template_shift
                # level stdev
                #level_E_o = Ek_x[1] * npRead.template_var

                # noise mean
                noise_E_u = Ek_x[2] * npRead.template_scale_sd
                # noise lambda
                #noise_lam = calculate_lambda(Ek_x[2], Ek_x[3]) * npRead.template_var_sd
                # noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
                #noise_E_o = np.sqrt(power(noise_E_u, 3) / noise_lam)

            if strand == "c":
                e_u = npRead.complement_event_table[y][0]   # event mean
                e_o = npRead.complement_event_table[y][2]   # event noise
                e_t = npRead.complement_event_table[y][3]   # event duration

                # adjust the model
                Ek_x = complement_model[k_x]                # unscaled model

                # level mean
                level_E_u = Ek_x[0] * npRead.complement_scale + npRead.complement_shift
                # level stdev
                #level_E_o = Ek_x[1] * npRead.complement_var

                # noise mean
                noise_E_u = Ek_x[2] * npRead.complement_scale_sd
                # noise lambda
                #noise_lam = calculate_lambda(Ek_x[2], Ek_x[3]) * npRead.complement_var_sd
                # noise_sd = sqrt(adjusted_noise_mean**3 / adjusted_noise_lambda);
                #noise_E_o = np.sqrt(power(noise_E_u, 3) / noise_lam)

            # print the output
            print(orientation, x, k_x, name, strand, e_u, e_o, e_t, p, level_E_u, noise_E_u,
                  sep="\t", end="\n", file=sys.stdout)

            # close files
        npRead.close()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
