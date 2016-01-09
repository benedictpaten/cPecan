#!/usr/bin/env python

from __future__ import division, print_function
import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from random import sample


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--set1', '-s1', action='store',
                        dest='set1', required=True, type=str,
                        help="set1 directory")
    parser.add_argument('--set2', '-s2', action='store',
                        dest='set2', required=True, type=str,
                        help="set2 directory")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold',
                        required=False, default=0.2, type=float)
    parser.add_argument('--out', '-o', action='store',
                        dest='dump', required=False, type=str, default="./",
                        help="place to put the results")
    parser.add_argument('--set1_label', '-s1l', action='store', dest='set1_label', required=False,
                        default=None)
    parser.add_argument('--set2_label', '-s2l', action='store', dest='set2_label', required=False,
                        default=None)

    args = parser.parse_args()
    return args

def hash_reads_to_paths(path_to_files, forward):
    if forward:
        file_suffix = ".forward.tsv"
    else:
        file_suffix = ".backward.tsv"

    read_hash = {}

    for f in os.listdir(path_to_files):
            if f.endswith(file_suffix) and os.stat(path_to_files + f).st_size != 0:
                read_hash[f.split("_strand")[0]] = path_to_files + f

    return read_hash


def get_paired_alignments(set1_dir, set2_dir, forward):
    # collect the files
    tsvs_1 = hash_reads_to_paths(set1_dir, forward)
    tsvs_2 = hash_reads_to_paths(set2_dir, forward)

    read_intersection = set(tsvs_1.keys()).intersection(set(tsvs_2.keys()))

    paired_alignments = []
    add_alignment_pair = paired_alignments.append
    for read in read_intersection:
        add_alignment_pair((tsvs_1[read], tsvs_2[read]))

    return paired_alignments


def eventalign_coordinate_correction(forward, strand):
    if strand == "t":
        if forward:
            return 0
        else:
            return -891
    if strand == "c":
        if forward:
            return -891
        else:
            return 0


def parse_aligned_pairs(table, strand, signalalign_bool, forward, threshold=0.2):
    if signalalign_bool:
        pairs = dict()
    else:
        pairs = []

    for line in table:
        if line[4] == strand:
            if signalalign_bool:
                x = int(line[0])
                y = int(line[1])
                if float(line[8]) >= threshold:
                    pairs[(x, y)] = float(line[8])
            else:
                x = abs(eventalign_coordinate_correction(forward, strand) + int(line[1]))
                y = int(line[5])
                pairs.append((x, y))
        else:
            pass
    return pairs


def aggregate_aligned_pairs2(set1_dir, set2_dir, threshold=0.01, out_path="./", set1_label=None, set2_label=None):
    total_intersecting_pairs = 0
    total_pairs = 0
    set1_intersection_posteriors = pd.DataFrame()
    set2_intersection_posteriors = pd.DataFrame()
    set1_unique_posteriors = pd.DataFrame()
    set2_unique_posteriors = pd.DataFrame()

    for forward in [True, False]:
        for x, y in get_paired_alignments(set1_dir=set1_dir, set2_dir=set2_dir,
                                          forward=forward):

            set1_df = pd.read_table(x, usecols=(1, 4, 5, 12),
                                    dtype={'ref_pos': np.int32,
                                           'event_idx': np.int32,
                                           'strand': np.str,
                                           'prob': np.float64,
                                           },
                                    header=None,
                                    names=['ref_pos', 'strand', 'event_idx', 'prob']
                                    )

            set2_df = pd.read_table(y, usecols=(1, 4, 5, 12),
                                    dtype={'ref_pos': np.int32,
                                           'event_idx': np.int32,
                                           'strand': np.str,
                                           'prob': np.float64,
                                           },
                                    header=None,
                                    names=['ref_pos', 'strand', 'event_idx', 'prob']
                                    )

            set1_df = set1_df[set1_df.prob >= threshold]
            set2_df = set2_df[set2_df.prob >= threshold]

            union = pd.merge(set1_df, set2_df, how='outer', on=['ref_pos', 'event_idx', 'strand'],
                             indicator=True)

            intersect = union[union._merge == 'both']
            x_only_posteriors = union[union._merge == "left_only"].drop(['prob_y', 'ref_pos',
                                                                         'event_idx', 'strand',
                                                                         '_merge'], 1)
            y_only_posteriors = union[union._merge == "right_only"].drop(['prob_x', 'ref_pos',
                                                                         'event_idx', 'strand',
                                                                         '_merge'], 1)
            set1_unique_posteriors = pd.concat([set1_unique_posteriors, x_only_posteriors])
            set2_unique_posteriors = pd.concat([set2_unique_posteriors, y_only_posteriors])
            x_inter_probs = intersect.prob_x
            y_inter_probs = intersect.prob_y
            set1_intersection_posteriors = pd.concat([set1_intersection_posteriors, x_inter_probs])
            set2_intersection_posteriors = pd.concat([set2_intersection_posteriors, y_inter_probs])
            total_pairs += union.shape[0]
            total_intersecting_pairs += intersect.shape[0]

    jaccard_idx = total_intersecting_pairs / total_pairs

    report = """
    Comparing alignments.
    SignalAlign alignments: {set1_dir}
    EventAlign alignments: {set2_dir}
    Postior prob threshold: {threshold}
    Jaccard Index: {jaccard}
    Intersection: {inter} pairs
    SignalAlign has {set1_uniq} unique pairs
    EventAlign has {set2_uniq} unique pairs
    Posteriors written to: {posteriors}
    """.format(set1_dir=set1_dir, set2_dir=set2_dir, jaccard=jaccard_idx,
               inter=total_intersecting_pairs, set1_uniq=set1_unique_posteriors.shape[0],
               set2_uniq=set2_unique_posteriors.shape[0],
               posteriors=out_path, threshold=threshold)
    print(report)


    if set1_label is None:
        set1_label = "set1"
    if set2_label is None:
        set2_label = "set2"

    set1_unique_posteriors.to_pickle(out_path + "{}_unique_posteriors.pkl".format(set1_label))
    set1_intersection_posteriors.to_pickle(out_path + "{}_intersection_posteriors.pkl".format(set1_label))

    set2_unique_posteriors.to_pickle(out_path + "{}_unique_posteriors.pkl".format(set2_label))
    set2_intersection_posteriors.to_pickle(out_path + "{}_intersection_posteriors.pkl".format(set2_label))


def main(args):
    args = parse_args()
    aggregate_aligned_pairs2(set1_dir=args.set1, set2_dir=args.set2,
                             set1_label=args.set1_label, set2_label=args.set2_label,
                             out_path=args.dump)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

