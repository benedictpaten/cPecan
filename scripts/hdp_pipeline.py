#!/usr/bin/env python
"""Master pipeline script for generating trained HDPs for MinION signal data

Input: alignments using non-HDP model, desired HDP type
Output: trained HDP and model

The objective of this pipeline is to:
    1. use input alignments to make 'build alignment' for generating the initial HDP
    2. generates the initial HDP
    3. trains the HDP on MinION reads
    4. outputs distributions for all kmers from the HDP
"""

import os
from argparse import ArgumentParser
from subprocess import check_call, Popen
from shutil import copyfile


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # build alignment
    parser.add_argument('--C_alignments', '-C', action='store',
                        dest='C_alns', required=False, type=str, default=None,
                        help="C files")
    parser.add_argument('--mC_alignments', '-mC', action='store',
                        dest='mC_alns', required=False, type=str, default=None,
                        help="mC files")
    parser.add_argument('--hmC_alignments', '-hmC', action='store',
                        dest='hmC_alns', required=False, type=str, default=None,
                        help="hmC files")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    # initial HDP
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.2, dest='threshold')
    parser.add_argument('--hdp_type', action='store', type=str, required=True, dest='hdp_type',
                        help="Build Hdp, specify type, options: "
                             "singleLevelFixed, singleLevelPrior, multisetFixed, multisetPrior")
    # train models
    parser.add_argument('--file_directory', '-d', action='append', default=None,
                        dest='files_dir', required=False, type=str,
                        help="directories with fast5 files to train on")
    parser.add_argument('--ref', '-r', action='store', default=None,
                        dest='ref', required=False, type=str,
                        help="location of refrerence sequence in FASTA")
    parser.add_argument('--iterations', '-i', action='store', dest='iter', default=10,
                        required=False, type=int)
    parser.add_argument('--train_amount', '-a', action='store', dest='amount', default=15,
                        required=False, type=int,
                        help="limit the total length of sequence to use in training.")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False, default=4,
                        type=int, help="number of jobs to run concurrently")
    parser.add_argument('--cytosine_substitution', '-cs', action='append', default=None,
                        dest='cytosine_sub', required=False, type=str,
                        help="mutate cytosines to this letter in the reference")

    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


# Pipeline Script
args = parse_args()  # parse arguments
working_directory = args.out  # this is the directory we will use for everything
pipeline_log = open(working_directory + "pipeline.log", 'a')

# build alignment
cPecan_dir = "../../cPecan/"
build_alignment_location = working_directory + "buildAlignment.tsv"
build_alignment_command = "{cP}scripts/makeBuildAlignments.py -o={bA} -t={threshold} " \
                          "".format(cP=cPecan_dir, C=args.C_alns, mC=args.mC_alns, threshold=args.threshold,
                                    hmC=args.hmC_alns, bA=build_alignment_location)
if args.C_alns is not None: build_alignment_command += "-C={C} ".format(C=args.C_alns)
if args.mC_alns is not None: build_alignment_command += "-mC={mC} ".format(mC=args.mC_alns)
if args.hmC_alns is not None: build_alignment_command += "-mC={hmC} ".format(hmC=args.hmC_alns)
pipeline_log.write("[pipeline] NOTICE: Making build alignment using files from:\n\t{C}\n\t{mC}\n\t{hmC}\n"
                   "".format(C=args.C_alns, mC=args.mC_alns, hmC=args.hmC_alns))
check_call(build_alignment_command.split(), stderr=pipeline_log, stdout=pipeline_log)

# initial HDP
pipeline_log.write("[pipeline] NOTICE: Making initial HDP of type {}\n".format(args.hdp_type))
template_hdp_location = working_directory + "template." + args.hdp_type + ".nhdp"
complement_hdp_location = working_directory + "complement." + args.hdp_type + ".nhdp"
initial_hdp_build_out = open(working_directory + "build_initial_hdp.out", 'w')
initial_hdp_build_err = open(working_directory + "build_initial_hdp.err", 'w')
assert (os.path.isfile(build_alignment_location)), "ERROR: Didn't find build alignment"
build_initial_hdp_command = "./trainModels --buildHdp={hdpType} -tH={tHdpLoc} -cH={cHdpLoc} -al={buildAln}" \
                            "".format(hdpType=args.hdp_type, tHdpLoc=template_hdp_location,
                                      cHdpLoc=complement_hdp_location, buildAln=build_alignment_location)
check_call(build_initial_hdp_command.split(), stdout=initial_hdp_build_out, stderr=initial_hdp_build_err)
initial_hdp_build_out.close()
initial_hdp_build_err.close()

# trainModels
pipeline_log.write("[pipeline] NOTICE: Training HDP models.\n")
template_trained_hdp_location = working_directory + "template_trained." + args.hdp_type + ".nhdp"
complement_trained_hdp_location = working_directory + "complement_trained." + args.hdp_type + ".nhdp"
# make a copy of the HDP files so we can compare before and after training
copyfile(template_hdp_location, template_trained_hdp_location)
copyfile(complement_hdp_location, complement_trained_hdp_location)
train_hdp_out = open(working_directory + "train_hdp.out", 'w')
train_hdp_err = open(working_directory + "train_hdp.err", 'w')
train_models_command = "./trainModels -r={ref} -i={iter} -a={amount} -smt=threeStateHdp -tH={tHdp} " \
                       "-cH={cHdp} -o={wd} -t={threshold} " \
                       "".format(ref=args.ref, iter=args.iter, amount=args.amount, tHdp=template_trained_hdp_location,
                                 cHdp=complement_trained_hdp_location, wd=working_directory, threshold=args.threshold)
for directory in args.files_dir:
    train_models_command += "-d={dir} ".format(dir=directory)
if args.cytosine_sub is not None:
    for substitution in args.cytosine_sub:
        train_models_command += "-cs={sub} ".format(sub=substitution)
check_call(train_models_command.split(), stdout=train_hdp_out, stderr=train_hdp_err)

# get HDP distributions
pipeline_log.write("[pipeline] NOTICE: running compareDistributions.\n")
template_distr_dir = working_directory + "template_distrs/"
complement_distr_dir = working_directory + "complement_distrs/"
os.makedirs(template_distr_dir)
os.makedirs(complement_distr_dir)
compare_distributions_commands = [
    "./compareDistributions {tHdp} {tDir}".format(tHdp=template_trained_hdp_location, tDir=template_distr_dir),
    "./compareDistributions {cHdp} {cDir}".format(cHdp=complement_trained_hdp_location, cDir=complement_distr_dir)
]
procs = [Popen(x.split(), stdout=pipeline_log, stderr=pipeline_log) for x in compare_distributions_commands]
status = [p.wait() for p in procs]

pipeline_log.write("[pipeline] DONE.\n")
pipeline_log.close()





