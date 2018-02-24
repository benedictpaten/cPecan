from __future__ import print_function

from argparse import ArgumentParser
import os
import sys
import pysam
import subprocess
from datetime import datetime
import glob
import numpy as np
import string

# nucleotides
NUC_A = "A"
NUC_C = "C"
NUC_G = "G"
NUC_T = "T"
NUC_GAP = "_"
NUCLEOTIDES = [NUC_A, NUC_C, NUC_G, NUC_T, NUC_GAP]

# misc
POS = "#"

# metadata for output
META_READ_ID = "read_id"
META_ALIGN_FILE = "alignment_file"
META_REFERNCE_FILE = "reference_file"
META_CONTIG = "contig"
META_PECAN_CMD = "pecan_cmd"
META_FORWARD = 'forward'
META_DATE_GENERATED = "date_generated"
META_REF_START = 'reference_start_pos'
METADATA = [META_READ_ID, META_ALIGN_FILE, META_REFERNCE_FILE, META_CONTIG, META_PECAN_CMD, META_FORWARD,
            META_DATE_GENERATED, META_DATE_GENERATED]

# validation
VALID_IDENTITY_RATIO = .85
V_LENGTH = "length"
V_CONSENSUS_IDENTITY = "consensus_identity"
V_POSTERIOR_IDENTITY = "posterior_identity"
V_NO_GAP_POST_IDENT = "no_gap_posterior_identity"
V_METADATA = "metadata"

# tmp signalAlign porting
READ_NAME_KEY = "read_name"
FAST5_NAME_KEY = "fast5_name"
POSITIVE_STRAND_KEY = "positive_strand"
ALIGNED_IDENTITY_KEY="aligned_identity"
NO_GAP_ALIGNED_IDENTITY_KEY="no_gap_aligned_identity"


def parse_args(args=None):
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--ref', '-r', action='store', dest='ref', required=True,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument("--alignment_file", '-a', action='store', dest="aln", default="*.?am",
                        help="glob to aligned reads (.sam or .bam)")
    parser.add_argument('--output_directory', '-o', action='store', dest='out', default=".",
                        help="directory to put the alignments")
    parser.add_argument('--workdir_directory', '-w', action='store', dest='workdir', default="tmp",
                        help="workdir for temporary files")
    parser.add_argument('--validate', '-v', action='store_true', dest='validate', default=False,
                        help="validate files in output directory after normal execution")
    parser.add_argument('--validate_directory', action='store', dest='validate_dir', default=None,
                        help="validate files in directory instead of normal execution")
    parser.add_argument('--lastz_exe', action='store', dest='lastz', default='cPecanLastz',
                        help="location of the cPecanLastz executable")
    parser.add_argument('--realign_exe', action='store', dest='realign', default='cPecanRealign',
                        help="location of the cPecanRealign executable")
    parser.add_argument('--keep_temp_files', action='store_true', dest='keep_temp', default=False,
                        help="Keep temporary files")
    parser.add_argument('--failed_align_directory', action='store', dest='fail', default=None,
                        help="Special location to store failed alignment files")

    return parser.parse_args()


def log(msg):
    print(msg, file=sys.stderr)


def get_reference_map(ref_location):
    reference_map = dict()
    with open(ref_location, 'r') as ref_in:
        contig = None
        fasta = None
        for line in ref_in:
            if len(line.strip()) == 0: continue
            if line.startswith("#"): continue #shouldn't happen, but just in casies
            # new contig
            if line.startswith(">"):
                # save old contig
                if contig is not None:
                    reference_map[contig] = "".join(fasta)
                # read new contig
                contig = line.strip().lstrip(">")
                fasta = list()
            else:
                # sanity check
                assert fasta is not None, "Malformed reference: {}".format(ref_location)
                fasta.append(line.strip())
    # save last fasta
    reference_map[contig] = "".join(fasta)
    return reference_map


def write_probabilities_out(out_loc, probabilities, contig, ref_start, metadata, args):
    with open(out_loc, 'w') as output:
        # header
        metadata_keys = metadata.keys()
        metadata_keys.sort()
        for key in metadata_keys:
            output.write("##{}:{}\n".format(key, metadata[key]))
        output.write("#CHROM\tPOS\tpA\tpC\tpG\tpT\tp_\n")

        # probabilities
        for prob in probabilities:
            line = [contig, (ref_start + prob[POS]), prob[NUC_A], prob[NUC_C], prob[NUC_G], prob[NUC_T], prob[NUC_GAP]]
            output.write("\t".join(map(str,line)) + "\n")


def calculate_nucleotide_probs(aln_loc, read_str, args):
    # prep
    pos_alignments = {}

    # read alignment
    with open(aln_loc, 'r') as aln_in:
        for line in aln_in:
            if line.startswith("#") or len(line.strip()) == 0: continue
            line = line.split()
            ref_pos = int(line[0])
            read_pos = int(line[1])
            prob = float(line[2])
            if ref_pos not in pos_alignments:
                pos_alignments[ref_pos] = {n:0.0 for n in NUCLEOTIDES}
            read_char = read_str[read_pos].upper()
            pos_alignments[ref_pos][read_char] += prob

    # update missing probs (which should be treated as 'gap')
    max_pos = 0
    for ref_pos in pos_alignments.keys():
        probs = pos_alignments[ref_pos]
        total = sum(probs.values())
        if total > 1.0:
            if total >= 1.1: log("\tAlignment file {} has prob {} at reference pos {}".format(aln_loc, total, ref_pos))
            probs = {n:(probs[n]/total) for n in NUCLEOTIDES}
        probs[NUC_GAP] += max(0.0, 1.0 - sum(probs.values()))
        if ref_pos > max_pos: max_pos = ref_pos

    # get ordered probs
    nucleotide_probs = list()
    for pos in xrange(0, max_pos + 1):
        if pos not in pos_alignments:
            probs = {n:(0.0 if n != NUC_GAP else 1.0) for n in NUCLEOTIDES}
        else:
            probs = pos_alignments[pos]
        # save probs
        probs[POS] = pos
        nucleotide_probs.append(probs)

    # return it
    return nucleotide_probs


def run_pecan(read, reference_map, alignment_file, args):
    #prep
    workdir = args.workdir
    #todo assert soft clipping done appropriately

    # reference
    contig = read.reference_name
    assert contig in reference_map, "Contig {} not found in {} (from {})".format(
        contig, [reference_map.keys()], args.ref)
    full_reference = reference_map[contig]
    ref_start = read.reference_start
    ref_end = read.reference_end
    ref_str = full_reference[ref_start:ref_end + 1].upper()

    # read
    read_id = read.query_name
    read_str = read.query_alignment_sequence.upper()
    # if read.is_reverse:
    #     read_str = reverse_complement(read_str, reverse=True, complement=True)

    # write files
    ref_filename = "{}_ref.fa".format(read_id)
    read_filename = "{}_read.fa".format(read_id)
    ref_loc = os.path.join(workdir, ref_filename)
    read_loc = os.path.join(workdir, read_filename)
    cigar_loc = os.path.join(workdir, "{}_cig.txt".format(read_id))
    aln_loc = os.path.join(workdir, '{}_out.txt'.format(read_id))
    with open(ref_loc, 'w') as ref_out:
        ref_out.write(">{}_{}_REF\n{}\n".format(contig, read_id, ref_str))
    with open(read_loc, 'w') as read_out:
        read_out.write(">{}\n{}\n".format(read_id, read_str))

    # generate command
    lastz_args = "--format=cigar --ambiguous=iupac --chain"
    realign_args = "--rescoreByPosteriorProb --splitMatrixBiggerThanThis 100 --diagonalExpansion 20"
    lastz_command = "%s %s %s[nameparse=darkspace] %s[nameparse=darkspace]"\
                   % (args.lastz, lastz_args, ref_loc, read_loc)
    tee_command = "tee %s 2>/dev/null" % (cigar_loc) # gets rid of error when multiple alignments are found
    realign_command = "%s %s --outputAllPosteriorProbs %s %s %s"\
                     % (args.realign, realign_args, aln_loc, ref_loc, read_loc)
    # final command
    # cmd = "%s | %s | %s" % (lastz_command, tee_command, realign_command)
    cmd = "%s | %s | head -n1 | %s" % (lastz_command, tee_command, realign_command)

    # run realign command
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=-1)
    output, _ = process.communicate()
    sts = process.wait()
    if sts != 0: raise RuntimeError( "Command exited with non-zero status %i: %s" % (sts, cmd))

    # get alignment
    success = True
    if not os.path.isfile(aln_loc):
        log("\tAlignment output not created: {}".format(aln_loc))
        log("\t\tread_len:%6d\tref_len:%6d\ttags:%s" % (len(read_str), len(ref_str), read.get_tags()))
        if args.fail is not None:
            os.rename(ref_loc, os.path.join(args.fail, ref_filename))
            os.rename(read_loc, os.path.join(args.fail, read_filename))
        success = False
    else:
        nucleotide_probs = calculate_nucleotide_probs(aln_loc, read_str, args)

        # write it out
        out_loc = os.path.join(args.out, "{}.tsv".format(read_id))
        metadata = {META_READ_ID: read_id, META_ALIGN_FILE: alignment_file, META_REFERNCE_FILE: args.ref,
                    META_CONTIG: contig, META_PECAN_CMD: cmd, META_REF_START: ref_start,
                    META_DATE_GENERATED: str(datetime.now()),
                    META_FORWARD: None if read.is_reverse is None else not read.is_reverse}
        write_probabilities_out(out_loc, nucleotide_probs, contig, ref_start, metadata, args)
        assert os.path.isfile(out_loc), "Output file {} not created".format(out_loc)

    # cleanup
    if not args.keep_temp:
        if os.path.exists(ref_loc): os.remove(ref_loc)
        if os.path.exists(read_loc): os.remove(read_loc)
        if os.path.exists(cigar_loc): os.remove(cigar_loc)
        if os.path.exists(aln_loc): os.remove(aln_loc)

    return success


def realign_alignments(alignment_filename, reference_map, args):
    # prep
    samfile = None
    read_count = 0
    failure_count = 0

    # read alignments
    try:
        log("Alignment file {}:".format(alignment_filename))
        samfile = pysam.AlignmentFile(alignment_filename, 'rb' if alignment_filename.endswith("bam") else 'r')
        for read in samfile.fetch():
            # init
            read_count += 1

            # invoke pecan
            success = run_pecan(read, reference_map, alignment_filename, args)

            if not success:
                failure_count += 1

    # close
    finally:
        if samfile is not None: samfile.close()

    log("\tHandled {} reads".format(read_count))
    if failure_count > 0: log("\t%d failures (%2.3f%%)" % (failure_count, 100.0 * failure_count / read_count))

    return read_count



def reverse_complement(dna, reverse=True, complement=True):
    # Make translation table
    trans_table = string.maketrans('ATGCatgc', 'TACGtacg')

    # Make complement to DNA
    comp_dna = dna.translate(trans_table)

    # Output all as strings
    if reverse and complement:
        return comp_dna[::-1]
    if reverse and not complement:
        return dna[::-1]
    if complement and not reverse:
        return comp_dna
    if not complement and not reverse:
        return dna


def validate_snp_directory(snp_directory, reference_map, alignment_sam_location=None, print_indiv_summary=False):
    # prep
    all_consensus_identities = list()
    all_posterior_identities = list()
    all_consensus_iden_ratio = list()
    all_posterior_iden_ratio = list()
    no_gap_post_iden_ratio = list()
    all_lengths = list()
    all_summaries = list()
    files_glob = os.path.join(snp_directory, "*.tsv")
    files = glob.glob(files_glob)
    if len(files) == 0:
        log("\tNo files matching '{}'\n".format(files_glob))
        return

    log("\tValidating {} files in '{}'".format(len(files), snp_directory))

    for file in files:
        summary, problem = validate_snp_file(file, reference_map, print_summary=print_indiv_summary)
        consensus_identity = summary[V_CONSENSUS_IDENTITY]
        posterior_identity = summary[V_POSTERIOR_IDENTITY]
        all_consensus_identities.append(consensus_identity)
        all_posterior_identities.append(posterior_identity)
        length = summary[V_LENGTH]
        consensus_ratio = 1.0 * consensus_identity / length
        posterior_ratio = 1.0 * posterior_identity / length
        all_consensus_iden_ratio.append(consensus_ratio)
        all_posterior_iden_ratio.append(posterior_ratio)
        no_gap_post_iden_ratio.append(summary[V_NO_GAP_POST_IDENT])
        all_lengths.append(length)
        all_summaries.append(summary)

    # printing results
    log("\nSummary of {} files:".format(len(files)))
    log("\tAVG Identity Ratio:       {}".format(np.mean(all_consensus_identities)))
    log("\tAVG Length:               {}".format(np.mean(all_lengths)))
    log("\tOverall P Identity Ratio: {}".format(1.0 * sum(all_posterior_identities) / sum(all_lengths)))
    log("\tOverall C Identity Ratio: {}".format(1.0 * sum(all_consensus_identities) / sum(all_lengths)))
    log("\tPosterior Identity Avg:   {}".format(np.mean(all_posterior_iden_ratio)))
    log("\tConsensus Identity Avg:   {}".format(np.mean(all_consensus_iden_ratio)))
    log("\tNo Gap Post Identity Avg: {}".format(np.mean(no_gap_post_iden_ratio)))


def validate_snp_file(snp_file, reference_map, print_summary=False, print_sequences=False):
    identifier = os.path.basename(snp_file)
    consensus_sequence = list()
    all_probabilities = list()
    header_positions = list()
    header_characters = list()
    first_pos = None
    last_pos = None
    problem = False
    duplicated_positions = 0
    unspecified_positions = 0
    full_reference_sequence = None

    metadata = {}
    with open(snp_file, 'r') as snp:
        for line in snp:
            if line.startswith("##"):
                for key in METADATA:
                    if key in line: metadata[key] = line.split(":")[1].strip()
                continue
            elif line.startswith("#"):
                # identifier = "{}/{}".format(read_name, fast5_name)
                line = line.split("\t")
                i = 0
                for l in line:
                    if l.startswith("p"):
                        header_positions.append(i)
                        header_characters.append(l[1].upper())
                    i += 1

                # validation and prep
                if META_CONTIG not in metadata: raise Exception("Missing {} header for {}" .format(META_CONTIG, identifier))
                if metadata[META_CONTIG] not in reference_map: raise Exception("Contig {} not found in reference" .format(metadata[META_CONTIG]))
                if META_READ_ID in metadata: identifier = metadata[META_READ_ID]
                if META_FORWARD in metadata: metadata[META_FORWARD] = \
                    True if metadata[META_FORWARD] == 'True' else (False if metadata[META_FORWARD] == 'False' else None)
                full_reference_sequence = reference_map[metadata[META_CONTIG]]
            else:
                line = line.split("\t")
                # positions
                pos = int(line[1])
                # set first_position (for reference matching)
                if first_pos is None: first_pos = pos
                # for cases where positions are duplicated or missing
                if last_pos is not None:
                    if last_pos >= pos:
                        duplicated_positions += 1
                        continue
                    while last_pos + 1 < pos:
                        probabilities = dict()
                        for char in header_characters:
                            probabilities[char] = 0.0
                        all_probabilities.append(probabilities)
                        consensus_sequence.append("-")
                        unspecified_positions += 1
                        last_pos += 1
                last_pos = pos
                #probabilities
                probabilities = dict()
                #consensus
                max_prob = -1.0
                max_prob_idx = None
                idx = 0
                for pos in header_positions:
                    prob = float(line[pos].strip())
                    probabilities[header_characters[idx]] = prob
                    if prob > max_prob:
                        max_prob = prob
                        max_prob_idx = idx
                    idx += 1
                consensus_sequence.append(header_characters[max_prob_idx])
                all_probabilities.append(probabilities)

    # get sequences
    consensus_sequence = "".join(consensus_sequence)
    reference_sequence = full_reference_sequence[min(first_pos, last_pos):max(first_pos, last_pos)+1].upper()

    # this is our quality metric
    no_gap_length = 0
    no_gap_post_iden = 0.0
    length = 0
    identity = 0
    posterior_identity = 0.0
    for c, r, p in zip(consensus_sequence, reference_sequence, all_probabilities):
        # calculate with gaps
        length += 1
        if c == r: identity += 1
        if r != 'N': posterior_identity += p[r]
        # remove gap calc
        if c != '_':
            no_gap_length += 1
            no_gap_post_iden += p[r]

    # for separating results
    if print_sequences: log("")

    # sanity check
    if duplicated_positions > 0:
        log("{}: Found {} duplicated positions!"
              .format(identifier, duplicated_positions))
    if unspecified_positions * 100 > length:
        log("{}: Found {} unspecified positions ({}% of total length)"
              .format(identifier, unspecified_positions, int(100.0 * unspecified_positions / length)))
        problem = True

    if print_summary:
        strand_char = " " if META_FORWARD not in metadata else ("+" if metadata[META_FORWARD] else "-")
        log("%s:\tstrand: %s\tlen:%6d\tconsensus_iden:%6d (%.3f)\tposterior_iden:%6d (%.3f)\t"
            "gap_ratio: %.3f\tno_gap_length:%6d\tno_gap_post_iden:%6d (%.3f)" %
              (identifier, strand_char, length, identity, 1.0*identity/length,
               int(posterior_identity), posterior_identity / length, 1.0*(length-no_gap_length)/length,
               no_gap_length, int(no_gap_post_iden), 0.0 if no_gap_length == 0 else (no_gap_post_iden/no_gap_length)))

        # printing full sequences
        if print_sequences:
            log("\treference:  {}".format(reference_sequence))
            log("\tconsensus:  {}".format(consensus_sequence))

    summary = {
        V_LENGTH: length,
        V_CONSENSUS_IDENTITY: identity,
        V_POSTERIOR_IDENTITY: posterior_identity,
        V_NO_GAP_POST_IDENT: 0.0 if no_gap_length == 0 else (no_gap_post_iden/no_gap_length),
        V_METADATA: metadata,

    }
    return summary, problem


def main():
    args = parse_args()
    assert os.path.isfile(args.ref), "--ref argument does not exist: {}".format(args.ref)

    # validate action
    if args.validate_dir is not None:
        reference_map = get_reference_map(args.ref)
        validate_snp_directory(args.validate_dir, reference_map, print_indiv_summary=True)
        return

    # not validate action - continue sanity checks
    if not os.path.isdir(args.out): os.mkdir(args.out)
    assert os.path.isdir(args.out), "--output_directory argument could not be found/created: " \
                                                "{}".format(args.out)
    if not os.path.isdir(args.workdir): os.mkdir(args.workdir)
    assert os.path.isdir(args.workdir), "--workdir_directory argument could not be found/created: " \
                                                "{}".format(args.workdir)
    if args.fail is not None:
        if not os.path.isdir(args.fail): os.mkdir(args.fail)
        assert os.path.isdir(args.fail), "--failed_align_directory argument could not be found/created: " \
                                                    "{}".format(args.fail)


    alignment_files = glob.glob(args.aln)
    assert len(alignment_files) > 0, "Could not find files matching {}".format(args.aln)

    # load reference
    reference_map = get_reference_map(args.ref)

    # generate alignments
    for alignment_file in alignment_files:
        realign_alignments(alignment_file, reference_map, args)

    # maybe validate
    if args.validate:
        log("Validating output")
        validate_snp_directory(args.out, reference_map, print_indiv_summary=False)

    # cleanup
    if len(os.listdir(args.workdir)) == 0:
        log("Removing empty workdir")
        os.rmdir(args.workdir)
    log("Fin.")



if __name__ == "__main__":
    main()