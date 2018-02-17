from __future__ import print_function

from argparse import ArgumentParser
import os
import sys
import pysam
import subprocess
from datetime import datetime
import glob
import shutil

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
META_DATE_GENERATED = "date_generated"



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
    parser.add_argument('--lastz_exe', action='store', dest='lastz', default='cPecanLastz',
                        help="location of the cPecanLastz executable")
    parser.add_argument('--realign_exe', action='store', dest='realign', default='cPecanRealign',
                        help="location of the cPecanRealign executable")
    parser.add_argument('--keep_temp_files', action='store_true', dest='keep_temp', default=False,
                        help="Keep temporary files")

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
        if total >= 1.0:
            if total >= 1.01: log("Alignment file {} has prob {} at reference pos {}".format(aln_loc, total, ref_pos))
            probs = {n:(probs[n]/total) for n in NUCLEOTIDES}
        probs[NUC_GAP] += max(0.0, 1.0 - total)
        if ref_pos > max_pos: max_pos = ref_pos

    # get ordered probs
    nucleotide_probs = list()
    for pos in xrange(0, max_pos + 1):
        if pos not in pos_alignments:
            probs = {n:(0.0 if n != NUC_GAP else 1.0) for n in NUCLEOTIDES}
        else:
            probs = pos_alignments[pos]
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
    ref_str = full_reference[ref_start:ref_end + 1]

    # read
    read_id = read.query_name
    read_str = read.query_alignment_sequence

    # write files
    ref_loc = os.path.join(workdir, "{}_ref.fa".format(read_id))
    read_loc = os.path.join(workdir, "{}_read.fa".format(read_id))
    aln_loc = os.path.join(workdir, '{}_out.txt'.format(read_id))
    with open(ref_loc, 'w') as ref_out:
        ref_out.write(">{}_{}_REF\n{}\n".format(contig, read_id, ref_str))
    with open(read_loc, 'w') as read_out:
        read_out.write(">{}\n{}\n".format(read_id, read_str))

    # generate command
    lastz_args = "--format=cigar --ambiguous=iupac"
    realign_args = "--rescoreByPosteriorProb"
    lastz_command = "%s %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace]"\
                   % (args.lastz, lastz_args, ref_loc, read_loc)
    realign_command = "%s %s --outputAllPosteriorProbs %s %s %s"\
                     % (args.realign, realign_args, aln_loc, ref_loc, read_loc)
    cmd = "%s | %s" % (lastz_command, realign_command)

    # run realign command
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=-1)
    output, _ = process.communicate()
    sts = process.wait()
    if sts != 0: raise RuntimeError( "Command exited with non-zero status %i: %s" % (sts, cmd))

    # get alignment
    success = True
    if not os.path.isfile(aln_loc):
        log("\tAlignment output not created: {}".format(aln_loc))
        success = False
    else:
        nucleotide_probs = calculate_nucleotide_probs(aln_loc, read_str, args)

        # write it out
        out_loc = os.path.join(args.out, "{}.tsv".format(read_id))
        metadata = {META_READ_ID: read_id, META_ALIGN_FILE: alignment_file, META_REFERNCE_FILE: args.ref,
                    META_CONTIG: contig, META_PECAN_CMD: cmd, META_DATE_GENERATED: str(datetime.now())}
        write_probabilities_out(out_loc, nucleotide_probs, contig, ref_start, metadata, args)
        assert os.path.isfile(out_loc), "Output file {} not created".format(out_loc)

    # cleanup
    if not args.keep_temp:
        os.remove(ref_loc)
        os.remove(read_loc)
        os.remove(aln_loc)

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


def main():
    args = parse_args()
    assert os.path.isfile(args.ref), "--ref argument does not exist: {}".format(args.ref)
    if not os.path.isdir(args.out): os.mkdir(args.out)
    assert os.path.isdir(args.out), "--output_directory argument could not be found/created: " \
                                                "{}".format(args.out)
    if not os.path.isdir(args.workdir): os.mkdir(args.workdir)
    assert os.path.isdir(args.workdir), "--workdir_directory argument could not be found/created: " \
                                                "{}".format(args.workdir)

    alignment_files = glob.glob(args.aln)
    assert len(alignment_files) > 0, "Could not find files matching {}".format(args.aln)

    # load reference
    reference_map = get_reference_map(args.ref)

    # generate alignments
    for alignment_file in alignment_files:
        realign_alignments(alignment_file, reference_map, args)

    # cleanup
    if len(os.listdir(args.workdir)) == 0:
        log("Removing empty workdir")
        os.rmdir(args.workdir)
    log("Fin.")



if __name__ == "__main__":
    main()