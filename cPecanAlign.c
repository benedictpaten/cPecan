#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "commonC.h"

static void usage(char *argv[]) {
    fprintf(stderr, "%s fasta_target fasta_query\n", argv[0]);
}

// Returns a hash mapping from sequence header to sequence data.
static stHash *readFastaFile(char *filename) {
    FILE *fasta = fopen(filename, "r");
    if (fasta == NULL) {
        st_errnoAbort("Could not open fasta file %s", filename);
    }
    stHash *headerToData = stHash_construct3(stHash_stringKey,
                                             stHash_stringEqualKey,
                                             free,
                                             free);
    struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
    fastaRead(fasta, seqs, seqLengths, headers);

    for (int64_t i = 0; i < seqs->length; i++) {
        char *fullHeader = headers->list[i];
        stList *headerTokens = stString_splitByString(fullHeader, " ");
        char *usableHeader = stString_copy(stList_get(headerTokens, 0));
        stHash_insert(headerToData, usableHeader, seqs->list[i]);
        stList_destruct(headerTokens);
    }
    destructList(seqs);
    destructList(seqLengths);
    destructList(headers);

    return headerToData;
}

// copied from cPecanRealign, which is sloppy.
void *convertToAnchorPair(void *aPair, void *extraArg) {
    stIntTuple *i = stIntTuple_construct2(stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2));
    stIntTuple_destruct(aPair);
    return i;
}

// copied from cPecanRealign
struct PairwiseAlignment *convertAlignedPairsToPairwiseAlignment(char *seqName1, char *seqName2, double score,
        int64_t length1, int64_t length2, stList *alignedPairs) {
    //Make pairwise alignment
    int64_t pX = -1, pY = -1, mL = 0;
    //Create an end matched pair, which is used to ensure the alignment has the correct end indels.
    struct List *opList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
    stList_append(alignedPairs, stIntTuple_construct2(length1, length2));
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(alignedPair, 0);
        int64_t y = stIntTuple_get(alignedPair, 1);
        assert(x - pX > 0);
        assert(y - pY > 0);
        if (x - pX > 0 && y - pY > 0) { //This is a hack for filtering
            if (x - pX > 1) { //There is an indel.
                if (mL > 0) {
                    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
                    mL = 0;
                }
                listAppend(opList, constructAlignmentOperation(PAIRWISE_INDEL_X, x - pX - 1, 0));
            }
            if (y - pY > 1) {
                if (mL > 0) {
                    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
                    mL = 0;
                }
                listAppend(opList, constructAlignmentOperation(PAIRWISE_INDEL_Y, y - pY - 1, 0));
            }
            mL++;
            pX = x;
            pY = y;
        }
    }
    //Deal with a trailing match, but exclude the final match
    if (mL > 1) {
        listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL - 1, 0));
    }
    stIntTuple_destruct(stList_pop(alignedPairs));
    //Construct the alignment
    struct PairwiseAlignment *pA = constructPairwiseAlignment(seqName1, 0, length1, 1, seqName2, 0, length2, 1, score,
            opList);
    return pA;
}

int main(int argc, char *argv[]) {
    // Parse arguments
    if (argc != 3) {
        usage(argv);
        return 1;
    }

    // You would load a custom HMM here if you wanted using
    // hmm_getStateMachine (see the realign code)
    StateMachine *stateMachine  = stateMachine5_construct(fiveState);

    PairwiseAlignmentParameters *parameters = pairwiseAlignmentBandingParameters_construct();

    stHash *targetSequences = readFastaFile(argv[1]);
    stHash *querySequences = readFastaFile(argv[2]);

    // For each query sequence, align it against all target sequences.
    stHashIterator *queryIt = stHash_getIterator(querySequences);
    char *queryHeader;
    while ((queryHeader = stHash_getNext(queryIt)) != NULL) {
        char *querySeq = stHash_search(querySequences, queryHeader);
        stHashIterator *targetIt = stHash_getIterator(targetSequences);
        char *targetHeader;
        while ((targetHeader = stHash_getNext(targetIt)) != NULL) {
            char *targetSeq = stHash_search(targetSequences, targetHeader);
            // Here we should try both the target sequence and its
            // reverse-complemented version


            // Aligns the sequences.
            // If you have alignment constraints (anchors) you should
            // replace this with getAlignedPairsUsingAnchors.
            stList *alignedPairs = getAlignedPairs(stateMachine, targetSeq,
                                                   querySeq, parameters,
                                                   true, true);
            // Takes into account the probability of aligning to a
            // gap, by transforming the posterior probability into the
            // AMAP objective function (see Schwartz & Pachter, 2007).
            alignedPairs = reweightAlignedPairs2(alignedPairs, strlen(targetSeq),
                                                 strlen(querySeq),
                                                 parameters->gapGamma);
            // I think this calculates the optimal ordered set of
            // alignments from the unordered set of aligned pairs, not
            // completely sure.
            alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs,
                                                                     targetSeq,
                                                                     querySeq,
                                                                     // This parameter says that the minimum posterior probability we will accept has to be at least 0.9.
                                                                     0.9);

            // After this the "aligned pairs" data structure changes,
            // which is a little sketchy. It's just so that the
            // alignment can be printed properly.
            stList_mapReplace(alignedPairs, convertToAnchorPair, NULL);
            stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
            struct PairwiseAlignment *alignment = convertAlignedPairsToPairwiseAlignment(targetHeader, queryHeader,
                                                                                  0, strlen(targetSeq), strlen(querySeq), alignedPairs);
            // Output the cigar string
            cigarWrite(stdout, alignment, 0);

            stList_destruct(alignedPairs);
            destructPairwiseAlignment(alignment);
        }
        stHash_destructIterator(targetIt);
    }
    stHash_destructIterator(queryIt);

    // Clean up
    stHash_destruct(targetSequences);
    stHash_destruct(querySequences);

    pairwiseAlignmentBandingParameters_destruct(parameters);
    stateMachine_destruct(stateMachine);
}
