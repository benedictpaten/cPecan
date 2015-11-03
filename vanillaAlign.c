#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "nanopore.h"

typedef enum _strand {
    template = 0,
    complement = 1
} Strand;

void writePosteriorProbs(char *posteriorProbsFile, stList *alignedPairs, Strand strand) {
    /*
     * Writes the posterior match probabibilities to a tab separated file, each line being X coordinate, Y coordinate, Match probability
     */
    char *strandLabel;
    if (strand == template) {
        strandLabel = "t";
    }
    if (strand == complement) {
        strandLabel = "c";
    }

    FILE *fH = fopen(posteriorProbsFile, "a");
    for(int64_t i=0;i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        fprintf(fH, "%s " "%" PRIi64 "\t%" PRIi64 "\t%f\n", strandLabel, stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2), ((double)stIntTuple_get(aPair, 0))/PAIR_ALIGNMENT_PROB_1);
    }
    fclose(fH);
}

stList *getRemappedAnchorPairs(char *target, char *query, int64_t *eventMap, PairwiseAlignmentParameters *p) {
    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(target, query, p);

    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, eventMap);

    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);
    // destruct remapped pairs?
    return filteredRemappedAnchors;
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp) {
    StateMachine *sM = getSignalStateMachine3Vanilla(modelFile);
    emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
    return sM;
}

stList *performSignalAlignment(StateMachine *sM, double *events, int64_t nbEvents, int64_t *eventMap, char *target,
                               PairwiseAlignmentParameters *p, stList *unmappedAnchors) {
    // remap anchor pairs
    stList *remapedAnchors = nanopore_remapAnchorPairs(unmappedAnchors, eventMap);

    stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);

    // make sequences
    int64_t lX = sequence_correctSeqLength(strlen(target), event);
    Sequence *sX = sequence_construct(lX, target, sequence_getKmer2);
    Sequence *sY = sequence_construct(nbEvents, events, sequence_getEvent);

    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       1, 1);

    return alignedPairs;
}

int main(int argc, char *argv[]) {
    // usage:
    // ./vanillaAlign target.txt read.npRead posteriorProbs.tsv

    // load target sequence (reference sequence)
    FILE *targetFile = fopen(argv[1], "r");
    char *targetSeq = stFile_getLineFromFile(targetFile);
    char *rc_targetSeq = stString_reverseComplementString(targetSeq);

    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(argv[2]);

    // make some params
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // get anchors
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(targetSeq, npRead->twoDread, p);

    //_____template Strand alignment_______//
    // get anchors
    char *templateModelFile = stString_print("../../cPecan/models/template_median68pA.model");
    StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams);
    stList *templateAlignedPairs = performSignalAlignment(sMt, npRead->templateEvents, npRead->nbTemplateEvents,
                                                          npRead->templateEventMap, targetSeq, p,
                                                          anchorPairs);
    writePosteriorProbs(argv[3], templateAlignedPairs, template);
    // intermediate cleanup
    //stateMachine_destruct(sMt);

    //_____complement Strand alignment_______//
    char *complementModelFile = stString_print("../../cPecan/models/complement_median68pA_pop2.model");
    StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams);
    stList *complementAlignedPairs = performSignalAlignment(sMc, npRead->complementEvents, npRead->nbComplementEvents,
                                                            npRead->complementEventMap, rc_targetSeq, p,
                                                            anchorPairs);
    //stateMachine_destruct(sMc);
    writePosteriorProbs(argv[3], complementAlignedPairs, complement);

    return 0;
}
