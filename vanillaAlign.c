#include <getopt.h>
#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "emissionMatrix.h"
#include "stateMachine.h"
#include "nanopore.h"
#include "continuousHmm.h"


void usage() {
    fprintf(stderr, "vanillaAlign binary, meant to be used through the signalAlign program.\n");
    fprintf(stderr, "See doc for signalAlign for help\n");
}

void writePosteriorProbs(char *posteriorProbsFile, char *readFile, double *matchModel,
                         double *events, char *target,
                         int64_t eventSequenceOffset, int64_t referenceSequenceOffset,
                         stList *alignedPairs, Strand strand) {
    // label for tsv output
    char *strandLabel;
    if (strand == template) {
        strandLabel = "t";
    }
    if (strand == complement) {
        strandLabel = "c";
    }

    FILE *fH = fopen(posteriorProbsFile, "a");

    st_uglyf("refoffset ---> %lld\n", referenceSequenceOffset);

    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        // unpack it to make this easier
        int64_t x = stIntTuple_get(aPair, 1) + referenceSequenceOffset;         // target index
        int64_t y = stIntTuple_get(aPair, 2) + eventSequenceOffset;             // event index
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;  // posterior prob

        // get the observations from the events
        double eventMean = events[(y * NB_EVENT_PARAMS)];
        double eventNoise = events[(y * NB_EVENT_PARAMS) + 1];
        double eventDuration = events[(y * NB_EVENT_PARAMS) + 2];

        // make the kmer string at the target index,
        char *k_i = st_malloc(KMER_LENGTH * sizeof(char));
        for (int64_t k = 0; k < KMER_LENGTH; k++) {
            k_i[k] = *(target + (x + k));
        }
        // get the kmer index
        int64_t targetKmerIndex = emissions_discrete_getKmerIndexFromKmer(k_i);

        // get the expected event mean amplitude and noise
        double E_levelu = matchModel[1 + (targetKmerIndex * MODEL_PARAMS)];
        double E_noiseu = matchModel[1 + (targetKmerIndex * MODEL_PARAMS + 2)];

        // write to disk
        fprintf(fH, "%lld\t%lld\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, k_i, readFile, strandLabel, eventMean,
                eventNoise, eventDuration, p, E_levelu, E_noiseu);

        // cleanup
        free(k_i);
    }
    fclose(fH);
}

stList *getRemappedAnchorPairs(stList *unmappedAnchors, int64_t *eventMap, int64_t mapOffset) {
    stList *remapedAnchors = nanopore_remapAnchorPairsWithOffset(unmappedAnchors, eventMap, mapOffset);

    stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);
    return filteredRemappedAnchors;
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type,
                                Strand strand) {
    if ((type != threeState) && (type != vanilla) && (type != echelon) && (type != fourState)) {
        st_errAbort("vanillaAlign - incompatable stateMachine type request");
    }

    if (type == vanilla) {
        StateMachine *sM = getSignalStateMachine3Vanilla(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        stateMachine3Vanilla_setStrandTransitionsToDefaults(sM, strand);
        return sM;
    }
    if (type == threeState) {
        StateMachine *sM = getStrawManStateMachine3(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    if (type == fourState) {
        StateMachine *sM = getStateMachine4(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    if (type == echelon) {
        StateMachine *sM = getStateMachineEchelon(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    return 0; // else return 1 to check?
}

void loadHmmRoutine(const char *hmmFile, StateMachine *sM, StateMachineType type) {
    Hmm *hmm = hmmContinuous_loadSignalHmm(hmmFile, type);
    hmmContinuous_loadExpectations(sM, hmm, type);
    hmmContinuous_destruct(hmm, type);
}

static double totalScore(stList *alignedPairs) {
    double score = 0.0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        score += stIntTuple_get(aPair, 0);
    }
    return score;
}

double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs) {
    /*
     * Gives the average posterior match probability per base of the two sequences, ignoring indels.
     */
    return 100.0 * totalScore(alignedPairs) / ((double) stList_length(alignedPairs) * PAIR_ALIGNMENT_PROB_1);
}

stList *performSignalAlignmentP(StateMachine *sM, Sequence *sY, int64_t *eventMap, int64_t mapOffset, char *target,
                                PairwiseAlignmentParameters *p, stList *unmappedAnchors,
                                void *(*targetGetFcn)(void *, int64_t),
                                void (*posteriorProbFcn)(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                                         DpMatrix *backwardDpMatrix, Sequence* sX, Sequence* sY,
                                                         double totalProbability, PairwiseAlignmentParameters *p,
                                                         void *extraArgs),
                                bool banded) {
    int64_t lX = sequence_correctSeqLength(strlen(target), event);

    if (banded) {
        fprintf(stderr, "vanillaAlign - doing banded alignment\n");

        // remap anchor pairs
        stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

        // make sequences
        Sequence *sX = sequence_construct2(lX, target, targetGetFcn, sequence_sliceNucleotideSequence2);

        if (sM->type == echelon) {
            sequence_padSequence(sX);
        }

        // do alignment
        stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, filteredRemappedAnchors, p,
                                                           posteriorProbFcn, 1, 1);
        return alignedPairs;
    } else {
        fprintf(stderr, "vanillaAlign - doing non-banded alignment\n");

        stList *alignedPairs = getAlignedPairsWithoutBanding(sM, target, sY->elements, lX, sY->length, p, targetGetFcn,
                                                             sequence_getEvent, posteriorProbFcn, 1, 1);
        return alignedPairs;
    }
}

stList *performSignalAlignment(StateMachine *sM, const char *hmmFile, Sequence *eventSequence, int64_t *eventMap,
                               int64_t mapOffset,
                               char *target, PairwiseAlignmentParameters *p, stList *unmappedAncors, bool banded) {
    if ((sM->type != threeState) && (sM->type != vanilla) && (sM->type != echelon) && (sM->type != fourState)) {
        st_errAbort("vanillaAlign - You're trying to do the wrong king of alignment");
    }

    // load HMM if given
    if (hmmFile != NULL) {
        fprintf(stderr, "loading HMM from file, %s\n", hmmFile);
        loadHmmRoutine(hmmFile, sM, sM->type);
    }

    // decision tree for different stateMachine types
    if ((sM->type == vanilla) || (sM->type == echelon)) {
        if (sM->type == vanilla) {
            stList *alignedPairs = performSignalAlignmentP(sM, eventSequence, eventMap, mapOffset,
                                                           target, p, unmappedAncors, sequence_getKmer2,
                                                           diagonalCalculationPosteriorMatchProbs, banded);
            return alignedPairs;
        } else {
            stList *alignedPairs = performSignalAlignmentP(sM, eventSequence, eventMap, mapOffset,
                                                           target, p, unmappedAncors, sequence_getKmer2,
                                                           diagonalCalculationMultiPosteriorMatchProbs, banded);
            return alignedPairs;
        }
    }
    if ((sM->type == threeState) || (sM->type == fourState)) {
        stList *alignedPairs = performSignalAlignmentP(sM, eventSequence, eventMap, mapOffset, target, p,
                                                       unmappedAncors, sequence_getKmer,
                                                       diagonalCalculationPosteriorMatchProbs, banded);
        return alignedPairs;
    }
    return 0;
}

char *getSubSequence(char *seq, int64_t start, int64_t end, bool strand) {
    if (strand) {
        seq = stString_getSubString(seq, start, end - start);
        return seq;
    }
    seq = stString_getSubString(seq, end, start - end);
    return seq;
}

void rebasePairwiseAlignmentCoordinates(int64_t *start, int64_t *end, int64_t *strand,
                                        int64_t coordinateShift, bool flipStrand) {
    *start += coordinateShift;
    *end += coordinateShift;
    if (flipStrand) {
        *strand = *strand ? 0 : 1;
        int64_t i = *end;
        *end = *start;
        *start = i;
    }
}

stList *guideAlignmentToRebasedAnchorPairs(struct PairwiseAlignment *pA, PairwiseAlignmentParameters *p) {
    // check if we need to flip the reference
    bool flipStrand1 = !pA->strand1;
    int64_t refCoordShift = (pA->strand1 ? pA->start1 : pA->end1);

    // rebase the reference alignment to (0), but not the nanopore read, this is corrected when remapping the
    // anchorPairs
    rebasePairwiseAlignmentCoordinates(&(pA->start1), &(pA->end1), &(pA->strand1), -refCoordShift, flipStrand1);
    checkPairwiseAlignment(pA);

    //Convert input alignment into anchor pairs
    stList *unfilteredAnchorPairs = convertPairwiseForwardStrandAlignmentToAnchorPairs(
            pA, p->constraintDiagonalTrim);

    // sort
    stList_sort(unfilteredAnchorPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);

    // filter
    stList *anchorPairs = filterToRemoveOverlap(unfilteredAnchorPairs);

    return anchorPairs;
}

Sequence *makeEventSequenceFromPairwiseAlignment(double *events, int64_t queryStart, int64_t queryEnd,
                                                 int64_t *eventMap) {
    // find the event mapped to the start and end of the 2D read alignment
    int64_t startIdx = eventMap[queryStart];
    int64_t endIdx = eventMap[queryEnd];

    // move the event pointer to the first event
    size_t elementSize = sizeof(double);
    void *elements = (char *)events + ((startIdx * NB_EVENT_PARAMS) * elementSize);

    // make the eventSequence
    Sequence *eventS = sequence_construct2(endIdx - startIdx, elements, sequence_getEvent,
                                           sequence_sliceEventSequence2);

    return eventS;
}

void getSignalExpectations(const char *model, const char *inputHmm, Hmm *hmmExpectations, StateMachineType type,
                           NanoporeReadAdjustmentParameters npp, Sequence *eventSequence,
                           int64_t *eventMap, int64_t mapOffset, char *trainingTarget, PairwiseAlignmentParameters *p,
                           stList *unmappedAnchors, Strand strand) {
    // load match model, build stateMachine
    StateMachine *sM = buildStateMachine(model, npp, type, strand);

    // load HMM transitions, if needed
    if (inputHmm != NULL) {
        fprintf(stderr, "vanillaAlign - loading HMM from file, %s\n", inputHmm);
        loadHmmRoutine(inputHmm, sM, type);
    }

    // correct sequence length
    int64_t lX = sequence_correctSeqLength(strlen(trainingTarget), event);

    // remap the anchors
    stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap, mapOffset);

    // make sequence objects, seperate the target sequences based on HMM type, also implant the match model if we're
    // using a conditional model
    //Sequence *eventS = sequence_construct2(nbEvents, events, sequence_getEvent, sequence_sliceEventSequence2);

    if (type == vanilla) {
        Sequence *target = sequence_construct2(lX, trainingTarget, sequence_getKmer2,
                                               sequence_sliceNucleotideSequence2);
        vanillaHmm_implantMatchModelsintoHmm(sM, hmmExpectations);
        // get expectations
        getExpectationsUsingAnchors(sM, hmmExpectations, target, eventSequence, filteredRemappedAnchors, p,
                                    diagonalCalculation_signal_Expectations, 1, 1);
    } else {
        Sequence *target = sequence_construct2(lX, trainingTarget, sequence_getKmer,
                                               sequence_sliceNucleotideSequence2);
        // get expectations
        getExpectationsUsingAnchors(sM, hmmExpectations, target, eventSequence, filteredRemappedAnchors, p,
                                    diagonalCalculation_signal_Expectations, 1, 1);
    }
}

void printPairwiseAlignmentSummary(struct PairwiseAlignment *pA) {
    st_uglyf("contig 1: %s\n", pA->contig1);
    st_uglyf("strand 1: %lld\n", pA->strand1);
    st_uglyf("start  1: %lld\n", pA->start1);
    st_uglyf("end    1: %lld\n", pA->end1);
    st_uglyf("contig 2: %s\n", pA->contig2);
    st_uglyf("strand 2: %lld\n", pA->strand2);
    st_uglyf("start  2: %lld\n", pA->start2);
    st_uglyf("end    2: %lld\n", pA->end2);
}

int main(int argc, char *argv[]) {
    StateMachineType sMtype = vanilla;
    bool banded = FALSE;
    int64_t j;
    int64_t diagExpansion = 20;
    double threshold = 0.01;
    int64_t constraintTrim = 14;
    char *templateModelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *complementModelFile = stString_print("../../cPecan/models/complement_median68pA_pop2.model");
    char *readLabel = NULL;
    char *npReadFile = NULL;
    char *targetFile = NULL;
    char *posteriorProbsFile = NULL;
    char *templateHmmFile = NULL;
    char *complementHmmFile = NULL;
    char *templateExpectationsFile = NULL;
    char *complementExpectationsFile = NULL;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"strawMan",                no_argument,        0,  's'},
                {"fourState",               no_argument,        0,  'f'},
                {"echelon",                 no_argument,        0,  'e'},
                {"banded",                  no_argument,        0,  'b'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"readLabel",               required_argument,  0,  'L'},
                {"npRead",                  required_argument,  0,  'q'},
                {"reference",               required_argument,  0,  'r'},
                {"posteriors",              required_argument,  0,  'u'},
                {"inTemplateHmm",           required_argument,  0,  'y'},
                {"inComplementHmm",         required_argument,  0,  'z'},
                {"templateExpectations",    required_argument,  0,  't'},
                {"complementExpectations",  required_argument,  0,  'c'},
                {"diagonalExpansion",       required_argument,  0,  'x'},
                {"threshold",               required_argument,  0,  'd'},
                {"constraintTrim",          required_argument,  0,  'm'},

                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:s:f:e:b:T:C:L:q:r:u:y:z:t:c:i:x:d:m:", long_options, &option_index);

        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 0;
            case 's':
                sMtype = threeState;
                break;
            case 'f':
                sMtype = fourState;
                break;
            case 'e':
                sMtype = echelon;
                break;
            case 'b':
                banded = TRUE;
                break;
            case 'T':
                templateModelFile = stString_copy(optarg);
                break;
            case 'C':
                complementModelFile = stString_copy(optarg);
                break;
            case 'L':
                readLabel = stString_copy(optarg);
                break;
            case 'q':
                npReadFile = stString_copy(optarg);
                break;
            case 'r':
                targetFile = stString_copy(optarg);
                break;
            case 'u':
                posteriorProbsFile = stString_copy(optarg);
                break;
            case 't':
                templateExpectationsFile = stString_copy(optarg);
                break;
            case 'c':
                complementExpectationsFile = stString_copy(optarg);
                break;
            case 'y':
                templateHmmFile = stString_copy(optarg);
                break;
            case 'z':
                complementHmmFile = stString_copy(optarg);
                break;
            case 'x':
                j = sscanf(optarg, "%" PRIi64 "", &diagExpansion);
                assert (j == 1);
                assert (diagExpansion >= 0);
                diagExpansion = (int64_t)diagExpansion;
                break;
            case 'd':
                j = sscanf(optarg, "%lf", &threshold);
                assert (j == 1);
                assert (threshold >= 0);
                break;
            case 'm':
                j = sscanf(optarg, "%" PRIi64 "", &constraintTrim);
                assert (j == 1);
                assert (constraintTrim >= 0);
                constraintTrim = (int64_t)constraintTrim;
                break;
            default:
                usage();
                return 1;
        }
    }

    if (sMtype == threeState) {
        fprintf(stderr, "vanillaAlign - using strawMan model\n");
    }
    if (sMtype == vanilla) {
        fprintf(stderr, "vanillaAlign - using vanilla model\n");
    }
    if (sMtype == fourState) {
        fprintf(stderr, "vanillaAlign - using four-state PairHMM model\n");
    }
    if (sMtype == echelon) {
        fprintf(stderr, "vanillaAlign - using echelon model\n");
    }

    // load reference sequence (reference sequence)
    FILE *reference = fopen(targetFile, "r");
    char *referenceSequence = stFile_getLineFromFile(reference);

    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // make some params
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = threshold;
    p->constraintDiagonalTrim = constraintTrim;
    p->diagonalExpansion = diagExpansion;

    // get pairwise alignment from stdin, in exonerate CIGAR format
    FILE *fileHandleIn = stdin;

    // parse input
    struct PairwiseAlignment *pA;
    pA = cigarRead(fileHandleIn);

    // todo put in if debug:
    //printPairwiseAlignmentSummary(pA);

    // slice out the section of the reference we're aligning to
    char *trimmedRefSeq = getSubSequence(referenceSequence, pA->start1, pA->end1, pA->strand1);

    // reverse complement for complement event sequence
    char *rc_trimmedRefSeq = stString_reverseComplementString(trimmedRefSeq);

    // constrain the event sequence to the positions given by the guide alignment
    Sequence *tEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->templateEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->templateEventMap);

    Sequence *cEventSequence = makeEventSequenceFromPairwiseAlignment(npRead->complementEvents,
                                                                      pA->start2, pA->end2,
                                                                      npRead->complementEventMap);


    // the alignments start at (0,0) so we need to correct them based on the guide alignment later.
    // record this shift here
    int64_t tCoordinateShift = npRead->templateEventMap[pA->start2];
    int64_t cCoordinateShift = npRead->complementEventMap[pA->start2];

    // orig
    int64_t rCoordinateShift_t = (pA->strand1 ? pA->start1 : pA->end1);
    int64_t rCoordinateShift_c = (int64_t)strlen(referenceSequence) - (pA->strand1 ? pA->end1 : pA->start1);

    // next try works for ecoli not for zymo
    //int64_t rCoordinateShift_t = (pA->strand1 ? pA->start1 : (int64_t)strlen(referenceSequence) - pA->start1);
    //int64_t rCoordinateShift_c = (pA->strand1 ? (int64_t)strlen(referenceSequence) - pA->end1 : pA->end1);

    //int64_t rCoordinateShift_t = (pA->strand1 ? pA->start1 : (int64_t)strlen(referenceSequence) - pA->end1);
    //int64_t rCoordinateShift_c = (pA->strand1 ? (int64_t)strlen(referenceSequence) - pA->end1 : pA->start1);

    st_uglyf("tOffset: %lld, cOffset %lld \n", rCoordinateShift_t, rCoordinateShift_c);

    stList *anchorPairs = guideAlignmentToRebasedAnchorPairs(pA, p);

    if ((templateExpectationsFile != NULL) && (complementExpectationsFile != NULL)) {
        // Expectation Routine //

        if ((sMtype != threeState) && (sMtype != vanilla)) {
            st_errAbort("vanillaAlign - getting expectations not allowed for this HMM type, yet");
        }

        // make empty HMM to collect expectations
        Hmm *templateExpectations = hmmContinuous_getEmptyHmm(sMtype, 0.0001);
        Hmm *complementExpectations = hmmContinuous_getEmptyHmm(sMtype, 0.0001);

        // get expectations for template
        fprintf(stderr, "vanillaAlign - getting expectations for template\n");
        getSignalExpectations(templateModelFile, templateHmmFile, templateExpectations, sMtype, npRead->templateParams,
                              tEventSequence, npRead->templateEventMap, pA->start2,
                              trimmedRefSeq, p,
                              anchorPairs, template);

        // write to file
        fprintf(stderr, "vanillaAlign - writing expectations to file: %s\n\n", templateExpectationsFile);
        hmmContinuous_writeToFile(templateExpectationsFile, templateExpectations, sMtype);

        // get expectations for the complement
        fprintf(stderr, "vanillaAlign - getting expectations for complement\n");
        getSignalExpectations(complementModelFile, complementHmmFile, complementExpectations, sMtype,
                              npRead->complementParams, cEventSequence, npRead->complementEventMap, pA->start2,
                              rc_trimmedRefSeq, p, anchorPairs, complement);

        // write to file
        fprintf(stderr, "vanillaAlign - writing expectations to file: %s\n\n", complementExpectationsFile);
        hmmContinuous_writeToFile(complementExpectationsFile, complementExpectations, sMtype);
        // todo make hmm destruct, test
        return 0;
    } else {
        // Alignment Procedure //

        // Template alignment
        fprintf(stderr, "vanillaAlign - starting template alignment\n");

        // make template stateMachine
        StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype, template);

        // get aligned pairs
        stList *templateAlignedPairs = performSignalAlignment(sMt, templateHmmFile, tEventSequence,
                                                              npRead->templateEventMap, pA->start2, trimmedRefSeq,
                                                              p, anchorPairs, banded);

        double templatePosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(templateAlignedPairs);

        fprintf(stdout, "%s %lld\t%lld(%f)\t", readLabel, stList_length(anchorPairs),
                stList_length(templateAlignedPairs), templatePosteriorScore);

        // sort
        stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            writePosteriorProbs(posteriorProbsFile, readLabel, sMt->EMISSION_MATCH_PROBS,
                                npRead->templateEvents, referenceSequence, // use non-trimmed objects here
                                tCoordinateShift, rCoordinateShift_t,        // correction is in function
                                templateAlignedPairs, template);
        }

        // intermediate clean up
        stateMachine_destruct(sMt);
        sequence_sequenceDestroy(tEventSequence);
        stList_destruct(templateAlignedPairs);

        // Complement alignment
        fprintf(stderr, "vanillaAlign - starting complement alignment\n");
        StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype, complement);

        // get aligned pairs
        stList *complementAlignedPairs = performSignalAlignment(sMc, complementHmmFile, cEventSequence,
                                                                npRead->complementEventMap, pA->start2,
                                                                rc_trimmedRefSeq, p, anchorPairs, banded);

        double complementPosteriorScore = scoreByPosteriorProbabilityIgnoringGaps(complementAlignedPairs);

        fprintf(stdout, "%lld(%f)\n", stList_length(complementAlignedPairs), complementPosteriorScore);

        // sort
        stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing

        // write to file
        if (posteriorProbsFile != NULL) {
            writePosteriorProbs(posteriorProbsFile, readLabel, sMc->EMISSION_MATCH_PROBS,
                                npRead->complementEvents, stString_reverseComplementString(referenceSequence),
                                cCoordinateShift, rCoordinateShift_c,  // see above
                                complementAlignedPairs, complement);
        }

        // final alignment clean up
        stateMachine_destruct(sMc);
        sequence_sequenceDestroy(cEventSequence);
        stList_destruct(complementAlignedPairs);
        fprintf(stderr, "vanillaAlign - SUCCESS: finished alignment of query %s, exiting\n", readLabel);
    }

    return 0;
}
