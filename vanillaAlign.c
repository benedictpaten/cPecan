#include <getopt.h>
#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "nanopore.h"
#include "continuousHmm.h"
typedef enum _strand {
    template = 0,
    complement = 1
} Strand;

void usage() {
    fprintf(stderr, "./vanillaAlign target.txt read.npRead posteriorProbs.tsv\n");
}

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

stList *getRemappedAnchorPairs(stList *unmappedAnchors, int64_t *eventMap) {
    // get anchors using lastz
    stList *remapedAnchors = nanopore_remapAnchorPairs(unmappedAnchors, eventMap);

    stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);
    return filteredRemappedAnchors;
}

StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type) {
    assert((type == threeState) || (type == vanilla));
    if (type == vanilla) {
        StateMachine *sM = getSignalStateMachine3Vanilla(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    if (type == threeState) {
        StateMachine *sM = getStrawManStateMachine3(modelFile);
        emissions_signal_scaleModel(sM, npp.scale, npp.shift, npp.var, npp.scale_sd, npp.var_sd);
        return sM;
    }
    return 0;
}

void loadHmmRoutine(const char *hmmFile, StateMachine *sM, StateMachineType type) {
    Hmm *hmm = hmmContinuous_loadSignalHmm(hmmFile, type);
    hmmContinuous_loadExpectations(sM, hmm, type);
    hmmContinuous_destruct(hmm, type);
}

stList *performSignalAlignmentP(StateMachine *sM, double *events, int64_t nbEvents, int64_t *eventMap, char *target,
                                PairwiseAlignmentParameters *p, stList *unmappedAnchors,
                                void *(*targetGetFcn)(void *, int64_t)) {
    // remap anchor pairs
    //stList *remapedAnchors = nanopore_remapAnchorPairs(unmappedAnchors, eventMap);
    //stList *filteredRemappedAnchors = filterToRemoveOverlap(remapedAnchors);
    stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap);


    // make sequences
    int64_t lX = sequence_correctSeqLength(strlen(target), event);
    Sequence *sX = sequence_construct(lX, target, targetGetFcn);
    Sequence *sY = sequence_construct(nbEvents, events, sequence_getEvent);

    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       1, 1);

    return alignedPairs;
}

stList *performSignalAlignment(const char *model, const char *hmmFile, StateMachineType type,
                               NanoporeReadAdjustmentParameters npp, double *events, int64_t nbEvents,
                               int64_t *eventMap, char *target, PairwiseAlignmentParameters *p,
                               stList *unmappedAncors) {
    if ((type != threeState) && (type != vanilla)) {
        st_errAbort("you're trying to do the wrong kind of alignment\n");
    }

    // load model into stateMachine
    StateMachine *sM = buildStateMachine(model, npp, type);

    // load HMM if given
    if (hmmFile != NULL) {
        fprintf(stderr, "loading HMM from file, %s\n", hmmFile);
        loadHmmRoutine(hmmFile, sM, type);
    }

    // do alignment
    if (type == vanilla) {
        stList *alignedPairs = performSignalAlignmentP(sM, events, nbEvents, eventMap, target, p, unmappedAncors,
                                                       sequence_getKmer2);
        return alignedPairs;
    }
    if (type == threeState) {
        stList *alignedPairs = performSignalAlignmentP(sM, events, nbEvents, eventMap, target, p, unmappedAncors,
                                                       sequence_getKmer);
        return alignedPairs;
    }
    return 0;
}

Hmm *performBaumWelchTrainingP(const char *model, const char *inputHmm, StateMachineType type,
                              NanoporeReadAdjustmentParameters npp, double *events, int64_t nbEvents,
                              int64_t *eventMap, char *trainingTarget, PairwiseAlignmentParameters *p,
                              stList *unmappedAnchors, int64_t iterations, void *(*getFcn)(void *, int64_t)) {
    // load model into stateMachine
    StateMachine *sM = buildStateMachine(model, npp, type);

    // load HMM if given
    if (inputHmm != NULL) {
        fprintf(stderr, "vanillaAlign - loading HMM from file, %s\n", inputHmm);
        loadHmmRoutine(inputHmm, sM, type);
    }

    // correct sequence length, once
    int64_t lX = sequence_correctSeqLength(strlen(trainingTarget), event);

    // remap the anchors
    stList *filteredRemappedAnchors = getRemappedAnchorPairs(unmappedAnchors, eventMap);

    // initial setup
    double pLikelihood = -INFINITY;
    Hmm *hmm = hmmContinuous_getEmptyHmm(type);
    int64_t i = 0;

    // EM
    fprintf(stderr, "vanillaAlign - about to start EM, set to %lld iterations\n", iterations);
    fprintf(stderr, "vanillaAlign - likelihoods:\n");

    while ((pLikelihood < hmm->likelihood) || i <= iterations) {
        hmmContinuous_destruct(hmm, type);
        hmm = hmmContinuous_getEmptyHmm(type);
        // make sequence objects
        Sequence *target = sequence_construct(lX, trainingTarget, getFcn);
        Sequence *eventS = sequence_construct(nbEvents, events, sequence_getEvent);

        // implant if using vanilla
        if (type == vanilla) {
            vanillaHmm_implantMatchModelsintoHmm(sM, hmm);
        }

        // E Step
        getExpectationsUsingAnchors(sM, hmm, target, eventS, filteredRemappedAnchors,
                                    p, diagonalCalculation_signal_Expectations, 1, 1);

        // normalize
        hmmContinuous_normalize(hmm, type);

        fprintf(stderr, "%f ", hmm->likelihood);

        // M Step
        hmmContinuous_loadExpectations(sM, hmm, type);

        // check if we're done
        if (pLikelihood == hmm->likelihood) {
            break;
        }
        // update
        pLikelihood = hmm->likelihood;
        i++;
    }
    fprintf(stderr, "\nvanillaAlign - finished iterations\n\n");

    return hmm;
}

Hmm *performBaumWelchTraining(const char *model, const char *inputHmm, StateMachineType type,
                              NanoporeReadAdjustmentParameters npp, double *events, int64_t nbEvents,
                              int64_t *eventMap, char *trainingTarget, PairwiseAlignmentParameters *p,
                              stList *unmappedAnchors, int64_t iterations) {
    if (type == vanilla) {
        Hmm *hmm = performBaumWelchTrainingP(model, inputHmm, type, npp, events, nbEvents, eventMap, trainingTarget,
                                             p, unmappedAnchors, iterations, sequence_getKmer2);
        return hmm;
    }
    if (type == threeState) {
        Hmm *hmm = performBaumWelchTrainingP(model, inputHmm, type, npp, events, nbEvents, eventMap, trainingTarget,
                                             p, unmappedAnchors, iterations, sequence_getKmer);
        return hmm;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    StateMachineType sMtype = vanilla;
    int64_t j;
    int64_t iter = 0;
    char *templateModelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *complementModelFile = stString_print("../../cPecan/models/complement_median68pA_pop2.model");
    char *npReadFile = NULL;
    char *targetFile = NULL;
    char *posteriorProbsFile = NULL;
    char *templateHmmFile = NULL;
    char *complementHmmFile = NULL;
    char *templateTrainedHmmFile = NULL;
    char *complementTrainedHmmFile = NULL;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"strawMan",                no_argument,        0,  's'},
                {"templateModel",           required_argument,  0,  'T'},
                {"complementModel",         required_argument,  0,  'C'},
                {"npRead",                  required_argument,  0,  'q'},
                {"target",                  required_argument,  0,  'r'},
                {"posteriors",              required_argument,  0,  'u'},
                {"loadTemplateHmm",         required_argument,  0,  'y'},
                {"loadComplementHmm",       required_argument,  0,  'z'},
                {"templateTrainedHmm",      required_argument,  0,  't'},
                {"complementTrainedHmm",    required_argument,  0,  'c'},
                {"iterations",              required_argument,  0,  'i'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:s:T:C:q:r:u:y:z:t:c:i:", long_options, &option_index);

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
            case 'T':
                templateModelFile = stString_copy(optarg);
                break;
            case 'C':
                complementModelFile = stString_copy(optarg);
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
                templateTrainedHmmFile = stString_copy(optarg);
                break;
            case 'c':
                complementTrainedHmmFile = stString_copy(optarg);
                break;
            case 'y':
                templateHmmFile = stString_copy(optarg);
                break;
            case 'z':
                complementHmmFile = stString_copy(optarg);
                break;
            case 'i':
                j = sscanf(optarg, "%" PRIi64 "", &iter);
                assert(j == 1);
                assert(iter >= 0);
                iter = (int64_t)iter;
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

    // load target sequence (reference sequence)
    FILE *target = fopen(targetFile, "r");
    char *targetSeq = stFile_getLineFromFile(target);
    char *rc_targetSeq = stString_reverseComplementString(targetSeq);

    // load nanopore read
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    // make some params
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

    // get anchors
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(targetSeq, npRead->twoDread, p);

    // EM training routine //
    if ((templateTrainedHmmFile != NULL) || (complementTrainedHmmFile != NULL)) {
        if (templateTrainedHmmFile != NULL) {
            fprintf(stderr, "vanillaAlign - starting training for template hmm\n");
            Hmm *templateTrainedHmm = performBaumWelchTraining(templateModelFile, templateHmmFile, sMtype,
                                                               npRead->templateParams, npRead->templateEvents,
                                                               npRead->nbTemplateEvents, npRead->templateEventMap,
                                                               targetSeq, p, anchorPairs, iter);
            hmmContinuous_writeToFile(templateTrainedHmmFile, templateTrainedHmm, sMtype);
        }
        if (complementTrainedHmmFile != NULL) {
            fprintf(stderr, "vanillaAlign - starting training for complement hmm\n");
            Hmm *complementTrainedHmm = performBaumWelchTraining(complementModelFile, complementHmmFile, sMtype,
                                                                 npRead->complementParams, npRead->complementEvents,
                                                                 npRead->nbComplementEvents, npRead->complementEventMap,
                                                                 rc_targetSeq, p, anchorPairs, iter);
            hmmContinuous_writeToFile(complementTrainedHmmFile, complementTrainedHmm, sMtype);
        }
    } else {
        fprintf(stderr, "vanillaAlign - starting alignment\n");
        stList *templateAlignedPairs = performSignalAlignment(templateModelFile, templateHmmFile, sMtype,
                                                              npRead->templateParams, npRead->templateEvents,
                                                              npRead->nbTemplateEvents, npRead->templateEventMap,
                                                              targetSeq, p, anchorPairs);

        stList *complementAlignedPairs = performSignalAlignment(complementModelFile, complementHmmFile, sMtype,
                                                                npRead->complementParams, npRead->complementEvents,
                                                                npRead->nbComplementEvents, npRead->complementEventMap,
                                                                rc_targetSeq, p, anchorPairs);
        if (posteriorProbsFile != NULL) {
            fprintf(stderr, "vanillaAlign - writing %lld template aligned pairs", stList_length(templateAlignedPairs));
            fprintf(stderr, "vanillaAlign - writing %lld complement aligned pairs\n", stList_length(complementAlignedPairs));
            writePosteriorProbs(posteriorProbsFile, templateAlignedPairs, template);
            writePosteriorProbs(posteriorProbsFile, complementAlignedPairs, complement);
        }
    }

    return 0;
}
