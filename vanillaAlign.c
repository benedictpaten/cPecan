#include <getopt.h>
#include "pairwiseAlignment.h"
#include "pairwiseAligner.h"
#include "emissionMatrix.h"
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

void writePosteriorProbs(char *posteriorProbsFile, char *readFile, double *matchModel, double *events, char *target,
                         stList *alignedPairs, Strand strand) {
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
    // x, k_x, name, strand, e_u, e_o, e_t, p, level_E_u, noise_E_u
    // c   c    c     c       c    c    c    c    c         c  CHECK!
    FILE *fH = fopen(posteriorProbsFile, "a");
    for(int64_t i = 0; i < stList_length(alignedPairs); i++) {
        // grab the aligned pair
        stIntTuple *aPair = stList_get(alignedPairs, i);

        // unpack it to make this easier
        int64_t x = stIntTuple_get(aPair, 1);         // target index
        int64_t y = stIntTuple_get(aPair, 2);         // event index
        double p = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1;  // posterior prob

        // make the kmer at the target index
        char *k_i = st_malloc(KMER_LENGTH * sizeof(char));
        for (int64_t k = 0; k < KMER_LENGTH; k++) {
            k_i[k] = *(target + (x + k));
        }
        // get the kmer at the target index
        int64_t targetKmerIndex = emissions_discrete_getKmerIndexFromKmer(k_i);

        // get the observations from the events
        //double eventMean = *(events + y);
        double eventMean = events[(y * NB_EVENT_PARAMS)];
        //double eventNoise = *(events + (y + 1));
        double eventNoise = events[(y * NB_EVENT_PARAMS) + 1];
        //double eventDuration = *(events + (y + 2));
        double eventDuration = events[(y * NB_EVENT_PARAMS) + 2];

        // get the expected event mean amplitude and noise
        double E_levelu = matchModel[1 + (targetKmerIndex * MODEL_PARAMS)];
        double E_noiseu = matchModel[1 + (targetKmerIndex * MODEL_PARAMS + 2)];

        fprintf(fH, "%lld\t%lld\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, k_i, readFile, strandLabel, eventMean, eventNoise,
                eventDuration, p, E_levelu, E_noiseu);
        //fprintf(fH, "%s\t" "%" PRIi64 "\t%" PRIi64 "\t%f\n", strandLabel, stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2), ((double)stIntTuple_get(aPair, 0))/PAIR_ALIGNMENT_PROB_1);
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

stList *performSignalAlignment(StateMachine *sM, const char *hmmFile, double *events, int64_t nbEvents,
                               int64_t *eventMap, char *target, PairwiseAlignmentParameters *p,
                               stList *unmappedAncors) {
    if ((sM->type != threeState) && (sM->type != vanilla)) {
        st_errAbort("you're trying to do the wrong kind of alignment\n");
    }

    // load model into stateMachine
    //StateMachine *sM = buildStateMachine(model, npp, type);

    // load HMM if given
    if (hmmFile != NULL) {
        fprintf(stderr, "loading HMM from file, %s\n", hmmFile);
        loadHmmRoutine(hmmFile, sM, sM->type);
    }

    // do alignment
    st_uglyf("SENTINAL - about to start alignmentp\n");
    if (sM->type == vanilla) {
        stList *alignedPairs = performSignalAlignmentP(sM, events, nbEvents, eventMap, target, p, unmappedAncors,
                                                       sequence_getKmer2);
        return alignedPairs;
    }
    if (sM->type == threeState) {
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
    fprintf(stderr, "\nvanillaAlign - finished iterations\n");

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
    char *readLabel = NULL;
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
                {"readLabel",               required_argument,  0,  'L'},
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

        key = getopt_long(argc, argv, "h:s:T:C:L:q:r:u:y:z:t:c:i:", long_options, &option_index);

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
            fprintf(stderr, "vanillaAlign - writing hmm to file: %s\n\n", templateTrainedHmmFile);
            hmmContinuous_writeToFile(templateTrainedHmmFile, templateTrainedHmm, sMtype);
        }
        if (complementTrainedHmmFile != NULL) {
            fprintf(stderr, "vanillaAlign - starting training for complement hmm\n");
            Hmm *complementTrainedHmm = performBaumWelchTraining(complementModelFile, complementHmmFile, sMtype,
                                                                 npRead->complementParams, npRead->complementEvents,
                                                                 npRead->nbComplementEvents, npRead->complementEventMap,
                                                                 rc_targetSeq, p, anchorPairs, iter);
            fprintf(stderr, "vanillaAlign - writing hmm to file: %s\n\n", complementTrainedHmmFile);
            hmmContinuous_writeToFile(complementTrainedHmmFile, complementTrainedHmm, sMtype);
        }
    } else {
        fprintf(stderr, "vanillaAlign - starting alignment\n");
        // Template alignment
        // load template stateMachine
        StateMachine *sMt = buildStateMachine(templateModelFile, npRead->templateParams, sMtype);
        // get aligned pairs
        stList *templateAlignedPairs = performSignalAlignment(sMt, templateHmmFile, npRead->templateEvents,
                                                              npRead->nbTemplateEvents, npRead->templateEventMap,
                                                              targetSeq, p, anchorPairs);
        // sort
        stList_sort(templateAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing
        // write to file
        if (posteriorProbsFile != NULL) {
            fprintf(stderr, "vanillaAlign - writing %lld template aligned pairs\n",
                    stList_length(templateAlignedPairs));
            writePosteriorProbs(posteriorProbsFile, readLabel, sMt->EMISSION_MATCH_PROBS, npRead->templateEvents,
                                targetSeq, templateAlignedPairs, template);
        }
        // Complement alignment
        // load complement stateMachine
        StateMachine *sMc = buildStateMachine(complementModelFile, npRead->complementParams, sMtype);
        // get aligned pairs
        stList *complementAlignedPairs = performSignalAlignment(sMc, complementHmmFile, npRead->complementEvents,
                                                                npRead->nbComplementEvents, npRead->complementEventMap,
                                                                rc_targetSeq, p, anchorPairs);
        // sort
        stList_sort(complementAlignedPairs, sortByXPlusYCoordinate2); //Ensure the coordinates are increasing
        // write to file
        if (posteriorProbsFile != NULL) {
            fprintf(stderr, "vanillaAlign - writing %lld complement aligned pairs\n",
                    stList_length(complementAlignedPairs));
            writePosteriorProbs(posteriorProbsFile, readLabel, sMt->EMISSION_MATCH_PROBS, npRead->complementEvents,
                                rc_targetSeq, complementAlignedPairs, complement);
        }
    }

    return 0;
}
