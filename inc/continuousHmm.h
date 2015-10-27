#ifndef CONTINUOUS_HMM_H
#define CONTINUOUS_HMM_H

typedef struct _hmmContinuous {
    Hmm baseHmm;
    // these are just for storage not EM, yet
    double *matchModel;
    double *extraEventMatchModel;
} HmmContinuous;

typedef struct _strawManHmm {
    HmmContinuous baseContinuousHmm;
    double *transitions;
    double *individualKmerSkipProbs; // for learning skip probs/kmer
} ContinuousPairHmm;

typedef struct _vanillaHmm {
    HmmContinuous baseContinuousHmm;
    double *kmerSkipBins;
    int64_t (*getKmerSkipBin)(double *matchModel, void *cX);
} VanillaHmm;

Hmm *continuousPairHmm_constructEmpty(double pseudocount, int64_t stateNumber,
                                      int64_t symbolSetSize, StateMachineType type,
                                      void (*addToTransitionExpFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                      void (*setTransitionFcn)(Hmm *hmm, int64_t from, int64_t to, double p),
                                      double (*getTransitionsExpFcn)(Hmm *hmm, int64_t from, int64_t to),
                                      int64_t (*getElementIndexFcn)(void *));


#endif