#include "pairwiseAligner.h"
#include "stateMachine.h"
#include "nanopore.h"

void writePosteriorProbs(char *posteriorProbsFile, stList *alignedPairs) {
    /*
     * Writes the posterior match probabibilities to a tab separated file, each line being X coordinate, Y coordinate, Match probability
     */
    FILE *fH = fopen(posteriorProbsFile, "w");
    for(int64_t i=0;i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        fprintf(fH, "%" PRIi64 "\t%" PRIi64 "\t%f\n", stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2), ((double)stIntTuple_get(aPair, 0))/PAIR_ALIGNMENT_PROB_1);
    }
    fclose(fH);
}

int main(int argc, char *argv[]) {
    printf("Hello from signalAlign\n");

    return 0;
}
