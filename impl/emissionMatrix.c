#include "emissionMatrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

void emissionsKmers_setMatchProbsToDefaults(double *emissionMatchProbs) {
    /*
     * This sets the match probabilities to default values for matching kmers
     */

    const double M=-2.1149196655034745; //log(0.12064298095701059);
    const double V=-4.5691014376830479; //log(0.010367271172731285);
    const double S=-3.9833860032220842; //log(0.01862247669752685);
    const double N=-2.772588722;        //log(0.25**2)

    //Symmetric matrix of emission probabilities.

    const double i[MATRIX_SIZE] = {
        M, V, S, V, N, 
        M, V*V, V*S, V*V, V*N, 
        M, S*V, S*S, S*V, S*N, 
        M, V*V, V*S, V*V, V*N, 
        M, N*V, N*S, N*V, N*N, 
        V, M, V, S, N, 
        V*V, M, V*V, V*S, V*N, 
        S*V, M, S*V, S*S, S*N, 
        V*V, M, V*V, V*S, V*N, 
        N*V, M, N*V, N*S, N*N, 
        S, V, M, V, N, 
        V*S, V*V, M, V*V, V*N, 
        S*S, S*V, M, S*V, S*N, 
        V*S, V*V, M, V*V, V*N, 
        N*S, N*V, M, N*V, N*N, 
        V, S, V, M, N, 
        V*V, V*S, V*V, M, V*N, 
        S*V, S*S, S*V, M, S*N, 
        V*V, V*S, V*V, M, V*N, 
        N*V, N*S, N*V, M, N*N, 
        N, N, N, N, N, 
        V*N, V*N, V*N, V*N, N, 
        S*N, S*N, S*N, S*N, N, 
        V*N, V*N, V*N, V*N, N, 
        N*N, N*N, N*N, N*N, N, 
        M, V*V, V*S, V*V, V*N, 
        M, V, S, V, N, 
        M, V*V, V*S, V*V, V*N, 
        M, S*V, S*S, S*V, S*N, 
        M, N*V, N*S, N*V, N*N, 
        V*V, M, V*V, V*S, V*N, 
        V, M, V, S, N, 
        V*V, M, V*V, V*S, V*N, 
        S*V, M, S*V, S*S, S*N, 
        N*V, M, N*V, N*S, N*N, 
        V*S, V*V, M, V*V, V*N, 
        S, V, M, V, N, 
        V*S, V*V, M, V*V, V*N, 
        S*S, S*V, M, S*V, S*N, 
        N*S, N*V, M, N*V, N*N, 
        V*V, V*S, V*V, M, V*N, 
        V, S, V, M, N, 
        V*V, V*S, V*V, M, V*N, 
        S*V, S*S, S*V, M, S*N, 
        N*V, N*S, N*V, M, N*N, 
        V*N, V*N, V*N, V*N, N, 
        N, N, N, N, N, 
        V*N, V*N, V*N, V*N, N, 
        S*N, S*N, S*N, S*N, N, 
        N*N, N*N, N*N, N*N, N, 
        M, S*V, S*S, S*V, S*N, 
        M, V*V, V*S, V*V, V*N, 
        M, V, S, V, N, 
        M, V*V, V*S, V*V, V*N, 
        M, N*V, N*S, N*V, N*N, 
        S*V, M, S*V, S*S, S*N, 
        V*V, M, V*V, V*S, V*N, 
        V, M, V, S, N, 
        V*V, M, V*V, V*S, V*N, 
        N*V, M, N*V, N*S, N*N, 
        S*S, S*V, M, S*V, S*N, 
        V*S, V*V, M, V*V, V*N, 
        S, V, M, V, N, 
        V*S, V*V, M, V*V, V*N, 
        N*S, N*V, M, N*V, N*N, 
        S*V, S*S, S*V, M, S*N, 
        V*V, V*S, V*V, M, V*N, 
        V, S, V, M, N, 
        V*V, V*S, V*V, M, V*N, 
        N*V, N*S, N*V, M, N*N, 
        S*N, S*N, S*N, S*N, N, 
        V*N, V*N, V*N, V*N, N, 
        N, N, N, N, N, 
        V*N, V*N, V*N, V*N, N, 
        N*N, N*N, N*N, N*N, N, 
        M, V*V, V*S, V*V, V*N, 
        M, S*V, S*S, S*V, S*N, 
        M, V*V, V*S, V*V, V*N, 
        M, V, S, V, N, 
        M, N*V, N*S, N*V, N*N, 
        V*V, M, V*V, V*S, V*N, 
        S*V, M, S*V, S*S, S*N, 
        V*V, M, V*V, V*S, V*N, 
        V, M, V, S, N, 
        N*V, M, N*V, N*S, N*N, 
        V*S, V*V, M, V*V, V*N, 
        S*S, S*V, M, S*V, S*N, 
        V*S, V*V, M, V*V, V*N, 
        S, V, M, V, N, 
        N*S, N*V, M, N*V, N*N, 
        V*V, V*S, V*V, M, V*N, 
        S*V, S*S, S*V, M, S*N, 
        V*V, V*S, V*V, M, V*N, 
        V, S, V, M, N, 
        N*V, N*S, N*V, M, N*N, 
        V*N, V*N, V*N, V*N, N, 
        S*N, S*N, S*N, S*N, N, 
        V*N, V*N, V*N, V*N, N, 
        N, N, N, N, N, 
        N*N, N*N, N*N, N*N, N, 
        M, N*V, N*S, N*V, N*N, 
        M, N*V, N*S, N*V, N*N, 
        M, N*V, N*S, N*V, N*N, 
        M, N*V, N*S, N*V, N*N, 
        M, V, S, V, N, 
        N*V, M, N*V, N*S, N*N, 
        N*V, M, N*V, N*S, N*N, 
        N*V, M, N*V, N*S, N*N, 
        N*V, M, N*V, N*S, N*N, 
        V, M, V, S, N, 
        N*S, N*V, M, N*V, N*N, 
        N*S, N*V, M, N*V, N*N, 
        N*S, N*V, M, N*V, N*N, 
        N*S, N*V, M, N*V, N*N, 
        S, V, M, V, N, 
        N*V, N*S, N*V, M, N*N, 
        N*V, N*S, N*V, M, N*N, 
        N*V, N*S, N*V, M, N*N, 
        N*V, N*S, N*V, M, N*N, 
        V, S, V, M, N, 
        N*N, N*N, N*N, N*N, N, 
        N*N, N*N, N*N, N*N, N, 
        N*N, N*N, N*N, N*N, N, 
        N*N, N*N, N*N, N*N, N, 
        N, N, N, N, N, 
        };

    memcpy(emissionMatchProbs, i, sizeof(double)*MATRIX_SIZE);
}


void emissionsKmers_setGapProbsToDefaults(double *emissionGapProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double G = -1.6094379124341003; //log(0.2)
    const double i[NUM_OF_KMERS] = {
        G, G, G, G, G, 
        G, G, G, G, G, 
        G, G, G, G, G, 
        G, G, G, G, G, 
        G, G, G, G, G, 
        };

    memcpy(emissionGapProbs, i, sizeof(double)*25);
}


