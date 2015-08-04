/*
 * randomSequences.c
 *
 *  Created on: 3 Mar 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sonLib.h"
#include "../../sonLib/lib/sonLibRandom.h"
#include "../../sonLib/lib/sonLibCommon.h"
#include "../../sonLib/lib/sonLibString.h"

char getRandomChar() {
    //char *positions = "AaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtN";
    char *positions = "AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTN";
    return positions[st_randomInt(0, strlen(positions))];
}

/*
 * Creates a random DNA sequence of the given length.
 */
char *getRandomSequence(int64_t length) {
    printf("1s. Running getRandomSequence\n");
    char *seq = st_malloc((length + 1) * sizeof(char));
    for (int64_t i = 0; i < length; i++) {
        seq[i] = getRandomChar();
    }
    seq[length] = '\0';
    return seq;
}

/*
 * Transfroms the given sequence into a different sequence.
 */
char *evolveSequence(const char *startSequence) {
    printf("2s. Running evolveSequence\n");
    //Copy sequence
    char *seq = stString_copy(startSequence);

    //Do substitutions
    for (int64_t i = 0; i < strlen(seq); i++) {
        if (st_random() > 0.8) {
            seq[i] = getRandomChar();
        }
    }

    //Do indels
    while (st_random() > 0.2) {
        char *toReplace = getRandomSequence(st_randomInt(2, 4));
        char *replacement = getRandomSequence(st_randomInt(0, 10));
        char *seq2 = stString_replace(seq, toReplace, replacement);
        free(seq);
        free(toReplace);
        free(replacement);
        seq = seq2;
    }

    return seq;
}
