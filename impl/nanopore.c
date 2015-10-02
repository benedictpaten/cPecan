#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "pairwiseAligner.h"
/*
#include "sonLibFile.h"
#include "sonLibTypes.h"
#include "sonLibString.h"
#include "sonLibList.h"
#include "sonLibCommon.h"
 */
#include "nanopore.h"

#define NB_EVENT_PARAMS 3

NanoporeRead *nanoporeRead_construct(int64_t length, int64_t nb_events, double scale, double shift) {
    // do a malloc
    NanoporeRead *npRead = st_malloc(sizeof(NanoporeRead));

    // set parts to size
    // initialize parts
    npRead->length = length;
    npRead->nb_events = nb_events;
    npRead->scale = scale,
    npRead->shift = shift;
    npRead->read = st_malloc(npRead->length * sizeof(char));
    npRead->eventMap = st_malloc(npRead->length * sizeof(int64_t));
    npRead->events = st_malloc(npRead->nb_events * NB_EVENT_PARAMS * sizeof(double));
    // return
    return npRead;
}

NanoporeRead *loadNanoporeReadFromFile(const char *nanoporeReadFile) {
    FILE *fH = fopen(nanoporeReadFile, "r");
    // the first line is the length of the read (in bases), the scale, and the shift
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    int64_t npRead_length, npRead_nbEvents;
    double npRead_scale, npRead_shift;

    int64_t j = sscanf(stList_get(tokens, 0), "%lld", &npRead_length);
    if (j != 1) {
        st_errAbort("error parsing nanopore read length\n");
    }
    int64_t h = sscanf(stList_get(tokens, 1), "%lld", &npRead_nbEvents);
    if (j != 1) {
        st_errAbort("error parsing nanopore read length\n");
    }
    int64_t s = sscanf(stList_get(tokens, 2), "%lf", &npRead_scale);
    if (s != 1) {
        st_errAbort("error parsing nanopore read scale\n");
    }
    int64_t n = sscanf(stList_get(tokens, 3), "%lf", &npRead_shift);
    if (s != 1) {
        st_errAbort("error parsing nanopore read shift\n");
    }

    NanoporeRead *npRead = nanoporeRead_construct(npRead_length, npRead_nbEvents,
                                                  npRead_scale, npRead_shift);

    // cleanup first line
    free(string);
    stList_destruct(tokens);

    // read sequence
    string = stFile_getLineFromFile(fH);
    j = sscanf(string, "%s", npRead->read);
    if (j != 1) {
        st_errAbort("error parsing read from npRead file\n");
    }
    free(string);

    // event map
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->length) {
        st_errAbort(
                "event map is not the correct length, should be %lld, got %lld",
                npRead->length,
                stList_length(tokens));
    }
    // load in the map
    for (int64_t i = 0; i < npRead->length; i++) {
        j = sscanf(stList_get(tokens, i), "%lld", &(npRead->eventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // events
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    // load in events
    for (int64_t i = 0; i < (npRead->nb_events * NB_EVENT_PARAMS); i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->events[i]));
        if (j != 1) {
            st_errAbort("error loading in events\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    fclose(fH);

    return npRead;
}