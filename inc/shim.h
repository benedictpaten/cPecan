// another attempt at making a generalized sequence object
// Art Rand

#ifndef SHIM_H
#define SHIM_H
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// Generalized sequence struct
typedef struct sequence {
    int64_t length;
    void *elements;
    void *(*get)(void *elements, int64_t index);
} Sequence;

Sequence* sequenceConstruct(int length, void *elements, void (*getFcn));

void sequenceDestroy(Sequence* seq);

typedef char base;

void* getBase(void *elements, int64_t index);

void* getKmer(void *elements, int64_t index);


typedef struct event {
    double mean;
    char *kmer;
    char base;
} Event;

void* getEvent(void* elements, int64_t index);

void* eventSequenceConstruct(int length, double* means, char* kmers);

void eventSequenceDestroy(void* eventArray);
#endif
