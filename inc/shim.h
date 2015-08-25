// another attempt at making a generalized sequence object
// Art Rand

#ifndef SHIM_H
#define SHIM_H
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// Types of sequences
typedef enum {
    nucleotide=0,
    kmer=1,
    event=2
} sequenceType;

// Generalized sequence struct
typedef struct sequence {
    int64_t length;
    void *elements;
    char* repr;
    sequenceType type;
    void* (*get)(void *elements, int64_t index);
} Sequence;


int64_t correctSeqLength(int64_t stringLength, sequenceType type);

Sequence* sequenceConstruct(int64_t length, void *elements, sequenceType type);

Sequence* sequence_getSubSequence(Sequence* wholeSequence, int64_t start, int64_t length, sequenceType t);

void sequenceDestroy(Sequence* seq);

typedef char base;

void* getBase(void *elements, int64_t index);

int64_t getBaseIndex(char base);

int64_t getKmerIndex(char* kmer);

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
