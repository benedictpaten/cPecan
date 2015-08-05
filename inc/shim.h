// another attempt at making a generalized sequence object
// Art Rand

#ifndef SHIM_H
#define SHIM_H
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>



// Generalized sequence struct
typedef struct sequence {
    int64_t length;
    void *elements;
    void* (*get)(void *elements, int64_t index);
    char n;
} Sequence;

// Types of sequences
typedef enum {
    nucleotide=0,
    kmer=1,
    event=2
} sequenceType;

int64_t correctSeqLength(sequenceType type, int64_t stringLength);

Sequence* sequenceConstruct(int length, void *elements, sequenceType type);

void sequenceDestroy(Sequence* seq);

typedef char base;

//typedef enum {
//    match=0,
//    transversion=1,
//    transition=2,
//} baseMatchTypes;

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
