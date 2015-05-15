// AR

#include "shim.h"
#include "../inc/shim.h"
#include "../inc/pairwiseAligner.h"
#include <assert.h>


// Sequence constructor function
Sequence* sequenceConstruct(int length, void *elements, void (*getfPtr)) {
    Sequence* self = malloc(sizeof(Sequence));
    
    self->length = length;
    self->elements = elements;
    self->get = getfPtr;
    //self->n = "n"; // TODO decide on a generic empty/Null character?
    return self;
}

// sequence destroying function
void sequenceDestroy(Sequence* seq) {
    assert(seq != NULL);
    free(seq);
}

// Get functions: retrieve a single element (base, kmer, or event) from an 
// array of elements.
 
// returns a pointer to base in a char array
void* getBase(void *elements, int64_t index) {
    char* n;
    n = "n";
    return index >= 0 ? &(((char *)elements)[index]) : n;
}

// hack that will exchange a char representing a base for an integer to make them
// compatable with other functions that usually take enum datatype
int64_t getBaseIndex(char base) {
    //printf("Running getBaseIndex\n");
    //char b;
    //b = *base;
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
    }
}


// returns a pointer to a kmer within a char array
void* getKmer(void *elements, int64_t index) {
    int64_t i = index;
    int kmerLength = 5;
    
    char *k_i = malloc((kmerLength+1) * sizeof(char));
    
    for (int x = 0; x < kmerLength; x++) {
        k_i[x] = *((char *)elements+(i+x));
    }
    
    return k_i;
}

// TODO make a function that compares kmers

// getEvent function, returns a pointer to a event in a sequence
void* getEvent(void* elements, int64_t index) {
    return (Event *) elements + index;
}


// Event-specific functions:
// event constructor function.
Event* event_construct(double mean, char kmer[6]) {
    Event *event = malloc(sizeof(Event));
    event->mean = mean;
    event->kmer = kmer;
    event->base = kmer[0];
    
    return event;
}

// eventSequenceConstruct function, takes arrays of means and kmers and
// returns an array of pointers to event objects
void* eventSequenceConstruct(int length, double* means, char* kmers) {
    Event* eventArray = malloc(length * sizeof(Event));
    for (int i = 0; i < length; i++) {
        eventArray[i] = *(event_construct(means[i], kmers+(i*6)));
    }
    
    return eventArray;
}

// eventSequence destroying function
void eventSequenceDestroy(void* eventArray) {
    assert(eventArray != NULL);
    free(eventArray);
}
