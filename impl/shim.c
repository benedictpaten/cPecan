// AR

#include "shim.h"
#include "../inc/shim.h"
#include "../inc/pairwiseAligner.h"
#include "../inc/emissionMatrix.h"
#include <assert.h>


// Sequence constructor function
Sequence* sequenceConstruct(int stringLength, void *elements, sequenceType type) {

    Sequence* self = malloc(sizeof(Sequence));
    self->length = correctSeqLength(stringLength, type);
    self->elements = elements;
    //self->get = getfPtr;
    switch (type) {
        case 0:
            self->get = getBase;
            break;
        case 1:
            self->get = getKmer;
            break;
        case 2:
            self->get = getEvent;
            break;
    }
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

// TODO need a proper unit test function for this
int64_t getKmerIndex(char* kmer) {
    //int64_t kmerLength = strlen(kmer);
    //printf("kmer length: %lld\n", kmerLength);
    int64_t kmerLen = strlen(kmer);
    assert(kmerLen == KMER_LENGTH);
    int64_t axisLength = 25; // for 2-mers
    int64_t l = axisLength/5;
    int64_t i = 0;
    int64_t x = 0; //instead of index
    while(l > 1) {
        //printf("at start x:%lld, l:%lld\n", x, l);
        x += l*getBaseIndex(kmer[i]);
        //printf("after math, x:%lld, gBI:%lld\n", x, getBaseIndex(kmer[i]));
        i += 1;
        l = l/5;
    }

    int64_t last = strlen(kmer)-1;
    //printf("last:%lld\n", last);
    x += getBaseIndex(kmer[last]);


    return x;
}

int64_t correctSeqLength(sequenceType type, int64_t stringLength) {
    if (stringLength == 0) {
        return 0;
    }
    if (stringLength > 0) {
        switch (type) {
            case 0: // nucleotide sequence
                return stringLength;
            case 1: // event and kmer sequence
            case 2:
                return stringLength - 1;
        }
    }
}

int64_t getKmerSeqLength(int64_t stringLength) {
    //int64_t l = stringLength;
    if (stringLength == 0) {
        return 0;
    }
    if (stringLength > 0) {
        return stringLength - 1;
    }
}


// returns a pointer to a kmer within a char array
void* getKmer(void *elements, int64_t index) {
    char* n;
    n = "NN"; // hardwired null kmer
    int64_t i = index;
    // change kmer length here, hardwired so far...
    int kmerLength = KMER_LENGTH;
    char *k_i = malloc((kmerLength+1) * sizeof(char)); // TODO remove this malloc
    for (int x = 0; x < kmerLength; x++) {
        k_i[x] = *((char *)elements+(i+x));
    }

    return index >= 0 ? k_i : n;

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
