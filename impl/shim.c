// AR

#include <inttypes.h>
#include "shim.h"
#include "../inc/shim.h"
#include "../inc/pairwiseAligner.h"
#include "../inc/emissionMatrix.h"
#include <assert.h>


/*
 * Sequence constructor function
 * stringLength should be the length of the sequence in bases ie ATGAC has
 * length 5.  *elements is a pointer to the char array, typically.  The
 * sequenceType is 0-nucleotide, 1-kmer, 2-event
 */
Sequence* sequenceConstruct(int64_t stringLength, void *elements, sequenceType t) {
    Sequence* self = malloc(sizeof(Sequence));
    // correct the sequence length for kmers/events
    self->length = correctSeqLength(stringLength, t);
    self->type = t;
    self->elements = elements;
    self->repr = (char*) elements;
    switch (t) {
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
    return self;
}

/*
 * Function to retrieve a sub sequence of a sequence object.
 */
Sequence* sequence_getSubSequence(Sequence* wholeSequence, int64_t start, int64_t length, sequenceType t) {
    char* wS_string = wholeSequence->repr;
    char* subString = stString_getSubString(wS_string, start, length);
    Sequence* subSequence = sequenceConstruct(length, subString, t);
    return subSequence;
}

void sequenceDestroy(Sequence* seq) {
    //assert(seq != NULL);
    free(seq);
}

/*
 * Returns a single base from a sequence object
 */
void* getBase(void *elements, int64_t index) {
    char* n;
    n = "n";
    return index >= 0 ? &(((char *)elements)[index]) : n;
}

/*
 * Returns a kmer from a sequence object
 */
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

/*
 * Returns the index for a base, for use with matrices and getKmerIndex
 */
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

/*
 * Returns the index of a kmer
 */
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

/*
 * Correct the sequence length for non-nucleotide sequences, eg. kmers/events.
 */
int64_t correctSeqLength(int64_t stringLength, sequenceType type) {
    // for trivial case
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

/*
 * Returns an event from an eventSequence
 */
void* getEvent(void* elements, int64_t index) {
    return (Event *) elements + index;
}

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
