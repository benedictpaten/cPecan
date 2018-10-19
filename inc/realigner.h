/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef REALIGNER_H_
#define REALIGNER_H_

#include "sonLib.h"
#include "pairwiseAligner.h"

/*
 * Basic data structures for representing a POA alignment.
 */
 
typedef struct _Poa {
	char *refString; // The reference string
	stList *nodes; 
} Poa;

typedef struct _poaNode {
	stList *inserts; // Inserts that happen immediately after this position
	stList *deletes; // Deletes that happen immediately after this position
	char base; // Char representing base, e.g. 'A', 'C', etc.
	double *baseWeights; // Array of length SYMBOL_NUMBER, encoding the weight given go each base, using the Symbol enum
} PoaNode;

typedef struct _poaInsert {
	char *insert; // String representing characters of insert e.g. "GAT", etc.
	double weight;
} PoaInsert;

typedef struct _poaDelete {
	int64_t length; // Length of delete
	double weight;
} PoaDelete;

/*
 * Creates a POA representing the given reference sequence, with one node for each reference base and a 
 * prefix 'N' base to represent place to add inserts/deletes that precede the first position of the reference.
 */
Poa *poa_getReferenceGraph(char *reference);

/*
 * Adds to given POA the matches, inserts and deletes from the alignment of the given read to the reference.
 * Adds the inserts and deletes so that they are left aligned.
 */
void poa_augment(Poa *poa, char *read, stList *matches, stList *inserts, stList *deletes);

/*
 * Creates a POA representing the reference and the expected inserts / deletes and substitutions from the 
 * alignment of the given set of reads aligned to the reference.
 */
Poa *poa_realign(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p);

/*
 * Prints representation of the POA.
 */
void poa_print(Poa *poa, FILE *fH, float indelSignificanceThreshold);

/*
 * Prints some summary stats on the POA.
 */
void poa_printSummaryStats(Poa *poa, FILE *fH);

/*
 * Left align indels and then sorts them by character
 */
void poa_normalize(Poa *poa);

/*
 * Creates a consensus reference sequence from the POA
 */
char *poa_getConsensus(Poa *poa);

/*
 * Iteratively used poa_realign and poa_getConsensus to refine the median reference sequence for the given reads and 
 * the starting reference.
 */
Poa *poa_realignIterative(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p,
				 uint64_t iterations);
				 
void poa_destruct(Poa *poa);

/*
 * Finds shift, expressed as a reference coordinate, that the given substring str can
 * be shifted left in the refString, starting from a match at refStart.
 */
int64_t getShift(char *refString, int64_t refStart, char *str, int64_t length);

/*
 * Get sum of weights for reference bases in poa - proxy to agreement of reads
 * with reference.
 */
double poa_getReferenceNodeTotalMatchWeight(Poa *poa);

/*
 * Get sum of weights for non-reference bases in poa - proxy to disagreement of read positions
 * aligned with reference.
 */
double poa_getReferenceNodeTotalDisagreementWeight(Poa *poa);

/*
 * Get sum of weights for delete in poa - proxy to delete disagreement of reads
 * with reference.
 */
double poa_getDeleteTotalWeight(Poa *poa);

/*
 * Get sum of weights for inserts in poa - proxy to insert disagreement of reads
 * with reference.
 */
double poa_getInsertTotalWeight(Poa *poa);

#endif /* REALIGNER_H_ */
