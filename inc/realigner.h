/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef REALIGNER_H_
#define REALIGNER_H_

#include "sonLib.h"
#include "pairwiseAligner.h"

typedef struct _poaNode {
	stList *inserts;
	stList *deletes;
	char base;
	double *baseWeights;
} PoaNode;

typedef struct _poaInsert {
	char *insert;
	double weight;
} PoaInsert;

typedef struct _poaDelete {
	int64_t length;
	double weight;
} PoaDelete;

typedef struct _Poa {
	stList *nodes;
} Poa;

void poa_augment(Poa *poa, char *read, stList *matches, stList *inserts, stList *deletes);

Poa *poa_realign(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p);

void poa_print(Poa *poa, FILE *fH);

char *poa_getConsensus(Poa *poa);

Poa *poa_realignIterative(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p,
				 uint64_t iterations);

#endif /* REALIGNER_H_ */
