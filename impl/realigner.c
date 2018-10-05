/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "pairwiseAligner.h"
#include "realigner.h"
#include <stdlib.h>
#include <math.h>
#include "stGraph.h"
#include <inttypes.h>

/*
 * Basic overview of realigner algorithm:
 *
 * Inputs: a reference sequence and set of reads, each with a pairwise alignment to the reference
 * Outputs: a modified reference sequence, for each read a pairwise alignment to the reference
 *
 * Steps:
 * For each read, using its guide alignment, generate posterior match and gap probabilities
 * Construct the posterior graph
 * Pick new reference as path through posterior graph
 */

PoaInsert *poaInsert_construct(char *insert, double weight) {
	PoaInsert *poaInsert = st_calloc(1, sizeof(PoaInsert));

	poaInsert->insert = insert;
	poaInsert->weight = weight;

	return poaInsert;
}

void poaInsert_destruct(PoaInsert *poaInsert) {
	free(poaInsert->insert);
	free(poaInsert);
}

PoaDelete *poaDelete_construct(int64_t length, double weight) {
	PoaDelete *poaDelete = st_calloc(1, sizeof(PoaDelete));

	poaDelete->length = length;
	poaDelete->weight = weight;

	return poaDelete;
}

void poaDelete_destruct(PoaDelete *poaDelete) {
	free(poaDelete);
}

PoaNode *poaNode_construct(char base) {
	PoaNode *poaNode = st_calloc(1, sizeof(poaNode));

	poaNode->inserts = stList_construct3(0, (void(*)(void *)) poaInsert_destruct);
	poaNode->deletes = stList_construct3(0, (void(*)(void *)) poaDelete_destruct);
	poaNode->base = base; // A = 0, C = 1, G = 2, T = 3, N = 4
	poaNode->baseWeights = st_calloc(SYMBOL_NUMBER, sizeof(double)); // A = 0, C = 1, G = 2, T = 3, N = 4

	return poaNode;
}

void poaNode_destruct(PoaNode *poaNode) {
	stList_destruct(poaNode->inserts);
	stList_destruct(poaNode->deletes);
	free(poaNode->baseWeights);
	free(poaNode);
}

Poa *poa_getEmptyGraph() {
	Poa *poa = st_calloc(1, sizeof(poa));

	poa->nodes = stList_construct3(0, (void (*)(void *))poaNode_destruct);

	return poa;
}

void poa_destruct(Poa *poa) {
	stList_destruct(poa->nodes);
	free(poa);
}

Poa *poa_getReferenceGraph(char *reference) {
	/*
	 *
	 */

	Poa *poa = poa_getEmptyGraph();

	int64_t refLength = strlen(reference);
	stList_append(poa->nodes, poaNode_construct(0)); // Add empty prefix node
	for(int64_t i=0; i<refLength; i++) {
		stList_append(poa->nodes, poaNode_construct(reference[i]));
	}

	return poa;
}

int cmpMatches(const void *a, const void *b) {
	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 0) < stIntTuple_get(two, 0) ? -1 : stIntTuple_get(one, 0) > stIntTuple_get(two, 0) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	}

	return i;
}

int cmpDeletes(const void *a, const void *b) {
	// As cmpMatches, but compares by read coordinate and then reference coordinate

	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 0) < stIntTuple_get(two, 0) ? -1 : stIntTuple_get(one, 0) > stIntTuple_get(two, 0) ? 1 : 0;
	}

	return i;
}

bool isMatch(stSortedSet *matchesSet, int64_t x, int64_t y) {
	return 1;
}

void poa_augment(Poa *poa, char *read, stList *matches, stList *inserts, stList *deletes) {
	/*
	 *
	 */

	// Add weights of matches to the POA graph

	// For each match in alignment subgraph identify its corresponding node in the POA graph
	// add the weight of the match to the POA node
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		PoaNode *node = stList_get(poa->nodes, stIntTuple_get(match, 0)+1); // Get POA node

		char base = read[stIntTuple_get(match, 1)]; // Get base

		node->baseWeights[symbol_convertCharToSymbol(base)] += stIntTuple_get(match, 2); // Add base weight to POA node
	}

	// Create a set of match coordinates

	stSortedSet *matchesSet = stSortedSet_construct3(cmpMatches, NULL);
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		assert(stSortedSet_search(matchesSet, match) == NULL);
		stSortedSet_insert(matchesSet, match); // Add to matches
	}
	
	// Add inserts to the POA graph

	// Sort the inserts first by the reference coordinate and then by read coordinate
	stList_sort(inserts, cmpMatches);

	// Let a complete-insert be a sequence of inserts with the same reference coordinate i
	// and consecutive read coordinates j, j+1, ..., j+n, such that the i,j-1 is a match
	// in the alignment subgraph and i+1,j+n is similarly a match.

	// Enumerate set of complete inserts, adding them to the graph
	for(int64_t i=0; i<stList_length(inserts); i++) {
		stIntTuple *insertStart = stList_get(inserts, i); // Start of putative complete-insert

		int64_t j=i+1;
		for(;j<stList_length(inserts); j++) {
			stIntTuple *insertEnd = stList_get(inserts, j); // End of putative complete-insert

			// If they don't have the same reference coordinate then not part of same complete-insert
			if(stIntTuple_get(insertStart, 0) != stIntTuple_get(insertEnd, 0)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(insertStart, 1) + j - i != stIntTuple_get(insertEnd, 1)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a putative maximal
		// complete insert sequence in inserts

		// Now enumerate complete-inserts in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match then can not be a complete insert
			if(!isMatch(matchesSet, stIntTuple_get(insertStart, 0), stIntTuple_get(insertStart, 1) + k - i - 1)) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceding match then can not be a complete insert
				if(!isMatch(matchesSet, stIntTuple_get(insertStart, 0) + 1, stIntTuple_get(insertStart, 1) + l - i + 1)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-insert

				// Calculate weight and label
				double insertWeight = 1.0;
				char insertLabel[l+2-k];
				insertLabel[l+1-k] = '\0';
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *insert = stList_get(inserts, m);
					insertLabel[m-k] = read[stIntTuple_get(insert, 1)];
					insertWeight = insertWeight < stIntTuple_get(insert, 2) ? insertWeight : stIntTuple_get(insert, 2);
				}

				// Get the leftmost node in the poa graph to which the insert will connect
				PoaNode *leftNode = stList_get(poa->nodes, stIntTuple_get(insertStart, 0)+1); // Get POA node with same reference coordinate

				// Check if the complete insert is already in the poa graph:
				bool found = 0;
				for(int64_t m=0; m<stList_length(leftNode->inserts); m++) {
					PoaInsert *poaInsert = stList_get(leftNode->inserts, m);
					if(stString_eq(poaInsert->insert, insertLabel)) {
						poaInsert->weight += insertWeight;
						found = 1; break;
					}
				}
				// otherwise add the insert to the poa graph
				if(!found) {
					stList_append(leftNode->inserts, poaInsert_construct(stString_copy(insertLabel), insertWeight));
				}
			}
		}
	}

	// Add deletes to the POA graph

	// Sort the deletes first by the read coordinate and then by reference coordinate
	stList_sort(deletes, cmpDeletes);

	// Analogous to a complete-insert, let a complete-delete be a sequence of deletes with the same read coordinate j
	// and consecutive reference coordinates i, i+1, ..., i+m, such that the i,j-1 is a match
	// in the alignment subgraph and i+1,j+n is similarly a match.

	// Enumerate set of complete-deletes, adding them to the graph
	for(int64_t i=0; i<stList_length(deletes); i++) {
		stIntTuple *deleteStart = stList_get(inserts, i); // Start of putative complete-delete

		int64_t j=i+1;
		for(;j<stList_length(deletes); j++) {
			stIntTuple *deleteEnd = stList_get(deletes, j); // End of putative complete-delete

			// If they don't have the same read coordinate then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 1) != stIntTuple_get(deleteEnd, 1)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 0) + j - i != stIntTuple_get(deleteEnd, 1)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a putative maximal
		// complete-delete sequence in deletes

		// Now enumerate complete-deletes in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match then can not be a complete-delete
			if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 0) + k - i - 1, stIntTuple_get(deleteStart, 1))) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceding match then can not be a complete-delete
				if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 0) + l - i + 1, stIntTuple_get(deleteStart, 1) + 1)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-insert

				// Calculate weight
				double deleteWeight = 1.0;
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *delete = stList_get(deletes, m);
					deleteWeight = deleteWeight < stIntTuple_get(delete, 2) ? deleteWeight : stIntTuple_get(delete, 2);
				}

				// Get the leftmost node in the poa graph to which the delete will connect
				PoaNode *leftNode = stList_get(poa->nodes, stIntTuple_get(deleteStart, 0)); // Get POA node with same reference coordinate

				// Check if the delete is already in the poa graph:
				bool found = 0;
				for(int64_t m=0; m<stList_length(leftNode->deletes); m++) {
					PoaDelete *poaDelete = stList_get(leftNode->deletes, m);
					if(poaDelete->length == l-k+1) {
						poaDelete->weight += deleteWeight;
						found = 1; break;
					}
				}
				// otherwise add the delete to the poa graph
				if(!found) {
					stList_append(leftNode->deletes, poaDelete_construct(l-k+1, deleteWeight));
				}
			}
		}
	}
}

Poa *poa_realign(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p) {
	/*
	 *
	 */

	// Build a reference graph with zero weights
	Poa *poa = poa_getReferenceGraph(reference);

	// For each read
	for(int64_t i=0; i<stList_length(reads); i++) {
		char *read = stList_get(reads, i);

		// Generate set of posterior probabilities for matches, read deletes and read inserts
		stList *matches = NULL, *inserts = NULL, *deletes = NULL;
		getAlignedPairsWithIndels(sM, reference, read, p, &matches, &inserts, &deletes, 0, 0);

		// Add weights, edges and nodes to the poa
		poa_augment(poa, read, matches, inserts, deletes);

		// Cleanup
		stList_destruct(matches);
		stList_destruct(inserts);
		stList_destruct(deletes);
	}

	return poa;
}

void poa_print(Poa *poa, FILE *fH) {
	/*
	 *
	 */

	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		fprintf(fH, "%" PRIi64 "\t%c", i, node->base);

		// Bases
		double totalWeight = 0.0;
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			fprintf(fH, "\t%c:%f", symbol_convertSymbolToChar(j), node->baseWeights[j]);
			totalWeight += node->baseWeights[j];
		}
		fprintf(fH, "\tTotal-weight:%f\n", totalWeight);

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			fprintf(fH, "Insert\tSeq:%s\tWeight:%f\n", insert->insert, insert->weight);
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			fprintf(fH, "Delete\tLength:%" PRIi64 "\tWeight:%f\n", delete->length, delete->weight);
		}
	}
}

char *poa_getConsensus(Poa *poa) {
	/*
	 *
	 */
	return NULL;
}

Poa *poa_realignIterative(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p,
				 uint64_t iterations) {
	/*
	 *
	 */
	assert(iterations > 0);

	int64_t i=0;
	while(1) {
		Poa *poa = poa_realign(reads, reference, sM, p);
		if(i > 0) {
			free(reference);
		}
		if(++i >= iterations) {
			return poa;
		}
		reference = poa_getConsensus(poa);
		free(poa);
	}
}
