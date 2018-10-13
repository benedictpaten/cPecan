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
#include <ctype.h>

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
	PoaNode *poaNode = st_calloc(1, sizeof(PoaNode));

	poaNode->inserts = stList_construct3(0, (void(*)(void *)) poaInsert_destruct);
	poaNode->deletes = stList_construct3(0, (void(*)(void *)) poaDelete_destruct);
	poaNode->base = base;
	poaNode->baseWeights = st_calloc(SYMBOL_NUMBER, sizeof(double)); // Encoded using Symbol enum

	return poaNode;
}

void poaNode_destruct(PoaNode *poaNode) {
	stList_destruct(poaNode->inserts);
	stList_destruct(poaNode->deletes);
	free(poaNode->baseWeights);
	free(poaNode);
}

Poa *poa_getEmptyGraph() {
	Poa *poa = st_calloc(1, sizeof(Poa));

	poa->nodes = stList_construct3(0, (void (*)(void *))poaNode_destruct);

	return poa;
}

void poa_destruct(Poa *poa) {
	stList_destruct(poa->nodes);
	free(poa);
}

Poa *poa_getReferenceGraph(char *reference) {
	Poa *poa = poa_getEmptyGraph();

	int64_t refLength = strlen(reference);
	stList_append(poa->nodes, poaNode_construct('N')); // Add empty prefix node
	for(int64_t i=0; i<refLength; i++) {
		stList_append(poa->nodes, poaNode_construct(toupper(reference[i])));
	}

	return poa;
}

int cmpAlignedPairsByCoordinates(const void *a, const void *b) {
	/*
	 * Compares aligned pairs, represented as stIntTuples of the form (weight, x, y) first by
	 * x coordinate and then by y coordinate, ignoring the weight.
	 */
	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1 : 0;
	}

	return i;
}

int cmpAlignedPairsByInvertedCoordinates(const void *a, const void *b) {
	/*
	 * As cmpAlignedPairsByCoordinates, but compares by y coordinate and then x coordinate.
	 */
	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	}

	return i;
}

bool isMatch(stSortedSet *matchesSet, int64_t x, int64_t y) {
	stIntTuple pair[4];
	pair[0] = 3;
	pair[2] = x;
	pair[3] = y;
	return stSortedSet_search(matchesSet, &pair) != NULL;
}

static void addToInserts(PoaNode *node, char *insert, double weight) {
	/*
	 * Add given insert to node.
	 */

	// Check if the complete insert is already in the poa graph:
	bool found = 0;
	for(int64_t m=0; m<stList_length(node->inserts); m++) {
		PoaInsert *poaInsert = stList_get(node->inserts, m);
		if(stString_eq(poaInsert->insert, insert)) {
			poaInsert->weight += weight;
			found = 1; break;
		}
	}
	// otherwise add the insert to the poa graph
	if(!found) {
		stList_append(node->inserts, poaInsert_construct(stString_copy(insert), weight));
	}
}

static void addToDeletes(PoaNode *node, int64_t length, double weight) {
	/*
	 * Add given deletion to node.
	 */

	// Check if the delete is already in the poa graph:
	bool found = 0;
	for(int64_t m=0; m<stList_length(node->deletes); m++) {
		PoaDelete *poaDelete = stList_get(node->deletes, m);
		if(poaDelete->length == length) {
			poaDelete->weight += weight;
			found = 1; break;
		}
	}
	// otherwise add the delete to the poa graph
	if(!found) {
		stList_append(node->deletes, poaDelete_construct(length, weight));
	}
}

void poa_augment(Poa *poa, char *read, stList *matches, stList *inserts, stList *deletes) {
	// Add weights of matches to the POA graph

	// For each match in alignment subgraph identify its corresponding node in the POA graph
	// add the weight of the match to the POA node
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		PoaNode *node = stList_get(poa->nodes, stIntTuple_get(match, 1)+1); // Get corresponding POA node

		// Add base weight to POA node
		node->baseWeights[symbol_convertCharToSymbol(read[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);
	}

	// Create a set of match coordinates

	stSortedSet *matchesSet = stSortedSet_construct3(cmpAlignedPairsByCoordinates, NULL);
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		assert(stSortedSet_search(matchesSet, match) == NULL);
		stSortedSet_insert(matchesSet, match); // Add to matches
	}
	
	// Add inserts to the POA graph

	// Sort the inserts first by the reference coordinate and then by read coordinate
	stList_sort(inserts, cmpAlignedPairsByCoordinates);

	// Let a complete-insert be a sequence of inserts with the same reference coordinate i
	// and consecutive read coordinates j, j+1, ..., j+n, such that the i,j-1 is a match or equal to (-1,-1) (the beginning)
	// in the alignment subgraph and i+1,j+n+1 is similarly a match or equal to (N, M) (the end of the alignment).

	// Enumerate set of complete inserts
	int64_t readLength = strlen(read), refLength = stList_length(poa->nodes)-1;
	for(int64_t i=0; i<stList_length(inserts);) {

		stIntTuple *insertStart = stList_get(inserts, i); // Start of putative complete-insert

		int64_t j=i+1;
		for(;j<stList_length(inserts); j++) {
			stIntTuple *insertEnd = stList_get(inserts, j); // End of putative complete-insert

			// If they don't have the same reference coordinate then not part of same complete-insert
			if(stIntTuple_get(insertStart, 1) != stIntTuple_get(insertEnd, 1)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(insertStart, 2) + j - i != stIntTuple_get(insertEnd, 2)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a maximal putative
		// complete insert sequence in inserts

		// Now enumerate complete-inserts in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match or the beginning then can not be a complete insert
			if(!isMatch(matchesSet, stIntTuple_get(insertStart, 1), stIntTuple_get(insertStart, 2) + k - i - 1) &&
					(stIntTuple_get(insertStart, 1) > -1 || stIntTuple_get(insertStart, 2) + k - i - 1 > -1)) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or the end then can not be a complete insert
				if(!isMatch(matchesSet, stIntTuple_get(insertStart, 1) + 1, stIntTuple_get(insertStart, 2) + l - i + 1) &&
						(stIntTuple_get(insertStart, 1) + 1 < refLength || stIntTuple_get(insertStart, 2) + l - i + 1 < readLength)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-insert

				// Calculate weight and label
				double insertWeight = UINT_MAX;
				char insertLabel[l+2-k];
				insertLabel[l+1-k] = '\0';
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *insert = stList_get(inserts, m);
					insertLabel[m-k] = toupper(read[stIntTuple_get(insert, 2)]);
					insertWeight = insertWeight < stIntTuple_get(insert, 0) ? insertWeight : stIntTuple_get(insert, 0);
				}

				// Get the leftmost node in the poa graph to which the insert will connect
				assert(stIntTuple_get(insertStart, 1) >= -1);
				PoaNode *leftNode = stList_get(poa->nodes, stIntTuple_get(insertStart, 1)+1); // Get POA node with same reference coordinate

				// Add insert to graph
				addToInserts(leftNode, insertLabel, insertWeight);
			}
		}

		// Increase i to start of next maximal complete-insert
		i = j;
	}

	// Add deletes to the POA graph

	// Sort the deletes first by the read coordinate and then by reference coordinate
	stList_sort(deletes, cmpAlignedPairsByInvertedCoordinates);

	// Analogous to a complete-insert, let a complete-delete be a sequence of deletes with the same read coordinate j
	// and consecutive reference coordinates i, i+1, ..., i+m, such that the i-1,j is a match or equal to (-1,-1) (the beginning)
	// in the alignment subgraph and i+m+1,j+1 is similarly a match or (N, M) (the alignment end).

	// Enumerate set of complete-deletes, adding them to the graph
	for(int64_t i=0; i<stList_length(deletes);) {
		stIntTuple *deleteStart = stList_get(deletes, i); // Start of putative complete-delete

		int64_t j=i+1;
		for(;j<stList_length(deletes); j++) {
			stIntTuple *deleteEnd = stList_get(deletes, j); // End of putative complete-delete

			// If they don't have the same read coordinate then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 2) != stIntTuple_get(deleteEnd, 2)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 1) + j - i != stIntTuple_get(deleteEnd, 1)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a putative maximal
		// complete-delete sequence in deletes

		// Now enumerate complete-deletes in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match or alignment beginning then can not be a complete-delete
			if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + k - i - 1, stIntTuple_get(deleteStart, 2)) &&
					((stIntTuple_get(deleteStart, 1) + k - i - 1 > -1 || stIntTuple_get(deleteStart, 2) > -1))) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or alignment end then can not be a complete-delete
				if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + l - i + 1, stIntTuple_get(deleteStart, 2) + 1) &&
					(stIntTuple_get(deleteStart, 1) + l - i + 1 < refLength || stIntTuple_get(deleteStart, 2) + 1 < readLength)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-delete

				// Calculate weight
				double deleteWeight = UINT_MAX;
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *delete = stList_get(deletes, m);
					deleteWeight = deleteWeight < stIntTuple_get(delete, 0) ? deleteWeight : stIntTuple_get(delete, 0);
				}

				// Get the leftmost node in the poa graph to which the delete will connect
				assert(stIntTuple_get(deleteStart, 1) >= 0);
				PoaNode *leftNode = stList_get(poa->nodes, stIntTuple_get(deleteStart, 1)); // Get POA node preceding delete

				// Add deletion to graph
				addToDeletes(leftNode, l-k+1, deleteWeight);
			}
		}

		// Increase i to start of next maximal complete-delete
		i = j;
	}

	// Cleanup
	stSortedSet_destruct(matchesSet);
}

Poa *poa_realign(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p) {
	// Build a reference graph with zero weights
	Poa *poa = poa_getReferenceGraph(reference);

	// For each read
	for(int64_t i=0; i<stList_length(reads); i++) {
		char *read = stList_get(reads, i);

		// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
		stList *matches = NULL, *inserts = NULL, *deletes = NULL;
		getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);

		// Add weights, edges and nodes to the poa
		poa_augment(poa, read, matches, inserts, deletes);

		// Cleanup
		stList_destruct(matches);
		stList_destruct(inserts);
		stList_destruct(deletes);
	}

	return poa;
}

char *poa_getReferenceSubstring(Poa *poa, int64_t startIndex, int64_t length) {
	/*
	 * Get substring of the reference, starting from given node.
	 */
	char *refSubString = st_malloc(sizeof(char) * (1 + length));
	refSubString[length] = '\0';
	for(int64_t j=0; j<length; j++) {
		PoaNode *node = stList_get(poa->nodes, startIndex+j);
		refSubString[j] = node->base;
	}
	return refSubString;
}

static bool matchesReferenceSubstring(char *refString, int64_t refStart, char *str, int64_t length) {
	/*
	 * Returns true if the given string str matches the given reference substring starting
	 * from the given reference position, refStart
	 */
	for(int64_t l=0; l<length; l++) {
		if(refString[refStart+l] != str[l]) {
			return 0;
		}
	}
	return 1;
}

static bool hasInternalRepeat(char *str, int64_t length, int64_t repeatLength) {
	/*
	 * Establishes if str has an internal repeat of length repeatLength.
	 * e.g. if ATATAT, internal repeat AT (length 2) is an internal repeat, but ATA is not.
	 */
	if(length % repeatLength != 0) { // If not divisible by repeatLength then can not be repeat
		return 0;
	}
	for(int64_t i=repeatLength; i<length; i+=repeatLength) {
		for(int64_t j=0; j<repeatLength; j++) {
			if(str[j] != str[j+i]) {
				return 0;
			}
		}
	}
	return 1;
}

int64_t getShift(char *refString, int64_t refStart, char *str, int64_t length) {
	// Walk back over reference sequence and see if indel can be shifted

	// Establish minimal internal repeat length
	// if ATATAT, minimal internal repeat is AT,
	// similarly if AAAAAAA then minimal internal repeat is A
	int64_t minRepeatLength = 0;
	while(minRepeatLength++ < length) {
		if(hasInternalRepeat(str, length, minRepeatLength)) {
			break;
		}
	}

	// Now walk back by multiples of minimal internal repeat length
	for(int64_t k=refStart-minRepeatLength; k>=0; k-=minRepeatLength) {
		if(!matchesReferenceSubstring(refString, k, str, minRepeatLength)) {
			break;
		}
		refStart = k;
	}

	return refStart;
}

void poa_leftAlignIndels(Poa *poa) {
	char *refString = poa_getReferenceSubstring(poa, 1, stList_length(poa->nodes)-1); // Could avoid this if
	// we store the reference string in the poa struct

	// For each poa node that is not the leftmost (because those indels can't be shifted further, by definition)
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);

		// Shift inserts

		// Process replaces old list of inserts
		stList *inserts = node->inserts;
		node->inserts = stList_construct3(0, (void (*)(void *))poaInsert_destruct);

		while(stList_length(inserts) > 0) {
			PoaInsert *insert = stList_pop(inserts);

			// Walk back over reference sequence and see if insert can be shifted
			int64_t insertPosition = getShift(refString, i, insert->insert, strlen(insert->insert));

			if(insertPosition < i) { // There is a left shift
				addToInserts(stList_get(poa->nodes, insertPosition), insert->insert, insert->weight);
				poaInsert_destruct(insert); // Cleanup old insert
			}
			else { // No left shift, just add back to filtered inserts list
				stList_append(node->inserts, insert);
			}
		}
		stList_destruct(inserts);

		// Shift deletes

		// Process replaces old list of deletes
		stList *deletes = node->deletes;
		node->deletes = stList_construct3(0, (void (*)(void *))poaDelete_destruct);

		while(stList_length(deletes) > 0) {
			PoaDelete *delete = stList_pop(deletes);

			// Get string being deleted
			char *deleteString = poa_getReferenceSubstring(poa, i+1, delete->length);

			// Walk back over reference sequence and see if delete can be shifted
			int64_t insertPosition = getShift(refString, i, deleteString, delete->length);

			free(deleteString); // Cleanup

			if(insertPosition < i) { // There is a left shift
				addToDeletes(stList_get(poa->nodes, insertPosition), delete->length, delete->weight);
				poaDelete_destruct(delete); // Cleanup old delete
			}
			else { // No left shift, just add back to filtered deletes list
				stList_append(node->deletes, delete);
			}
		}
		stList_destruct(deletes);
	}

	// Cleanup
	free(refString);
}

static int cmpInsertsByWeight(const void *a, const void *b) {
	/*
	 * Compares PoaInserts by weight in ascending order.
	 */
	PoaInsert *one = (PoaInsert *)a, *two = (PoaInsert *)b;
	return one->weight < two->weight ? -1 : one->weight > two->weight ? 1 : 0;
}

stList *poa_rankInserts(Poa *poa) {
	stList *allInserts = stList_construct();
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			stList_append(allInserts, stList_get(node->inserts, j));
		}
	}
	stList_sort(allInserts, cmpInsertsByWeight);
	return allInserts;
}

static int cmpDeletesByWeight(const void *a, const void *b) {
	/*
	 * Compares PoaDelete by weight in ascending order.
	 */
	PoaDelete *one = (PoaDelete *)a, *two = (PoaDelete *)b;
	return one->weight < two->weight ? -1 : one->weight > two->weight ? 1 : 0;
}

stList *poa_rankDeletes(Poa *poa) {
	stList *allDeletes = stList_construct();
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			stList_append(allDeletes, stList_get(node->deletes, j));
		}
	}
	stList_sort(allDeletes, cmpDeletesByWeight);
	return allDeletes;
}

void poa_print(Poa *poa, FILE *fH) {
	// Print info for each base in reference in turn
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		fprintf(fH, "%" PRIi64 "\t%c", i, node->base);

		// Bases
		double totalWeight = 0.0;
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			totalWeight += node->baseWeights[j];
		}
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			if(node->baseWeights[j]/totalWeight > 0.25) {
				fprintf(fH, "\t%c:%f (%f)", symbol_convertSymbolToChar(j), (float)node->baseWeights[j]/PAIR_ALIGNMENT_PROB_1, node->baseWeights[j]/totalWeight);
			}
		}
		fprintf(fH, "\tTotal-weight:%f\n", (float)totalWeight/PAIR_ALIGNMENT_PROB_1);

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			if(insert->weight/PAIR_ALIGNMENT_PROB_1 >= 4.0) {
				fprintf(fH, "Insert\tSeq:%s\tWeight:%f\n", insert->insert, (float)insert->weight/PAIR_ALIGNMENT_PROB_1);
			}
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			if(delete->weight/PAIR_ALIGNMENT_PROB_1 >= 4.0) {
				fprintf(fH, "Delete\tLength:%" PRIi64 "\tWeight:%f\n", delete->length, (float)delete->weight/PAIR_ALIGNMENT_PROB_1);
			}
		}
	}

	// Print top events
	fprintf(fH, "Top ranked insert events\n");
	stList *inserts = poa_rankInserts(poa);
	int64_t i=0;
	while(stList_length(inserts) > 0 && i++ < 10) {
		PoaInsert *insert = stList_pop(inserts);
		fprintf(fH, "INSERT %s %f\n", insert->insert, (float)insert->weight/PAIR_ALIGNMENT_PROB_1);
	}
	stList *deletes = poa_rankDeletes(poa);
	i=0;
	while(stList_length(deletes) > 0 && i++ < 10) {
		PoaDelete *delete = stList_pop(deletes);
		fprintf(fH, "DELETE %" PRIi64 " %f\n", delete->length, (float)delete->weight/PAIR_ALIGNMENT_PROB_1);
	}
	stList_destruct(inserts);
	stList_destruct(deletes);
}

char *poa_getConsensus(Poa *poa) {
	return NULL;
}

Poa *poa_realignIterative(stList *reads, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p,
				 uint64_t iterations) {
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
