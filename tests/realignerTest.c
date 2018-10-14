/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "realigner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randomSequences.h"
#include "stateMachine.h"

static char *nanoporeHmmFile = "./threeStateNanopore.hmm";

static void test_poa_getReferenceGraph(CuTest *testCase) {
	char *reference = "GATTACA";

	Poa *poa = poa_getReferenceGraph(reference);

	CuAssertTrue(testCase, stList_length(poa->nodes) == strlen(reference) + 1);
	for(int64_t i=0; i<strlen(reference); i++) {
		PoaNode *node = stList_get(poa->nodes, i+1);

		CuAssertTrue(testCase, node->base == reference[i]);
		CuAssertTrue(testCase, stList_length(node->inserts) == 0);
		CuAssertTrue(testCase, stList_length(node->deletes) == 0);
	}

	PoaNode *node = stList_get(poa->nodes, 0);
	CuAssertTrue(testCase, node->base == 'N');
	CuAssertTrue(testCase, stList_length(node->inserts) == 0);
	CuAssertTrue(testCase, stList_length(node->deletes) == 0);

	poa_destruct(poa);
}

static char *makeShiftedString(char *str, char *insert, int64_t insertPoint) {
	char *suffix = stString_copy(&str[insertPoint]);
	char *prefix = stString_copy(str);
	prefix[insertPoint] = '\0';
	char *shiftedStr = stString_print("%s%s%s", prefix, insert, suffix);
	free(suffix);
	free(prefix);
	return shiftedStr;
}

static void test_getShift(CuTest *testCase) {
	for(int64_t test=0; test<10000; test++) {
		// Make random string
		int64_t length = st_randomInt(1, 20);
		char *str = getRandomACGTSequence(length);

		// Make random insert of length m
		int64_t m = st_randomInt(1, 4);
		char *insert = getRandomACGTSequence(m);

		// Run get shift
		int64_t i = getShift(str, length, insert, m);

		//if(i + 2 < length) {
		//	fprintf(stderr, "Str: %s, str-length:%" PRIi64 " insert: %s, insert:%" PRIi64 "\n", str, length, insert, i);
		//}

		// Test resulting transplanted string is same as concatenated str+insert
		char *shiftedStr = makeShiftedString(str, insert, i);
		char *concatenatedStr = stString_print("%s%s", str, insert);

		CuAssertStrEquals(testCase, concatenatedStr, shiftedStr);

		// Cleanup
		free(shiftedStr);

		// Test no further left shift would work
		for(int64_t j=0; j<i; j++) {
			shiftedStr = makeShiftedString(str, insert, j);
			CuAssertTrue(testCase, !stString_eq(shiftedStr, concatenatedStr));
			free(shiftedStr);
		}

		// Cleanup
		free(concatenatedStr);
		free(str);
		free(insert);
	}
}


static void checkInserts(CuTest *testCase, Poa *poa, int64_t nodeIndex,
					     int64_t insertNumber, const char **inserts, const double *insertWeights, bool divideWeights) {
	PoaNode *node = stList_get(poa->nodes, nodeIndex);

	CuAssertIntEquals(testCase, stList_length(node->inserts), insertNumber);

	for(int64_t i=0; i<insertNumber; i++) {
		PoaInsert *poaInsert = stList_get(node->inserts, i);
		CuAssertStrEquals(testCase, inserts[i], poaInsert->insert);
		CuAssertDblEquals(testCase, insertWeights[i], poaInsert->weight / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkDeletes(CuTest *testCase, Poa *poa, int64_t nodeIndex,
					     int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights, bool divideWeights) {
	PoaNode *node = stList_get(poa->nodes, nodeIndex);

	CuAssertIntEquals(testCase, stList_length(node->deletes), deleteNumber);

	for(int64_t i=0; i<deleteNumber; i++) {
		PoaDelete *poaDelete = stList_get(node->deletes, i);
		CuAssertIntEquals(testCase, deleteLengths[i], poaDelete->length);
		CuAssertDblEquals(testCase, deleteWeights[i], poaDelete->weight / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkNode(CuTest *testCase, Poa *poa, int64_t nodeIndex, char base, const double *baseWeights,
		int64_t insertNumber, const char **inserts, const double *insertWeights,
		int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights) {

	PoaNode *node = stList_get(poa->nodes, nodeIndex);
	CuAssertTrue(testCase, node->base == base);

	// Matches
	for(int64_t i=0; i<SYMBOL_NUMBER; i++) {
		CuAssertDblEquals(testCase, node->baseWeights[i], baseWeights[i], 0.0);
	}

	// Inserts
	checkInserts(testCase, poa, nodeIndex, insertNumber, inserts, insertWeights, 0);

	// Deletes
	checkDeletes(testCase, poa, nodeIndex, deleteNumber, deleteLengths, deleteWeights, 0);
}

static void test_poa_augment_example(CuTest *testCase) {
	char *reference = "GATTACA";

	Poa *poa = poa_getReferenceGraph(reference);

	char *read = "GATACGGT";

	stList *matches = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
	stList *inserts = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
	stList *deletes = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

	stList_append(matches, stIntTuple_construct3(100, 0, 0));
	stList_append(matches, stIntTuple_construct3(100, 1, 1));
	stList_append(matches, stIntTuple_construct3(50, 2, 2));
	stList_append(matches, stIntTuple_construct3(50, 3, 2));
	stList_append(matches, stIntTuple_construct3(100, 4, 3));
	stList_append(matches, stIntTuple_construct3(100, 5, 4));
	stList_append(matches, stIntTuple_construct3(50, 6, 5));
	stList_append(matches, stIntTuple_construct3(25, 6, 6));
	stList_append(matches, stIntTuple_construct3(25, 6, 7));

	stList_append(inserts, stIntTuple_construct3(50, 5, 5));
	stList_append(inserts, stIntTuple_construct3(25, 5, 6));
	stList_append(inserts, stIntTuple_construct3(50, 6, 6));
	stList_append(inserts, stIntTuple_construct3(75, 6, 7));

	stList_append(deletes, stIntTuple_construct3(50, 2, 1));
	stList_append(deletes, stIntTuple_construct3(50, 3, 2));

	poa_augment(poa, read, matches, inserts, deletes);

	// Check POA graph is what we expect

	CuAssertTrue(testCase, stList_length(poa->nodes) == 8); // Length + prefix node

	checkNode(testCase, poa, 0, 'N', (const double[]){ 0.0, 0.0, 0.0, 0.0, 0.0 },
			0, (const char *[]){ "" }, (const double[]){ 0.0 },
			0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 1, 'G', (const double[]){ 0.0, 0.0, 100.0, 0.0, 0.0 },
				0, (const char *[]){ "" }, (const double[]){ 0.0 },
				0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 2, 'A', (const double[]){ 100.0, 0.0, 0.0, 0.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					1, (const int64_t[]){ 1 }, (const double[]){ 50.0 });

	checkNode(testCase, poa, 3, 'T', (const double[]){ 0.0, 0.0, 0.0, 50.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					1, (const int64_t[]){ 1 }, (const double[]){ 50.0 });

	checkNode(testCase, poa, 4, 'T', (const double[]){ 0.0, 0.0, 0.0, 50.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 5, 'A', (const double[]){ 100.0, 0.0, 0.0, 0.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 6, 'C', (const double[]){ 0.0, 100.0, 0.0, 0.0, 0.0 },
					2, (const char *[]){ "G", "GG" }, (const double[]){ 50.0, 25.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 7, 'A', (const double[]){ 0.0, 0.0, 75.0, 25.0, 0.0 },
					2, (const char *[]){ "GT", "T" }, (const double[]){ 50.0, 75.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });
}

static void test_poa_realign(CuTest *testCase) {
	for (int64_t test = 0; test < 100; test++) {

		//Make true reference
		char *trueReference = getRandomSequence(st_randomInt(1, 100));

		// Make starting reference
		char *reference = evolveSequence(trueReference);

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, free);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, evolveSequence(trueReference));
		}

		PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

		Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

		StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

		Poa *poa = poa_realign(reads, reference, sM, p);

		poa_leftAlignIndels(poa); // Shift all the indels

		// Generate the read alignments and check the matches
		// Currently don't check the insert and deletes

		double *baseWeights = st_calloc(SYMBOL_NUMBER*strlen(reference), sizeof(double));

		for(int64_t i=0; i<readNumber; i++) {
			char *read = stList_get(reads, i);

			// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
			stList *matches = NULL, *inserts = NULL, *deletes = NULL;
			getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);

			// Collate matches
			for(int64_t j=0; j<stList_length(matches); j++) {
				stIntTuple *match = stList_get(matches, j);
				baseWeights[stIntTuple_get(match, 1) * SYMBOL_NUMBER + symbol_convertCharToSymbol(read[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);
			}

			// Cleanup
			stList_destruct(matches);
			stList_destruct(inserts);
			stList_destruct(deletes);
		}

		// Check match weights tally
		for(int64_t i=0; i<strlen(reference); i++) {
			PoaNode *poaNode = stList_get(poa->nodes, i+1);
			for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
				CuAssertDblEquals(testCase, poaNode->baseWeights[j], baseWeights[i*SYMBOL_NUMBER + j], 0.0001);
			}
		}

		st_logInfo("True-reference:%s\n", trueReference);
		if (st_getLogLevel() >= info) {
			poa_print(poa, stderr);
		}

		//Cleanup
		free(baseWeights);
		stateMachine_destruct(sM);
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		pairwiseAlignmentBandingParameters_destruct(p);
		poa_destruct(poa);
		hmm_destruct(hmm);
	}
}

static void test_poa_realign_tiny_example1(CuTest *testCase) {

	char *reference = "GATACAGCGGG";
	char *read = "GATTACAGCG";

	stList *reads = stList_construct();
	stList_append(reads, read);

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	StateMachine *sM = stateMachine3_construct(threeState);

	/*
	// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
	stList *matches = NULL, *inserts = NULL, *deletes = NULL;
	getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);

	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *alignedPair = stList_get(matches, i);
		fprintf(stderr, "Match: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
				stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}

	for(int64_t i=0; i<stList_length(inserts); i++) {
		stIntTuple *alignedPair = stList_get(inserts, i);
		fprintf(stderr, "Insert: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
				stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}

	for(int64_t i=0; i<stList_length(deletes); i++) {
		stIntTuple *alignedPair = stList_get(deletes, i);
		fprintf(stderr, "Delete: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
					stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}*/

	Poa *poa = poa_realign(reads, reference, sM, p);

	// Check we get the set of inserts and deletes we expect

	// . = match
	// | = insert
	// - = delete
	// : = insert and match
	// % = delete and match

	//     Reference
	//     -1 0 1 2 3 4 5 6 7 8 9 10
	//      N G A T A C A G C G G G
	// -1 N .
	//  0 G   .
	//  1 A   | .
	//  2 T     : . -
	//  3 T       : . %
	//  4 A         :   .
	//  5 C           .   .
	//  6 A             . - %
	//  7 G               . - %
	//  8 C                 . - %
	//  9 G                   . - %

	// Check inserts

	// A after ref 0
	checkInserts(testCase, poa, 1, 1, (const char *[]){ "A" }, (const double[]){ 0.038656 }, 1);
	// T after ref 1
	checkInserts(testCase, poa, 2, 1, (const char *[]){ "T" }, (const double[]){ 0.436572 }, 1);
	// T after ref 2
	checkInserts(testCase, poa, 3, 1, (const char *[]){ "T" }, (const double[]){ 0.437963 }, 1);
	// A after ref 3
	checkInserts(testCase, poa, 4, 1, (const char *[]){ "A" }, (const double[]){ 0.038831 }, 1);

	checkInserts(testCase, poa, 0, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 5, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 6, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 7, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 8, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 9, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 10, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);

	// Check deletes

	/*
	Delete: (x:3) (y:2) (weight:0.021429)
	Delete: (x:4) (y:3) (weight:0.011958)

	Delete: (x:6) (y:6) (weight:0.039542)
	Delete: (x:7) (y:6) (weight:0.041150)

	Delete: (x:7) (y:7) (weight:0.039140)
	Delete: (x:8) (y:7) (weight:0.038841)

	Delete: (x:8) (y:8) (weight:0.440979)
	Delete: (x:9) (y:8) (weight:0.438735)

	Delete: (x:9) (y:9) (weight:0.437247)
	Delete: (x:10) (y:9) (weight:0.440302)*/

	// No deletes first three positions
	checkDeletes(testCase, poa, 0, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 1, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 2, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 5, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 10, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);

	// L1 after ref 2
	checkDeletes(testCase, poa, 3, 1, (const int64_t[]){ 1 }, (const double[]){ 0.021429 }, 1);
	// L1 after ref 3
	checkDeletes(testCase, poa, 4, 1, (const int64_t[]){ 1 }, (const double[]){ 0.011958 }, 1);
	// L2 after ref 5
	checkDeletes(testCase, poa, 6, 1, (const int64_t[]){ 2 }, (const double[]){ 0.039542 }, 1);
	// L2 after ref 6
	checkDeletes(testCase, poa, 7, 1, (const int64_t[]){ 2 }, (const double[]){ 0.038841 }, 1);
	// L2 after ref 7
	checkDeletes(testCase, poa, 8, 1, (const int64_t[]){ 2 }, (const double[]){ 0.438735 }, 1);
	// L2 after ref 8
	checkDeletes(testCase, poa, 9, 1, (const int64_t[]){ 2 }, (const double[]){ 0.437247 }, 1);

	st_logInfo("Read:%s\n", read);
	st_logInfo("Reference:%s\n", reference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr);
	}

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	stList_destruct(reads);
}

static char * runLengthEncode(char *str) {
	/*
	 * Returns a run length encoded version of the string str.
	 */

	int64_t length = strlen(str);
	char *rleStr = st_calloc(length+1, sizeof(char));

	int64_t j=0;
	for(int64_t i=0; i<length; i++) {
		if(i==0 || str[i] != str[i-1]) {
			rleStr[j++] = str[i];
		}
	}
	rleStr[j] = '\0';

	return rleStr;
}

static const char *readArrayExample1[] = {
		"CATTTTTCTCCTCCACCTGCAACAGAAGATAAAAACGCGCATCACAAACTACTTTATTG",
		"CATTTTTCTCTCCGTCACGTAATAGGAAAACAGATGAAAATGTGCACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCGTCACGACAGGAAACAGATGAAAATGGGCACAAGACCACAAACGCATTTTGAT",
		"CATTTTTCTCCGGTCATTTAATGAAAACAGATGGTACTGCGTATGTGACATAAACGCATTTTTATTT",
		"CATTTCCTCCGTCACTGCACAGGAAAACAGATGAAAATGCAAGTATGGACCCACAAAACGCATTTTATTT",
		"CATTTTTTCTCTCTCCGTCAGCTGCATTGAAAATGATGAAATGCGGGTATGACTATAAACGCATTTATTT",
		"CATTTTTTTTCTCTCCTCCACACACAGGAAACAGATGAAAAATGTATGTGACCATAAAACGCATTTTATTT",
		"TATTTTCTCCGTCATTGCAGGAAAACAGATGAAATGTAAAGTATGTGAATTACAAACGGTTTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAGTCAGATGAAAATGCGCATGTGACCATAACGCATTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATACCGTGAAACAGATGAAAAATGCGGGCATGGGACCATAAAACGCATTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAACAGATGAAAACGTGGGGCATGTGACCATAAACGCATTTTTATT",
		"CATTTTCTCTCCTCGTGTTGCACAGGAAAACAGATGAAAAATGCGAGATATGTGATCCACAAACATTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAATGATGAAAATGCGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCCTCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCATGTGACCATAAAACGCATTTTTT",
		"CATTTTCTCTCCTCCACATTGCACAGGAAAACAGATGAAAATGCGGCATGTGACCATAAAACGCATTTCTTTATTT",
		"CATTTTCTCCGTCAGTCAACAATATGAAAACAGATGAAACGCGGGCACGTGACCATAAAACGCATTTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCATTGTGGAACAGATGAAAATGCGGGGTATGTGAATCATAAAACGCATTTTATTT",
		"CATTTTTCTCTCCGTCATTGCATTAGAAAACAGGGATGAAAATGCGGGCATGTGACCATAAAAACGCATTTTTATTT",
		"CATTTTTCTCTCTCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGCGTGACTATAAAACGCATTTTTATTTT",
		"CATTTTCCTCTCCCTCCGTCATTTGCACAGGAAAACAGATGAAAAAATGCGGAATGGCTATTATAAACATTTTTAACT",
		"CATTTTTTTCTCCTCTGTCATTGCACAGGAAAACAGATGAAAAATGCGTATGTGACCATAAAATCCATTTCTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAATGATGAAAAAATGCGGGCATGTGACCATAAAACGTGCATTTTTTATTT",
		"CATTTTCTCTCTCCTCCGTGTTGCACAGGAAAACCAGATGAAAATGCGGAACATGTGTTCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCCTCTCGTCATTTGCACAGGAAGAGCAGATGAAAATGCAGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CCATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAAGCAGATGAAAAATGCGGGCATGTGACCATAAAACGCATTTTATTT",
		"CATTTTTTTCTCTCCTGTCATTGCACAGGAAACAAAGAGATGAAAAATGCGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCATTGCATAGGAAAACAGATGAAAATGCGGGGTATGTGGACCATAAAACGCATTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGGAAAACAGATGAAAATTGCGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCACTGCACAGGAAAAACAGATGAAAATGCGGGGCATGCATCATAAAACGTATTTTTATTGAATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATAAGAAAAATGCAGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTTTCACTACTCTCCCTCCGTCGTACTGGAAAACAAACAGATAAATGCAGGGCATGTGACCATAAAACATTTTTTTATTT",
		"CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATAAAAAAAAATGCAGGGGCATGTGACCATAAAACATTTTTATTT",
		"CATTTTTTCTCTCTCTCGTGTTGCACACAGGAAAACAGATGAAAAATGCCGGGGCATCATGACCATAAAACGCGTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCAGGGGCGTAACTGACCATAAAACGCATTTTTTATTT",
		"CATTTTTTCTCTCCTCCGTCATTGCACAGGAAAAATGTGATGAAAATGCGGGGTATGTGACCATAAAACGCATTTTTATGCTTCT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGAGGACATGTGACCATAAAACGCATTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCCATTGCACAGGAAAACAGATATAAAAAATGCAGGGCATCAAACCATAAAACATTTTTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCATTGCAATAGGAAAACAGATATTTTGGTGTACCGCAAGTATGTGACCATAAAACGTATTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACCAGATAGAAAAAATACAGGGCATGTGTTCATAAAGCACGCATTTTTATTT",
		"CATTTTTTTCTCTCTCCTCCGCTTTCACACACAGGAGTAAACAGATGAAAAATGTGGGCATGTGACCATAAAACGCATTTTTTATTT",
		"CATTTTCTCTCTCCCTCAAAATCATTTGCACAGGAAAACAGATAGAAAAATGCAACGGGGCATGTGATATAAAACGCATTTTTTATTT",
		"CATTTTCTACTCTCTCCCTCCGTCATTGCAGGAAAACAGATGAAAATGCAGGGAACATATATGACCATAAAACGCATTTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAAAGAGCTGGCATGCGGGGCATGTGACCATAAAACGCATTTTTTTGT" };
static char *referenceExample1 =     "CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT";
static char *trueReferenceExample1 = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT";

//                 000000   0000111111111122222222223333333333444 444444455555555556666666666777777 7777
//                 012345   6789012345678901234567890123456789012 345678901234567890123456789012345 6789
//reference =     "CATTTT   CTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAA TGCAGGGCAATAATGACCATAAAACGCATTTTT ATTT";
//trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATG  TGACCATAAAACGCATTTTTTATTT";

static void test_poa_realign_example1(CuTest *testCase) {
	stList *reads = stList_construct();
	for(int64_t i=0; i<45; i++) {
		stList_append(reads, (char *)readArrayExample1[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, referenceExample1, sM, p);

	poa_leftAlignIndels(poa); // Shift all the indels

	//st_logInfo("True-reference:%s\n", trueReference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr);
	}

	st_logInfo("True-reference:%s\n", trueReferenceExample1);
	st_logInfo("Reference:%s\n", referenceExample1);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	stList_destruct(reads);
	hmm_destruct(hmm);
}

static void test_poa_realign_rle_example1(CuTest *testCase) {
	stList *reads = stList_construct3(0, free);
	for(int64_t i=0; i<45; i++) {
		stList_append(reads, runLengthEncode((char *)readArrayExample1[i]));
	}
	char *referenceExample1RLE = runLengthEncode(referenceExample1);
	char *trueReferenceExample1RLE = runLengthEncode(trueReferenceExample1);

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, referenceExample1RLE, sM, p);

	poa_leftAlignIndels(poa); // Shift all the indels

	//st_logInfo("True-reference:%s\n", trueReference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr);
	}

	//                 **                              *    *
	//               00  00000000111111111122222222223333333333444444444455
	//               01  23456789012345678901234567890123456789012345678901
	//Reference:     CA  TCTCTCTCGTCATGCACAGACAGATGATGCAGCATATGACATACGCATAT
	//True-reference:CATCTCTCTCTCGTCATGCACAGACAGATGATGC GCATGTGACATACGCATAT

	st_logInfo("True-reference:%s\n", trueReferenceExample1RLE);
	st_logInfo("Reference:%s\n", referenceExample1RLE);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	stList_destruct(reads);
	hmm_destruct(hmm);
	free(referenceExample1RLE);
	free(trueReferenceExample1RLE);
}

static const char *readArrayExample2[] = {
			"GATGTAAAAATGACTGAGTTAGAACAGGCATAAATACATCTGT",
			"GATGTAAAAAAAAATGACAGAGAATAAAACTATCCTTATCTATT",
			"GATGTAAAAAGAAGCGGAAGTTAGAACAGGCATAAATACATCTGT",
			"GATGTAAAAAGAAATGACGGAAGAACAGAGCATAACACACATCTGT",
			"GATGTAAAAAAAGAATGATTTAGTTGAACAGAGCATAAATATCTGT",
			"GATGTAAAAAAAAGAAATGACGGAAGAACAGAGCATAACACATCTGT",
			"GATGTAAAAGAAATGGAGGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGATTTGGAAGAACAGAGCATAAATATCTGT",
			"GATGTAAAAGAAATGACGGAAGTTAGAATATATATAACACACATCTGT",
			"GATGTAAAAAAAGAATGGACGGTTAGAACAGAGCACAACACACATCTGT",
			"GATGTAAAAAAGAATGATAAAGTTAGAATAGAGCATAAATAACATCTGT",
			"GATGTAAAAGAAATGTGGAAGTTAGAACAGAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAAGAAATGAAGCTAGAACAGAGCATAAATACATCTGT",
			"GATGTAAAAAAAAAATGACCCGGAAGTTGAACAGAGCATAATACATCTGT",
			"GATGTAAAAAAGAAATGATTTAAAGGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGTGACGGAAGTTAAGACAGGCATAAATACACATCTGT",
			"GATGTAAAAAAGAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGATGGGTTAGAATAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGAAATGACGAGTTAGAACCAGAGCACCATCTACATCAT"
			"GATGTAAAAAAAAAAATGACGGAAGTTAGACAAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAGAATGATTTGAAGTTAGAACAGAGCATAACACATCTGT",
			"GATGTAAAAGAAAATCGACTGAAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGACGGAGTTAGAACAGAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAAGAAATGTGTGAGTTAAGACAGAGCATAAATACATCTAT",
			"GATGTAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTATT",
			"GATGTAAAAAAAAAAATGTGGAAGTTAAAACAGAGCATAAATACACATCTATT",
			"GATGCAAAAAAAAAGAAATGACGGAAGTTAAATTAGAGCATAAATACATCTGT",
			"GATGTAAAAAAAAGAAATGATTTGGAAGTTACAGAGCATAAATACACATCTGT",
			"GATGTAAAGAAAATGATTTTAGAAGTTAGAACAGAGCATAACACAATATCTGT",
			"GATAAAAAAAAAGGAATGATTGGAAGCTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAATGACGGAAGTTCTGAAACAGGCATAAATACACATCTGTAT",
			"GATGTAAAAAGAATGATTTGAAGTTAGAACAGAGTATATTAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTATT",
			"GATGTAAAAAAAAAGAAATGATTTGAAGCAGAACAGAGCATAAATACAAGATCTGT",
			"GATGTAAAAAAAAGAAGAAATGACGGAAGTTAAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAAAGAAATGATTGAAGTTAGAAATATACATAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGATTTTAAAGTGAACAGAGCATAAATACACACCTTGGT",
			"GATGTAAAAAAAAAAAAGAAATGACGGAAGTTGAACTAGGCTTATAAATACATCTGT",
			"GATGCCAAAAAAAAAAAGAAATGGCCAGAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGCGGATTTGGAAGTTAGAACAGTATATAAAGCACACATCCGT" };

//reference =      GATGTAAAAAA   GAAATGATTT     GCTAGAACAGAGCATAAATACACATCTGT
//trueReference =  GATGTAAAAAAAAAGAAATGA   CGGAAGTTAGAACAGAGCATAAATACACATCTGT
static char *referenceExample2 =     "GATGTAAAAAAGAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT";
static char *trueReferenceExample2 = "GATGTAAAAAAAAAGAAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT";

static void test_poa_realign_example2(CuTest *testCase) {
	stList *reads = stList_construct();
	for(int64_t i=0; i<41; i++) {
		stList_append(reads, (char *)readArrayExample2[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, referenceExample2, sM, p);

	poa_leftAlignIndels(poa);

	//st_logInfo("True-reference:%s\n", trueReference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr);
	}

	st_logInfo("True-reference:%s\n", trueReferenceExample2);
	st_logInfo("Reference:%s\n", referenceExample2);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	stList_destruct(reads);
	hmm_destruct(hmm);
}

static void test_poa_realign_rle_example2(CuTest *testCase) {
	stList *reads = stList_construct3(0, free);
	for(int64_t i=0; i<41; i++) {
		stList_append(reads, runLengthEncode((char *)readArrayExample2[i]));
	}
	char *referenceExample2RLE = runLengthEncode(referenceExample2);
	char *trueReferenceExample2RLE = runLengthEncode(trueReferenceExample2);

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, referenceExample2RLE, sM, p);

	poa_leftAlignIndels(poa);

	//st_logInfo("True-reference:%s\n", trueReference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr);
	}

	//                           * **
	//                000000000011111111112222222222333333333
	//                012345678901234567890123456789012345678
	//True-reference: GATGTAGATGACGAGTAGACAGAGCATATACACATCTGT
	//Reference:      GATGTAGATGATG CTAGACAGAGCATATACACATCTGT

	st_logInfo("True-reference:%s\n", trueReferenceExample2RLE);
	st_logInfo("Reference:%s\n", referenceExample2RLE);

	// Get alignments between true reference and original and consensus references
	// TODO

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	stList_destruct(reads);
	hmm_destruct(hmm);
	free(referenceExample2RLE);
	free(trueReferenceExample2RLE);
}

/*
 // Crufty code used to generate an initial model for nanopore alignment

typedef enum {
    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
} State;

static inline double *hmm_getTransition2(Hmm *hmm, int64_t from, int64_t to) {
    return &(hmm->transitions[from * hmm->stateNumber + to]);
}

static inline double *hmm_getEmissionsExpectation2(Hmm *hmm, int64_t state, Symbol x, Symbol y) {
    return &(hmm->emissions[state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y]);
}

static void test_hmm(CuTest *testCase) {
	Hmm *hmm = hmm_constructEmpty(0.0, threeState);

	//                 000000   0000111111111122222222223333333333444 444444455555555556666666666777777 7777
	//                 012345   6789012345678901234567890123456789012 345678901234567890123456789012345 6789
	//reference =     "CATTTT   CTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAA TGCAGGGCAATAATGACCATAAAACGCATTTTT ATTT";
	//trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATG  TGACCATAAAACGCATTTTTTATTT";

	//reference =      GATGTAAAAAA   GAAATGATTT     GCTAGAACAGAGCATAAATACACATCTGT
	//trueReference =  GATGTAAAAAAAAAGAAATGA   CGGAAGTTAGAACAGAGCATAAATACACATCTGT

	hmm_getTransition2(hmm, match, match)[0] = 0.9;
	hmm_getTransition2(hmm, match, shortGapX)[0] = 0.05;
	hmm_getTransition2(hmm, match, shortGapY)[0] = 0.05;

	hmm_getTransition2(hmm, shortGapX, match)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapX, shortGapX)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapX, shortGapY)[0] = 0.0;

	hmm_getTransition2(hmm, shortGapY, match)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapY, shortGapX)[0] = 0.0;
	hmm_getTransition2(hmm, shortGapY, shortGapY)[0] = 0.5;

	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		hmm_getEmissionsExpectation2(hmm, match, i, i)[0] = 0.98;
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			if(j != i) {
				hmm_getEmissionsExpectation2(hmm, match, i, j)[0] = 0.02/3;
			}
		}
	}

	// Gaps
	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			hmm_getEmissionsExpectation2(hmm, shortGapX, i, j)[0] = 0.025;
		}
	}

	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			hmm_getEmissionsExpectation2(hmm, shortGapY, i, j)[0] = 0.025;
		}
	}

	FILE *fH = fopen("./threeStateNanopore.hmm", "w");
	hmm_write(hmm, fH);
	fclose(fH);
}*/

CuSuite* realignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    //SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
    //SUITE_ADD_TEST(suite, test_poa_augment_example);
    //SUITE_ADD_TEST(suite, test_poa_realign_tiny_example1);
    //SUITE_ADD_TEST(suite, test_poa_realign_example1);
    //SUITE_ADD_TEST(suite, test_poa_realign_example2);
    SUITE_ADD_TEST(suite, test_poa_realign_rle_example1);
    //SUITE_ADD_TEST(suite, test_poa_realign_rle_example2);
    //SUITE_ADD_TEST(suite, test_poa_realign);
    //SUITE_ADD_TEST(suite, test_getShift);

    //SUITE_ADD_TEST(suite, test_hmm);

    return suite;
}
