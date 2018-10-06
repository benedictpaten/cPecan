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
	CuAssertTrue(testCase, stList_length(node->inserts) == insertNumber);
	for(int64_t i=0; i<stList_length(node->inserts); i++) {
		PoaInsert *insert = stList_get(node->inserts, i);
		CuAssertTrue(testCase, stString_eq(inserts[i], insert->insert));
		CuAssertTrue(testCase, insert->weight == insertWeights[i]);
	}

	// Deletes
	CuAssertTrue(testCase, stList_length(node->deletes) == deleteNumber);
	for(int64_t i=0; i<stList_length(node->deletes); i++) {
		PoaDelete *delete = stList_get(node->deletes, i);
		CuAssertTrue(testCase, delete->length == deleteLengths[i]);
		CuAssertTrue(testCase, delete->weight == deleteWeights[i]);
	}
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
					2, (const char *[]){ "G", "GT" }, (const double[]){ 50.0, 50.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });
}

static void test_poa_realign_mini_example1(CuTest *testCase) {
	const char *readArray[] = {
		"CATTTTTCTCCTCCACCTGCAACAGAAGATAAAAACGCGCATCACAAACTACTTTATTG",
		"CATTTTTCTCTCCGTCACGTAATAGGAAAACAGATGAAAATGTGCACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCGTCACGACAGGAAACAGATGAAAATGGGCACAAGACCACAAACGCATTTTGAT" };
	char *reference = "CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT";
	char *trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT";

	stList *reads = stList_construct();
	for(int64_t i=0; i<45; i++) {
		stList_append(reads, (char *)readArray[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	StateMachine *sM = stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, reference, sM, p);

	fprintf(stderr, "True-reference:%s\n", trueReference);
	poa_print(poa, stderr);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
}

static void test_poa_realign_example1(CuTest *testCase) {
	const char *readArray[] = {
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
	char *reference = "CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT";
	char *trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT";

	stList *reads = stList_construct();
	for(int64_t i=0; i<45; i++) {
		stList_append(reads, (char *)readArray[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	StateMachine *sM = stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, reference, sM, p);

	fprintf(stderr, "True-reference:%s\n", trueReference);
	poa_print(poa, stderr);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
}

static void test_poa_realign_example2(CuTest *testCase) {
	const char *readArray[] = {
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
	char *reference = "GATGTAAAAAAGAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT";
	char *trueReference = "GATGTAAAAAAAAAGAAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT";

	stList *reads = stList_construct();
	for(int64_t i=0; i<42; i++) {
		stList_append(reads, (char *)readArray[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	StateMachine *sM = stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, reference, sM, p);

	fprintf(stderr, "True-reference:%s\n", trueReference);
	poa_print(poa, stderr);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
}

static void test_poa_realign(CuTest *testCase) {
	for (int64_t test = 0; test < 100; test++) {

		//Make true reference
		char *trueReference = getRandomSequence(st_randomInt(0, 100));

		// Make starting reference
		char *reference = evolveSequence(reference);

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, free);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, evolveSequence(trueReference));
		}

		PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
		StateMachine *sM = stateMachine3_construct(threeState);

		Poa *poa = poa_realign(reads, reference, sM, p);

		fprintf(stderr, "True-reference:%s\n", trueReference);
		poa_print(poa, stderr);

		//Cleanup
		stateMachine_destruct(sM);
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		pairwiseAlignmentBandingParameters_destruct(p);
		poa_destruct(poa);
	}
}

CuSuite* realignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
    SUITE_ADD_TEST(suite, test_poa_augment_example);
    SUITE_ADD_TEST(suite, test_poa_realign_mini_example1);
    SUITE_ADD_TEST(suite, test_poa_realign_example1);
    SUITE_ADD_TEST(suite, test_poa_realign_example2);
    SUITE_ADD_TEST(suite, test_poa_realign);

    return suite;
}
