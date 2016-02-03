#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <nanopore.h>
#include "hdp_math_utils.h"
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "continuousHmm.h"
#include "discreteHmm.h"
#include "emissionMatrix.h"
#include "multipleAligner.h"
#include "randomSequences.h"
#include "nanopore_hdp.h"


void test_first_kmer_index(CuTest* ct) {
    int64_t kmer_id = 0;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);

    for (int64_t i = 0; i < length; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }

    free(kmer);
}

void test_second_kmer_index(CuTest* ct) {
    int64_t kmer_id = 1;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);

    for (int i = 0; i < length - 1; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }
    CuAssertIntEquals(ct, 1, kmer[length - 1]);

    free(kmer);
}

void test_sixth_kmer_index(CuTest* ct) {
    int64_t kmer_id = 6;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);

    for (int64_t i = 0; i < length - 2; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }
    CuAssertIntEquals(ct, 1, kmer[length - 2]);
    CuAssertIntEquals(ct, 2, kmer[length - 1]);

    free(kmer);
}

void test_multiset_creation(CuTest* ct) {
    int64_t length = 6;
    int64_t alphabet_size = 4;
    int64_t* multiset_1 = get_word_multiset(1, alphabet_size, length);
    int64_t* multiset_2 = get_word_multiset(4, alphabet_size, length);
    int64_t* multiset_3 = get_word_multiset(16, alphabet_size, length);

    for (int64_t i = 0; i < length; i++) {
        CuAssertIntEquals(ct, multiset_1[i], multiset_2[i]);
        CuAssertIntEquals(ct, multiset_2[i], multiset_3[i]);
    }

    free(multiset_1);
    free(multiset_2);
    free(multiset_3);
}

void test_word_id_to_multiset_id(CuTest* ct) {
    int64_t length = 8;
    int64_t alphabet_size = 4;

    CuAssertIntEquals(ct, word_id_to_multiset_id(0, alphabet_size, length), 0);
    CuAssertIntEquals(ct, word_id_to_multiset_id(1, alphabet_size, length), 1);
    CuAssertIntEquals(ct, word_id_to_multiset_id(2, alphabet_size, length), 2);
    CuAssertIntEquals(ct, word_id_to_multiset_id(3, alphabet_size, length), 3);
    CuAssertIntEquals(ct, word_id_to_multiset_id(4, alphabet_size, length), 1);
    CuAssertIntEquals(ct, word_id_to_multiset_id(5, alphabet_size, length), 4);
    CuAssertIntEquals(ct, word_id_to_multiset_id(6, alphabet_size, length), 5);
    CuAssertIntEquals(ct, word_id_to_multiset_id(7, alphabet_size, length), 6);
    CuAssertIntEquals(ct, word_id_to_multiset_id(8, alphabet_size, length), 2);
    CuAssertIntEquals(ct, word_id_to_multiset_id(10, alphabet_size, length), 7);
    CuAssertIntEquals(ct, word_id_to_multiset_id(11, alphabet_size, length), 8);
    CuAssertIntEquals(ct, word_id_to_multiset_id(12, alphabet_size, length), 3);
    CuAssertIntEquals(ct, word_id_to_multiset_id(13, alphabet_size, length), 6);
    CuAssertIntEquals(ct, word_id_to_multiset_id(14, alphabet_size, length), 8);
    CuAssertIntEquals(ct, word_id_to_multiset_id(15, alphabet_size, length), 9);
    CuAssertIntEquals(ct, word_id_to_multiset_id(16, alphabet_size, length), 1);

}

void test_kmer_id(CuTest* ct) {
    CuAssertIntEquals(ct, kmer_id("AAAC", "ACGT", 4, 4), 1);
    CuAssertIntEquals(ct, kmer_id("AAAT", "ACGT", 4, 4), 3);
    CuAssertIntEquals(ct, kmer_id("AAAT", "ACT", 3, 4), 2);
    CuAssertIntEquals(ct, kmer_id("GGGG", "ABCDEFG", 7, 4), power(7, 4) - 1);
    CuAssertIntEquals(ct, standard_kmer_id("AACAA", 5), 16);
}

void add_hdp_copy_tests(CuTest* ct, HierarchicalDirichletProcess* original_hdp,
                        HierarchicalDirichletProcess* copy_hdp) {

    CuAssertIntEquals_Msg(ct, "struct finalized fail, check output .hdp files",
                          (int) is_structure_finalized(original_hdp), (int) is_structure_finalized(copy_hdp));
    CuAssertIntEquals_Msg(ct, "sampling gamma fail, check output .hdp files",
                          (int) is_gamma_random(original_hdp), (int) is_gamma_random(copy_hdp));
    CuAssertIntEquals_Msg(ct, "distr finalized fail, check output .hdp files",
                          (int) is_sampling_finalized(original_hdp), (int) is_sampling_finalized(copy_hdp));
    CuAssertIntEquals_Msg(ct,  "num dir proc fail, check output .hdp files",
                          (int) get_num_dir_proc(original_hdp), (int) get_num_dir_proc(copy_hdp));
    CuAssertIntEquals_Msg(ct, "depth fail, check output .hdp files",
                          (int) get_depth(original_hdp), (int) get_depth(copy_hdp));
    CuAssertIntEquals_Msg(ct, "num data fail, check output .hdp files",
                          (int) get_num_data(original_hdp), (int) get_num_data(copy_hdp));
    CuAssertIntEquals_Msg(ct, "grid length fail, check output .hdp files",
                          (int) get_grid_length(original_hdp), (int) get_grid_length(copy_hdp));

    CuAssertDblEquals_Msg(ct, "mu fail, check output .hdp files",
                          get_mu(original_hdp), get_mu(copy_hdp), 0.000001);
    CuAssertDblEquals_Msg(ct, "nu fail, check output .hdp files",
                          get_nu(original_hdp), get_nu(copy_hdp), 0.0000001);
    CuAssertDblEquals_Msg(ct, "alpha fail, check output .hdp files",
                          get_alpha(original_hdp), get_alpha(copy_hdp), 0.000001);
    CuAssertDblEquals_Msg(ct, "beta fail, check output .hdp files",
                          get_beta(original_hdp), get_beta(copy_hdp), 0.000001);

    int64_t num_data = get_num_data(original_hdp);
    int64_t grid_length = get_grid_length(original_hdp);
    int64_t num_dps = get_num_dir_proc(original_hdp);

    double* original_data = get_data_copy(original_hdp);
    double* copy_data = get_data_copy(copy_hdp);

    int64_t* original_dp_ids = get_data_pt_dp_ids_copy(original_hdp);
    int64_t* copy_dp_ids = get_data_pt_dp_ids_copy(copy_hdp);

    for (int64_t i = 0; i < num_data; i++) {
        CuAssertIntEquals_Msg(ct, "data pt dp id fail, check output .hdp files",
                              (int) original_dp_ids[i], (int) copy_dp_ids[i]);
        CuAssertDblEquals_Msg(ct, "data pt fail, check output .hdp files",
                              original_data[i], copy_data[i], 0.000001);
    }

    free(original_data);
    free(copy_data);
    free(original_dp_ids);
    free(copy_dp_ids);

    int64_t original_num_dps;
    int64_t* original_num_dp_fctrs;
    int64_t original_num_gamma_params;
    double* original_gamma_params;
    double original_log_likelihood;
    double original_log_density;

    int64_t copy_num_dps;
    int64_t* copy_num_dp_fctrs;
    int64_t copy_num_gamma_params;
    double* copy_gamma_params;
    double copy_log_likelihood;
    double copy_log_density;

    take_snapshot(original_hdp, &original_num_dp_fctrs, &original_num_dps, &original_gamma_params,
                  &original_num_gamma_params, &original_log_likelihood, &original_log_density);

    take_snapshot(copy_hdp, &copy_num_dp_fctrs, &copy_num_dps, &copy_gamma_params,
                  &copy_num_gamma_params, &copy_log_likelihood, &copy_log_density);


    CuAssertIntEquals_Msg(ct, "num dps snapshot fail, check output .hdp files",
                          (int) original_num_dps, (int) copy_num_dps);
    CuAssertIntEquals_Msg(ct, "num gam params snapshot fail, check output .hdp files",
                          (int) original_num_gamma_params, (int) original_num_gamma_params);

    for (int64_t i = 0; i < original_num_dps; i++) {
        CuAssertIntEquals_Msg(ct, "num dp fctrs fail, check output .hdp files",
                              (int) original_num_dp_fctrs[i], (int) copy_num_dp_fctrs[i]);
    }

    for (int64_t i = 0; i < original_num_gamma_params; i++) {
        CuAssertDblEquals_Msg(ct, "gamma param fail, check output .hdp files",
                              original_gamma_params[i], copy_gamma_params[i], 0.000001);
    }

    CuAssertDblEquals_Msg(ct, "log likelihood fail, check output .hdp files",
                          original_log_likelihood, copy_log_likelihood, 0.000001);

    CuAssertDblEquals_Msg(ct, "log density fail, check output .hdp files",
                          original_log_density, copy_log_density, 0.000001);

    free(original_num_dp_fctrs);
    free(original_gamma_params);
    free(copy_num_dp_fctrs);
    free(copy_gamma_params);

    double* original_grid = get_sampling_grid_copy(original_hdp);
    double* copy_grid = get_sampling_grid_copy(copy_hdp);

    for (int64_t i = 0; i < grid_length; i++) {
        CuAssertDblEquals_Msg(ct, "grid fail, check output .hdp files",
                              original_grid[i], copy_grid[i], 0.000001);
    }

    if (is_sampling_finalized(original_hdp) && is_sampling_finalized(copy_hdp)) {
        int64_t test_grid_length = grid_length / 2;
        double* test_grid = linspace(1.1 * original_grid[0],
                                     1.1 * original_grid[grid_length - 1], test_grid_length);

        for (int64_t i = 0; i < test_grid_length; i++) {
            for (int64_t id = 0; id < num_dps; id++) {
                CuAssertDblEquals_Msg(ct, "dp density fail, check output .hdp files",
                                      dir_proc_density(original_hdp, test_grid[i], id),
                                      dir_proc_density(copy_hdp, test_grid[i], id),
                                      0.000001);
            }
        }

        free(test_grid);
    }

    free(original_grid);
    free(copy_grid);

    for (int64_t id = 0; id < num_dps; id++) {
        CuAssertIntEquals_Msg(ct, "dp size fail, check output .hdp files",
                              (int) get_dir_proc_num_factors(original_hdp, id),
                              (int) get_dir_proc_num_factors(copy_hdp, id));
        CuAssertIntEquals_Msg(ct, "dp parent fail, check output .hdp files",
                              (int) get_dir_proc_parent_id(original_hdp, id),
                              (int) get_dir_proc_parent_id(copy_hdp, id));
    }
}

void test_serialization(CuTest* ct) {

    FILE* data_file = fopen("../../cPecan/tests/test_hdp/data.txt","r");
    FILE* dp_id_file = fopen("../../cPecan/tests/test_hdp/dps.txt", "r");

    stList* data_list = stList_construct3(0, free);
    stList* dp_id_list = stList_construct3(0, free);

    char* data_line = stFile_getLineFromFile(data_file);
    char* dp_id_line = stFile_getLineFromFile(dp_id_file);

    double* datum_ptr;
    int64_t* dp_id_ptr;
    while (data_line != NULL) {

        datum_ptr = (double*) malloc(sizeof(double));
        dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));

        sscanf(data_line, "%lf", datum_ptr);
        sscanf(dp_id_line, "%"SCNd64, dp_id_ptr);

        if (*dp_id_ptr != 4) {
            stList_append(data_list, datum_ptr);
            stList_append(dp_id_list, dp_id_ptr);
        }

        free(data_line);
        data_line = stFile_getLineFromFile(data_file);

        free(dp_id_line);
        dp_id_line = stFile_getLineFromFile(dp_id_file);
    }

    int64_t data_length;
    int64_t dp_ids_length;

    double* data = stList_toDoublePtr(data_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);

    fclose(data_file);
    fclose(dp_id_file);

    int64_t num_dir_proc = 8;
    int64_t depth = 3;

    double mu = 0.0;
    double nu = 1.0;
    double alpha = 2.0;
    double beta = 10.0;

    int64_t grid_length = 250;
    double grid_start = -10.0;
    double grid_end = 10.0;

    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;

    HierarchicalDirichletProcess* original_hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
                                                                     gamma_beta, grid_start, grid_end,
                                                                     grid_length, mu, nu, alpha, beta);


    set_dir_proc_parent(original_hdp, 1, 0);
    set_dir_proc_parent(original_hdp, 2, 0);
    set_dir_proc_parent(original_hdp, 3, 1);
    set_dir_proc_parent(original_hdp, 4, 1);
    set_dir_proc_parent(original_hdp, 5, 1);
    set_dir_proc_parent(original_hdp, 6, 2);
    set_dir_proc_parent(original_hdp, 7, 2);
    finalize_hdp_structure(original_hdp);

    char* filepath = "../../cPecan/tests/test_hdp/test.hdp";
    char* copy_filepath = "../../cPecan/tests/test_hdp/test_copy.hdp";

    FILE* main_file;
    FILE* copy_file;

    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    HierarchicalDirichletProcess* copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    pass_data_to_hdp(original_hdp, data, dp_ids, data_length);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    execute_gibbs_sampling(original_hdp, 10, 10, 10, false);
    finalize_distributions(original_hdp);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    destroy_hier_dir_proc(original_hdp);

    data = stList_toDoublePtr(data_list, &data_length);
    dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);

    stList_destruct(dp_id_list);
    stList_destruct(data_list);

    double* gamma_params = (double*) malloc(sizeof(double) * depth);
    gamma_params[0] = 1.0; gamma_params[1] = 1.0; gamma_params[2] = 2.0;

    original_hdp = new_hier_dir_proc(num_dir_proc, depth, gamma_params, grid_start,
                                     grid_end, grid_length, mu, nu, alpha, beta);


    set_dir_proc_parent(original_hdp, 1, 0);
    set_dir_proc_parent(original_hdp, 2, 0);
    set_dir_proc_parent(original_hdp, 3, 1);
    set_dir_proc_parent(original_hdp, 4, 1);
    set_dir_proc_parent(original_hdp, 5, 1);
    set_dir_proc_parent(original_hdp, 6, 2);
    set_dir_proc_parent(original_hdp, 7, 2);
    finalize_hdp_structure(original_hdp);

    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    pass_data_to_hdp(original_hdp, data, dp_ids, data_length);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    execute_gibbs_sampling(original_hdp, 10, 10, 10, false);
    finalize_distributions(original_hdp);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);

    destroy_hier_dir_proc(original_hdp);

    remove(filepath);
}

void test_nhdp_serialization(CuTest* ct) {

    NanoporeHDP* nhdp = flat_hdp_model("ACGT", 4, 6, 4.0, 20.0, 0.0, 100.0, 100,
                                       "../../cPecan/models/template_median68pA.model");

    update_nhdp_from_alignment(nhdp, "../../cPecan/tests/test_alignments/simple_alignment.tsv", false);

    execute_nhdp_gibbs_sampling(nhdp, 10, 0, 1, false);
    finalize_nhdp_distributions(nhdp);

    serialize_nhdp(nhdp, "../../cPecan/tests/test_hdp/test.nhdp");
    NanoporeHDP* copy_nhdp = deserialize_nhdp("../../cPecan/tests/test_hdp/test.nhdp");

    double fakemean[1] = {65.0};
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, "AAAAAA", fakemean),
                          get_nanopore_kmer_density(copy_nhdp, "AAAAAA", fakemean),
                          0.000001);

    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, "GCACCA", fakemean),
                          get_nanopore_kmer_density(copy_nhdp, "GCACCA", fakemean),
                          0.000001);

    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, "CCTTAG", fakemean),
                          get_nanopore_kmer_density(copy_nhdp, "CCTTAG", fakemean),
                          0.000001);

    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, "ACTTCA", fakemean),
                          get_nanopore_kmer_density(copy_nhdp, "ACTTCA", fakemean),
                          0.000001);

    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, "GGAATC", fakemean),
                          get_nanopore_kmer_density(copy_nhdp, "GGAATC", fakemean),
                          0.000001);

    remove("../../cPecan/tests/test_hdp/test.nhdp");
}


static void test_sm3hdp_cell(CuTest *testCase) {
    // load model and make stateMachine
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    // make nanopore HDP
    NanoporeHDP *nHdp = flat_hdp_model("ACGT", 4, 6, 4.0, 20.0, 0.0, 100.0, 100,
                                       "../../cPecan/models/template_median68pA.model");
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 100, 0, 1, FALSE);
    finalize_nhdp_distributions(nHdp);

    StateMachine *sM = getHdpStateMachine3(nHdp);

    double lowerF[sM->stateNumber], middleF[sM->stateNumber], upperF[sM->stateNumber], currentF[sM->stateNumber];
    double lowerB[sM->stateNumber], middleB[sM->stateNumber], upperB[sM->stateNumber], currentB[sM->stateNumber];
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        middleF[i] = sM->startStateProb(sM, i);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = sM->endStateProb(sM, i);
    }

    char *referenceSeq = "ATGACACATT";
    double fakeEventSeq[15] = {
            60.032615, 0.791316, 0.005, //ATGACA
            60.332089, 0.620198, 0.012, //TGACAC
            61.618848, 0.747567, 0.008, //GACACA
            66.015805, 0.714290, 0.021, //ACACAT
            59.783408, 1.128591, 0.002, //CACATT
    };

    int64_t testLength = 5;

    int64_t correctedLength = sequence_correctSeqLength(strlen(referenceSeq), event);
    CuAssertIntEquals(testCase, testLength, correctedLength);

    // make sequence objects
    Sequence *eventSeq = sequence_construct(testLength, fakeEventSeq, sequence_getEvent);
    Sequence *referSeq = sequence_construct(correctedLength, referenceSeq, sequence_getKmer);

    // test sequence_getEvent
    for (int64_t i = 0; i < testLength; i++) {
        CuAssertDblEquals(testCase, *(double *)eventSeq->get(eventSeq->elements, i),
                          fakeEventSeq[i * NB_EVENT_PARAMS], 0.0);
    }

    // get one element from each sequence
    void *kX = referSeq->get(referSeq->elements, 1);
    void *eY = eventSeq->get(eventSeq->elements, 1);

    //Do forward
    cell_calculateForward(sM, lowerF, NULL, NULL, middleF, kX, eY, NULL);
    cell_calculateForward(sM, upperF, middleF, NULL, NULL, kX, eY, NULL);
    cell_calculateForward(sM, currentF, lowerF, middleF, upperF, kX, eY, NULL);

    //Do backward
    cell_calculateBackward(sM, currentB, lowerB, middleB, upperB, kX, eY, NULL);
    cell_calculateBackward(sM, upperB, middleB, NULL, NULL, kX, eY, NULL);
    cell_calculateBackward(sM, lowerB, NULL, NULL, middleB, kX, eY, NULL);
    double totalProbForward = cell_dotProduct2(currentF, sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, sM, sM->startStateProb);
    //st_uglyf("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);

    //Check the forward and back probabilities are about equal
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001);

    // cleanup
    sequence_sequenceDestroy(eventSeq);
    sequence_sequenceDestroy(referSeq);
}

static void test_sm3Hdp_dpDiagonal(CuTest *testCase) {
    // load model and make stateMachine
    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    // make nanopore HDP
    char *canonicalAlphabet = "ACGT\0";

    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, modelFile);
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 200, 10000, 50, FALSE);
    finalize_nhdp_distributions(nHdp);
    StateMachine *sM = getHdpStateMachine3(nHdp);

    Diagonal diagonal = diagonal_construct(3, -1, 1);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, sM->stateNumber);

    //Get cell
    double *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);

    double *c2 = dpDiagonal_getCell(dpDiagonal, 1);
    CuAssertTrue(testCase, c2 != NULL);

    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, sM, sM->endStateProb); //Test initialise values
    double totalProb = LOG_ZERO;
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        CuAssertDblEquals(testCase, c1[i], sM->endStateProb(sM, i), 0.0);
        CuAssertDblEquals(testCase, c2[i], sM->endStateProb(sM, i), 0.0);
        totalProb = logAdd(totalProb, 2 * c1[i]);
        totalProb = logAdd(totalProb, 2 * c2[i]);
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));

    //Check it runs
    CuAssertDblEquals(testCase, totalProb, dpDiagonal_dotProduct(dpDiagonal, dpDiagonal2), 0.001);

    // cleanup
    stateMachine_destruct(sM);
    dpDiagonal_destruct(dpDiagonal);
    dpDiagonal_destruct(dpDiagonal2);
}

static void test_sm3Hdp_diagonalDPCalculations(CuTest *testCase) {
    // make some DNA sequences and fake nanopore read data
    char *sX = "ACGATACGGACAT";
    double sY[21] = {
            58.743435, 0.887833, 0.0571, //ACGATA 0
            53.604965, 0.816836, 0.0571, //CGATAC 1
            58.432015, 0.735143, 0.0571, //GATACG 2
            63.684352, 0.795437, 0.0571, //ATACGG 3
            //63.520262, 0.757803, 0.0571, //TACGGA skip
            58.921430, 0.812959, 0.0571, //ACGGAC 4
            59.895882, 0.740952, 0.0571, //CGGACA 5
            61.684303, 0.722332, 0.0571, //GGACAT 6
    };

    // make variables for the (corrected) length of the sequences
    int64_t lX = sequence_correctSeqLength(strlen(sX), event);
    int64_t lY = 7;

    // make Sequence objects
    Sequence *SsX = sequence_construct(lX, sX, sequence_getKmer3);
    Sequence *SsY = sequence_construct(lY, sY, sequence_getEvent);

    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    // make nanopore HDP
    char *canonicalAlphabet = "ACGT\0";

    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, modelFile);
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 200, 10000, 50, FALSE);
    finalize_nhdp_distributions(nHdp);

    StateMachine *sM = getHdpStateMachine3(nHdp);

    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, SsX->length, SsY->length, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    // Initialize Matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);

    //Forward algorithm
    for (int64_t i = 1; i <= lX + lY; i++) {
        diagonalCalculationForward(sM, i, dpMatrixForward, SsX, SsY);
    }

    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        diagonalCalculationBackward(sM, i, dpMatrixBackward, SsX, SsY);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM, sM->startStateProb);
    //st_uglyf("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    // Test the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i,
                                                                       dpMatrixForward,
                                                                       dpMatrixBackward,
                                                                       SsX, SsY);
        //Check the forward and back probabilities are about equal
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01);
    }

    // Now do the posterior probabilities, get aligned pairs with posterior match probs above threshold
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.2;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, SsX, SsY,
                                               totalProbForward, p, extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    // Make a list of the correct anchor points
    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                          (void (*)(void *)) stIntTuple_destruct);

    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 2));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(3, 3));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(4, 3));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(5, 4));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(6, 5));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(7, 6));

    // make sure alignedPairs is correct
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        //st_uglyf("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        //CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2(x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, 12, (int) stList_length(alignedPairs));

    // clean up
    stateMachine_destruct(sM);
    sequence_sequenceDestroy(SsX);
    sequence_sequenceDestroy(SsY);
}

static void checkAlignedPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
                                                (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 3);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}
// alignment without banding

static void test_sm3Hdp_getAlignedPairsWithBanding(CuTest *testCase) {
    // load the reference sequence and the nanopore read
    char *ZymoReference = stString_print("../../cPecan/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(ZymoReference, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    char *npReadFile = stString_print("../../cPecan/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    nanopore_descaleEvents(npRead->nbTemplateEvents, npRead->templateEvents,
                           npRead->templateParams.scale,
                           npRead->templateParams.shift);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    // make nanopore HDP
    char *canonicalAlphabet = "ACGT\0";
    // todo try other hdp models
    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, modelFile);
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 200, 10000, 50, FALSE);
    finalize_nhdp_distributions(nHdp);
    StateMachine *sMt = getHdpStateMachine3(nHdp);


    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->threshold = 0.1;

    // get anchors using lastz
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

    // remap and filter
    stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
    stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

    // make Sequences for reference and template events
    Sequence *refSeq = sequence_construct2(lX, ZymoReferenceSeq, sequence_getKmer3,
                                           sequence_sliceNucleotideSequence2);
    Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                sequence_sliceEventSequence2);

    // do alignment of template events
    stList *alignedPairs = getAlignedPairsUsingAnchors(sMt, refSeq, templateSeq, filteredRemappedAnchors, p,
                                                       diagonalCalculationPosteriorMatchProbs,
                                                       0, 0);
    checkAlignedPairs(testCase, alignedPairs, lX, lY);

    // for ch1_file1 template there should be this many aligned pairs with banding
    //st_uglyf("got %lld alignedPairs with anchors\n", stList_length(alignedPairs));
    CuAssertTrue(testCase, stList_length(alignedPairs) == 2887);

    // clean
    pairwiseAlignmentBandingParameters_destruct(p);
    nanopore_nanoporeReadDestruct(npRead);
    sequence_sequenceDestroy(refSeq);
    sequence_sequenceDestroy(templateSeq);
    stList_destruct(alignedPairs);
    stateMachine_destruct(sMt);
}

static void test_hdpHmm(CuTest *testCase) {
    // make the hmm object
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";
    // make nanopore HDP
    NanoporeHDP *nHdp = flat_hdp_model("ACGT", 4, 6, 4.0, 20.0, 0.0, 100.0, 100,
                                       "../../cPecan/models/template_median68pA.model");
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 100, 0, 1, FALSE);
    finalize_nhdp_distributions(nHdp);

    Hmm *hmm = hdpHmm_constructEmpty(0.0, 3, NUM_OF_KMERS, threeState_hdp, 0.02,
                                        continuousPairHmm_addToTransitionsExpectation,
                                        continuousPairHmm_setTransitionExpectation,
                                        continuousPairHmm_getTransitionExpectation,
                                        continuousPairHmm_addToKmerGapExpectation,
                                        continuousPairHmm_setKmerGapExpectation,
                                        continuousPairHmm_getKmerGapExpectation,
                                        emissions_discrete_getKmerIndexFromKmer);
    HdpHmm *hdpHmm = (HdpHmm *) hmm;

    // Add some transition expectations
    int64_t nStates = hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.stateNumber;

    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double dummy = from * nStates + to;
            hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.addToTransitionExpectationFcn(
                    (Hmm *) hdpHmm, from, to, dummy);
        }
    }

    // Add some emissions to the kmer skip probs
    int64_t nSymbols = hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.symbolSetSize;
    double dummyTotal = 0.0;
    for (int64_t i = 0; i < nSymbols; i++) {
        double dummy = (nSymbols * nStates) + i;
        hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.setEmissionExpectationFcn(
                (Hmm *) hdpHmm, 0, i, 0, dummy);
        dummyTotal += dummy;
    }

    char *sequence = "ACGTCATACATGACTATA";
    double fakeMeans[3] = { 65.0, 64.0, 63.0 };
    for (int64_t a = 0; a < 3; a++) {
        hdpHmm->addToAssignments((Hmm *)hdpHmm, sequence + (a * KMER_LENGTH), fakeMeans + a);
    }
    CuAssertTrue(testCase, hdpHmm->numberOfAssignments == 3);

    // dump the HMM to a file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile)); //Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    hdpHmm_writeToFile((Hmm *) hdpHmm, fH);
    fclose(fH);
    hdpHmm_destruct((Hmm *) hdpHmm);
    /*

    //Load from a file
    hmm = continuousPairHmm_loadFromFile(tempFile);
    hdpHmm = (HdpHmm *)hmm;
    stFile_rmrf(tempFile);

    //CuAssertTrue(testCase, (nb_matches == hdpHmm->baseContinuousHmm.numberOfAssignments));

    // Check the transition expectations
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double retrievedProb = hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.getTransitionsExpFcn(
                    (Hmm *) hdpHmm, from, to);
            double correctProb = from * nStates + to;
            CuAssertTrue(testCase, retrievedProb == correctProb);
        }
    }

    // check the kmer skip probs
    for (int64_t i = 0; i < nSymbols; i++) {
        double received = hdpHmm->baseContinuousPairHmm.baseContinuousHmm.baseHmm.getEmissionExpFcn(
                (Hmm *) hdpHmm, 0, i, 0);
        double correct = (nSymbols * nStates) + i;
        CuAssertDblEquals(testCase, received, correct, 0.0);
    }

    // check the assignments
    //for (int64_t a = 0; a < hdpHmm->baseContinuousHmm.numberOfAssignments; a++) {
    //    stDoubleTuple *assignment = stList_get(hdpHmm->baseContinuousHmm.assignments, a);
    //    double retrievedMean = stDoubleTuple_getPosition(assignment, 0);
    //    int64_t retrievedKmerIdx = (int64_t )stDoubleTuple_getPosition(assignment, 1);
    //    CuAssertDblEquals(testCase, retrievedMean, ((a + 1) * b), 0.001);
    //    CuAssertTrue(testCase, retrievedKmerIdx == a);
    //}

    // normalize
    continuousPairHmm_normalize((Hmm *) hdpHmm);

    // Recheck transitions
    for (int64_t from = 0; from < nStates; from++) {
        for (int64_t to = 0; to < nStates; to++) {
            double z = from * nStates * nStates + (nStates * (nStates - 1)) / 2;
            double retrievedNormedProb = hdpHmm->baseContinuousHmm.baseHmm.getTransitionsExpFcn((Hmm *) hdpHmm, from, to);
            CuAssertDblEquals(testCase, (from * nStates + to) / z, retrievedNormedProb, 0.0);
        }
    }

    // Recheck kmerGap emissions
    for (int64_t i = 0; i < nSymbols; i++) {
        double received = hdpHmm->baseContinuousHmm.baseHmm.getEmissionExpFcn((Hmm *) hdpHmm, 0, i, 0);
        double correctNormed = ((nSymbols * nStates) + i) / dummyTotal;
        CuAssertDblEquals(testCase, received, correctNormed, 0.001);
    }
    continuousPairHmm_destruct((Hmm *) hdpHmm);
    */
}

static void test_hdpHmm_em(CuTest *testCase) {
    // load the reference sequence
    char *referencePath = stString_print("../../cPecan/tests/test_npReads/ZymoRef.txt");
    FILE *fH = fopen(referencePath, "r");
    char *ZymoReferenceSeq = stFile_getLineFromFile(fH);
    // parameters for pairwise alignment using defaults
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    // load the npRead
    char *npReadFile = stString_print("../../cPecan/tests/test_npReads/ZymoC_ch_1_file1.npRead");
    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(npReadFile);

    nanopore_descaleEvents(npRead->nbTemplateEvents, npRead->templateEvents,
                           npRead->templateParams.scale,
                           npRead->templateParams.shift);

    // get sequence lengths
    int64_t lX = sequence_correctSeqLength(strlen(ZymoReferenceSeq), event);
    int64_t lY = npRead->nbTemplateEvents;

    // start with random model
    double pLikelihood = -INFINITY;
    Hmm *hdpHmm = hdpHmm_constructEmpty(0.0, 3, NUM_OF_KMERS, threeState_hdp, p->threshold,
                                        continuousPairHmm_addToTransitionsExpectation,
                                        continuousPairHmm_setTransitionExpectation,
                                        continuousPairHmm_getTransitionExpectation,
                                        continuousPairHmm_addToKmerGapExpectation,
                                        continuousPairHmm_setKmerGapExpectation,
                                        continuousPairHmm_getKmerGapExpectation,
                                        emissions_discrete_getKmerIndexFromKmer);
    continuousPairHmm_randomize(hdpHmm);

    char *modelFile = stString_print("../../cPecan/models/template_median68pA.model");
    char *alignmentFile = stString_print("../../cPecan/tests/test_alignments/simple_alignment.tsv");
    char *strand = "t";

    // make nanopore HDP
    char *canonicalAlphabet = "ACGT\0";
    NanoporeHDP *nHdp = flat_hdp_model(canonicalAlphabet, 4, KMER_LENGTH, 1.0, 1.0, 40.0, 100.0, 100, modelFile);
    update_nhdp_from_alignment_with_filter(nHdp, alignmentFile, FALSE, strand);
    execute_nhdp_gibbs_sampling(nHdp, 200, 10000, 50, FALSE);
    finalize_nhdp_distributions(nHdp);
    StateMachine *sMt = getHdpStateMachine3(nHdp);

    // load (random) transitions into stateMachine
    continuousPairHmm_loadTransitionsAndKmerGapProbs(sMt, hdpHmm);

    // close hmm
    continuousPairHmm_destruct(hdpHmm);

    for (int64_t iter = 0; iter < 10; iter++) {
        hdpHmm = continuousPairHmm_constructEmpty(0.0, 3, NUM_OF_KMERS, threeState,
                                                  continuousPairHmm_addToTransitionsExpectation,
                                                  continuousPairHmm_setTransitionExpectation,
                                                  continuousPairHmm_getTransitionExpectation,
                                                  continuousPairHmm_addToKmerGapExpectation,
                                                  continuousPairHmm_setKmerGapExpectation,
                                                  continuousPairHmm_getKmerGapExpectation,
                                                  emissions_discrete_getKmerIndexFromKmer);
        // E step
        // get anchors using lastz
        stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(ZymoReferenceSeq, npRead->twoDread, p);

        // remap and filter
        stList *remappedAnchors = nanopore_remapAnchorPairs(anchorPairs, npRead->templateEventMap);
        stList *filteredRemappedAnchors = filterToRemoveOverlap(remappedAnchors);

        // make Sequences for reference and template events
        Sequence *refSeq = sequence_construct2(lX, ZymoReferenceSeq, sequence_getKmer3,
                                               sequence_sliceNucleotideSequence2);
        Sequence *templateSeq = sequence_construct2(lY, npRead->templateEvents, sequence_getEvent,
                                                    sequence_sliceEventSequence2);

        getExpectationsUsingAnchors(sMt, hdpHmm, refSeq, templateSeq, filteredRemappedAnchors,
                                    p, diagonalCalculation_signal_Expectations, 0, 0);

        continuousPairHmm_normalize(hdpHmm);

        //Log stuff
        for (int64_t from = 0; from < sMt->stateNumber; from++) {
            for (int64_t to = 0; to < sMt->stateNumber; to++) {
                st_logInfo("Transition from %" PRIi64 " to %" PRIi64 " has expectation %f\n", from, to,
                           hdpHmm->getTransitionsExpFcn(hdpHmm, from, to));
            }
        }
        for (int64_t x = 0; x < hdpHmm->symbolSetSize; x++) {
            st_logInfo("Emission x %" PRIi64 " has expectation %f\n", x, hdpHmm->getEmissionExpFcn(hdpHmm, 0, x, 0));

        }
        //st_uglyf("->->-> Got expected likelihood %f for iteration %" PRIi64 "\n", hdpHmm->likelihood, iter);

        // M step
        // make new HDP here

        // then load transitions
        continuousPairHmm_loadTransitionsAndKmerGapProbs(sMt, hdpHmm);

        // Tests
        assert(pLikelihood <= hdpHmm->likelihood * 0.95);
        CuAssertTrue(testCase, pLikelihood <= hdpHmm->likelihood * 0.95);
        // update
        pLikelihood = hdpHmm->likelihood;

        // per iteration clean up
        continuousPairHmm_destruct(hdpHmm);
        sequence_sequenceDestroy(refSeq);
        sequence_sequenceDestroy(templateSeq);
        stList_destruct(filteredRemappedAnchors);
    }
    nanopore_nanoporeReadDestruct(npRead);
    pairwiseAlignmentBandingParameters_destruct(p);
    stateMachine_destruct(sMt);
}


CuSuite *NanoporeHdpTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_first_kmer_index);
    SUITE_ADD_TEST(suite, test_second_kmer_index);
    SUITE_ADD_TEST(suite, test_sixth_kmer_index);
    SUITE_ADD_TEST(suite, test_multiset_creation);
    SUITE_ADD_TEST(suite, test_word_id_to_multiset_id);
    SUITE_ADD_TEST(suite, test_kmer_id);
    SUITE_ADD_TEST(suite, test_serialization);
    SUITE_ADD_TEST(suite, test_nhdp_serialization);
    SUITE_ADD_TEST(suite, test_sm3hdp_cell);
    SUITE_ADD_TEST(suite, test_sm3Hdp_dpDiagonal);
    SUITE_ADD_TEST(suite, test_sm3Hdp_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_sm3Hdp_getAlignedPairsWithBanding);
    //SUITE_ADD_TEST(suite, test_hdpHmm);
    //SUITE_ADD_TEST(suite, test_continuousPairHDPHmm_em);
    return suite;
}
