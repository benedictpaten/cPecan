#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <nanopore.h>
#include "stateMachine.h"
#include "CuTest.h"
#include "hdp.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "continuousHmm.h"
#include "discreteHmm.h"
#include "emissionMatrix.h"
#include "multipleAligner.h"
#include "randomSequences.h"


int64_t* stList_toIntPtr(stList* list, int64_t* length_out) {
    int64_t length = (int64_t) stList_length(list);
    int64_t* int_arr = (int64_t*) malloc(sizeof(int64_t) * length);
    int64_t* entry;
    for (int64_t i = 0; i < length; i++) {
        entry = (int64_t*) stList_get(list, i);
        int_arr[i] = *entry;
    }
    *length_out = length;
    return int_arr;
}

double* stList_toDoublePtr(stList* list, int64_t* length_out) {
    int64_t length  = stList_length(list);
    double* double_arr = (double*) malloc(sizeof(double) * length);
    double* entry;
    for (int64_t i = 0; i < length; i++) {
        entry = (double*) stList_get(list, i);
        double_arr[i] = *entry;
    }
    *length_out = length;
    return double_arr;
}

void load_data(const char* data_filepath, const char* dp_id_filepath,
               double** data_out, int64_t **dp_ids_out, int64_t *length_out) {

    FILE* data_file = fopen(data_filepath, "r");
    FILE* dp_id_file = fopen(dp_id_filepath, "r");

    stList* data_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);

    char* line = stFile_getLineFromFile(data_file);

    double* datum_ptr;
    while (line != NULL) {

        datum_ptr = (double*) malloc(sizeof(double));

        sscanf(line, "%lf", datum_ptr);

        stList_append(data_list, datum_ptr);

        free(line);
        line = stFile_getLineFromFile(data_file);
    }


    line = stFile_getLineFromFile(dp_id_file);
    int64_t *dp_id_ptr;
    while (line != NULL) {
        dp_id_ptr = (int64_t *) malloc(sizeof(int));

        sscanf(line, "%lld", dp_id_ptr);

        stList_append(dp_id_list, dp_id_ptr);

        free(line);
        line = stFile_getLineFromFile(dp_id_file);
    }

    int64_t data_length;
    int64_t dp_ids_length;

    double* data = stList_toDoublePtr(data_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);

    if (data_length != dp_ids_length) {
        fprintf(stderr, "Data and DP ID files have different lengths\n");
        exit(EXIT_FAILURE);
    }

    *data_out = data;
    *dp_ids_out = dp_ids;
    *length_out = data_length;

    stList_destruct(dp_id_list);
    stList_destruct(data_list);

    fclose(data_file);
    fclose(dp_id_file);
}

static void test_hdpBasic(CuTest *testCase) {
    /*
     * This vignette walks through creating a hierarchical Dirichlet process
     * with this structure:
     *           _
     *          |H|    base normal-inverse gamma distribution (not a DP)
     *          |
     *          |_            _______
     *          |0| <--------|gamma_0|
     *        __|__
     *      _|     |_         _______
     *     |1|     |2| <-----|gamma_1|
     *  ___|___     _|_
     * |_  |_  |_  |_  |_     _______
     * |3| |4| |5| |6| |7| <-|gamma_2|
     */
    srand(30);
    // the Dirichlet processes are the numbered boxes
    int64_t num_dir_proc = 8;
    // the depth of the tree
    int64_t depth = 3;

    // parameters of the normal inverse gamma base distribution
    double mu = 0.0;
    double nu = 10.0;
    double alpha = 10.0; // note: this parameter must be integer or half-integer valued
    double beta = 10.0;

    // the grid along which the HDP will record distribution samples
    int64_t grid_length = 250;
    double grid_start = -10.0;
    double grid_end = 10.0;

    // initialize an HDP
    fprintf(stderr, "Initializing HDP...\n");

    // choose whether to pre-define concentration parameters (gamma and alpha_0 in Teh et al)
    // or sample them from a Gamma distribution

    // the concentration parameters at each depth
    //double* gamma = (double*) malloc(sizeof(double) * depth);
    //gamma[0] = 5.0; gamma[1] = 5.0; gamma[2] = 20.0;
    //HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dir_proc, depth, gamma,
    //                                                      grid_start, grid_end, grid_length,
    //                                                      mu, nu, alpha, beta);


    // parameters for distributions of concentration parameters at each depth (gamma and alpha_0 in Teh et al)
    double *gamma_alpha = (double*) malloc(sizeof(double) * depth);
    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
    double *gamma_beta = (double*) malloc(sizeof(double) * depth);
    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
                                                            gamma_beta, grid_start, grid_end,
                                                            grid_length, mu, nu, alpha, beta);

    // establish the topology of the tree
    fprintf(stderr, "Establishing HDP tree topology...\n");
    set_dir_proc_parent(hdp, 1, 0);
    set_dir_proc_parent(hdp, 2, 0);
    set_dir_proc_parent(hdp, 3, 1);
    set_dir_proc_parent(hdp, 4, 1);
    set_dir_proc_parent(hdp, 5, 1);
    set_dir_proc_parent(hdp, 6, 2);
    set_dir_proc_parent(hdp, 7, 2);
    // note: the result must be a perfectly balanced tree (same depth at every leaf)
    // or you will get an error here
    finalize_hdp_structure(hdp);

    fprintf(stderr, "Loading data from disk...\n");
    double *data;
    int64_t *data_pt_dps;
    int64_t  data_length;
    load_data("../../hdp_mixture/test/data.txt",
              "../../hdp_mixture/test/dps.txt",
              &data, &data_pt_dps, &data_length);

    int64_t new_data_length = data_length / 2;
    double *new_data = (double*) malloc(sizeof(double) * new_data_length);
    int64_t *new_data_pt_dps = (int64_t*) malloc(sizeof(int64_t) * new_data_length);

    for (int64_t i = 0; i < new_data_length; i++) {
        new_data[i] = data[i] + 1.0;
        new_data_pt_dps[i] = data_pt_dps[i];
    }

    fprintf(stderr, "Giving HDP data...\n");
    pass_data_to_hdp(hdp, data, data_pt_dps, data_length);
    // note: you can also pass data before finalizing the structure
    // note: it is not necessary to observe every Dirichlet process in the data
    int64_t num_samples = 1000;
    int64_t burn_in = 10000;
    int64_t thinning = 50;

    // sample from the posterior distribution of distributions
    fprintf(stderr, "Executing Gibbs sampling...\n");
    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning);

    // calculate the mean a posteriori estimate of each distribution
    fprintf(stderr, "Computing mean a posteriori distributions...\n");
    finalize_distributions(hdp);

    // query with new values
    double x = 2.1;
    int dp_id = 4;
    double density = dir_proc_density(hdp, x, dp_id);
    st_uglyf("Got density %f\n", density);

    // reset the HDP without needing to re-initialize it and provide new data
    fprintf(stderr, "Reseting HDP data...\n");
    reset_hdp_data(hdp);

    fprintf(stderr, "Giving HDP new data...\n");
    pass_data_to_hdp(hdp, new_data, new_data_pt_dps, new_data_length);
    fprintf(stderr, "Executing Gibbs sampling...\n");
    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning);
    fprintf(stderr, "Computing mean a posteriori distributions...\n");
    finalize_distributions(hdp);

    // query density values with the new distributions
    density = dir_proc_density(hdp, 3.4, 6);
    st_uglyf("Got density %f\n", density);

    // free the memory
    fprintf(stderr, "Destroying HDP...\n");
    destroy_hier_dir_proc(hdp);

    printf("Completed!\n");
}


CuSuite *hdpTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    //SUITE_ADD_TEST(suite, test_hdpBasic);
    return suite;
}
