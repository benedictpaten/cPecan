//
//  nanopore_hdp.c
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "nanopore_hdp.h"
#include "sonLib.h"

// in 0-based index
//#declare KMER_COL 9
//#declare SIGNAL_COL 13
//#declare NUM_COLS 14

void normal_inverse_gamma_params_from_minION(const char* signal_lookup_table, double* mu_out, double* nu_out,
                                             double* alpha_out, double* beta_out) {
    
    FILE* lookup_table = fopen(signal_lookup_table, "r");
    
    double* means;
    int64_t length;
    
    //TODO: get the mean vector
    
    normal_inverse_gamma_params(means, length, mu_out, nu_out, alpha_out, beta_out);
    
}

// fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int num_dps, int depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int sampling_grid_length,
                                         const char* signal_lookup_table_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(signal_lookup_table_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc(num_dps, depth, gamma, sampling_grid_start, sampling_grid_stop,
                             sampling_grid_length, mu, nu, alpha, beta);
}

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int num_dps, int depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int sampling_grid_length,
                                           const char* signal_lookup_table_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(signal_lookup_table_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc_2(num_dps, depth, gamma_alpha, gamma_beta, sampling_grid_start,
                               sampling_grid_stop, sampling_grid_length, mu, nu, alpha, beta);
}

int64_t *stList_toIntPtr(stList* list, int64_t* length_out) {
    int64_t length = (int64_t) stList_length(list);
    int64_t* int_arr = (int64_t*) malloc(sizeof(int64_t) * length);
    int64_t* entry;
    for (int64_t i = 0; i < length; i++) {
        entry = (int64_t*) stList_get(list, i);
        int_arr[0] = *entry;
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
        double_arr[0] = *entry;
    }
    *length_out = length;
    return double_arr;
}
/*
void update_hdp_from_alignment(HierarchicalDirichletProcess* hdp, const char* alignment_filepath,
                               int (*kmer_to_dp_id_func) (char*), bool has_header) {
    
    stList* signal_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);
    
    FILE* align_file = fopen(alignment_filepath, "r");
    
    stList* tokens;
    int64_t line_length;
    char* kmer;
    char* signal_str;
    int* dp_id_ptr;
    double* signal_ptr;
    bool warned = false;
    
    char* line = stFile_getLineFromFile(align_file);
    if (has_header) {
        line = stFile_getLineFromFile(align_file);
    }
    while (line != NULL) {
        tokens = stString_split(line);
        line_length = stList_getLength(tokens);
        
        if (!warned) {
            if (line_length != NUM_COLS) {
                fprintf(stderr, "Input format has changed from design period, HDP may receive incorrect data.");
                warned = true;
            }
        }
        
        signal_str = (char*) stList_get(tokens, SIGNAL_COL);
        kmer = (char*) stList_get(tokens, KMER_COL);
        
        signal_ptr = (double*) malloc(sizeof(double));
        dp_id_ptr = (double*) malloc(sizeof(int));
        
        sscanf(signal_str, "%lf", signal_ptr);
        *dp_id_ptr = kmer_to_dp_id_func(kmer);
        
        stList_append(signal_list, signal_ptr);
        stList_append(dp_id_list, dp_id_ptr);
        
        stList_destruct(tokens);
        free(line);
        line = stFile_getLineFromFile(align_file);
    }
    
    fclose(align_file);
    
    int data_length;
    
    double *signal = stList_toDoublePtr(signal_list, &data_length);
    int64_t *dp_ids = stList_toIntPtr(dp_id_list, &data_length);
    
    stList_destruct(signal_list);
    stList_destruct(dp_id_list);
    
    reset_hdp_data(hdp);
    pass_data_to_hdp(hdp, signal, dp_ids, data_length);
}
*/






