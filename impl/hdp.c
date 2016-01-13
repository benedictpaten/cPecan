#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "sonLib.h"
#include "ranlib.h"

#define N_IG_NUM_PARAMS 4

typedef struct Factor Factor;
typedef struct DirichletProcess DirichletProcess;
typedef enum FactorType FactorType;

enum FactorType {
    BASE,
    MIDDLE,
    DATA_PT
};

struct Factor {
    FactorType factor_type;
    struct Factor* parent;
    stSet* children;
    double* factor_data;
    struct DirichletProcess* dp;
};

struct DirichletProcess {
    int64_t id;
    struct HierarchicalDirichletProcess* hdp;
    double* gamma;
    int64_t depth;

    struct DirichletProcess* parent;
    stList* children;
    stSet* factors;
    int64_t num_factor_children;

    double base_factor_wt;
    double* posterior_predictive;
    double* spline_slopes;

    double cached_factor_mean;
    double cached_factor_sum_sq_dev;
    int64_t cached_factor_size;

    bool observed;
};

struct HierarchicalDirichletProcess {
    bool finalized;
    double* data;
    int64_t *data_pt_dp_id;
    int64_t data_length;

    struct DirichletProcess* base_dp;
    struct DirichletProcess** dps;
    int64_t num_dps;

    // normal-inverse gamma parameters
    double mu;
    double nu;
    double two_alpha;
    double beta;

    double* sampling_grid;
    int64_t grid_length;

    int64_t samples_taken;
    bool splines_finalized;

    struct SumOfLogsMemo* log_sum_memo;
    
    int64_t depth;
    bool sample_gamma;
    double* gamma;
    
    double* gamma_alpha;
    double* gamma_beta;
    double* w_aux_vector;
    bool* s_aux_vector;
};

void cache_base_factor_params(Factor* fctr, double mu, double nu, double two_alpha, double beta, double log_post_term) {
    if (fctr->factor_type != BASE) {
        fprintf(stderr, "Can only cache parameters for base factors.\n");
        exit(EXIT_FAILURE);
    }

    double* param_array = fctr->factor_data;
    param_array[0] = mu;
    param_array[1] = nu;
    param_array[2] = two_alpha;
    param_array[3] = beta;
    param_array[4] = log_post_term;
}

Factor* new_base_factor(HierarchicalDirichletProcess* hdp) {
    //printf("$$$ making a new base factor\n");
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = BASE;

    fctr->factor_data = (double*) malloc(sizeof(double) * (N_IG_NUM_PARAMS + 1));
    cache_base_factor_params(fctr, hdp->mu, hdp->nu, hdp->two_alpha, hdp->beta, 1.0);

    fctr->parent = NULL;
    fctr->children = stSet_construct();

    DirichletProcess* base_dp = hdp->base_dp;
    fctr->dp = base_dp;
    stSet_insert(base_dp->factors, (void*) fctr);

    return fctr;
}

Factor* new_middle_factor(DirichletProcess* dp) {
    if (dp->parent == NULL) {
        fprintf(stderr, "Attempted to create middle factor in root Dirichlet process.\n");
        exit(EXIT_FAILURE);
    }

    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = MIDDLE;
    fctr->factor_data = NULL;

    // note: assigning to parent handled externally
    fctr->parent = NULL;
    fctr->children = stSet_construct();

    fctr->dp = dp;
    stSet_insert(dp->factors, (void*) fctr);
    return fctr;
}

Factor* new_data_pt_factor(HierarchicalDirichletProcess* hdp, int64_t data_pt_idx) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = DATA_PT;
    fctr->factor_data = &(hdp->data[data_pt_idx]);

    // note: assigning to parent handled externally
    fctr->parent = NULL;
    fctr->children = NULL;

    fctr->dp = NULL;

    return fctr;
}

void destroy_factor(Factor* fctr) {
    stSet* children = fctr->children;
    if (children != NULL) {
        if (stSet_size(children) > 0) {
            //fprintf(stderr, "Attempted to destroy factor that still has children.\n");
            exit(EXIT_FAILURE);
        }
        stSet_destruct(children);
    }

    Factor* parent = fctr->parent;
    if (parent != NULL) {
        stSet_remove(parent->children, (void*) fctr);
        (parent->dp->num_factor_children)--;
        if (stSet_size(parent->children) == 0) {
            destroy_factor(parent);
        }
    }

    if (fctr->factor_type == BASE) {
        //printf("$$$ destroying a base factor\n");
        free(fctr->factor_data);
    }
    
    DirichletProcess* dp = fctr->dp;
    if (dp != NULL) {
        stSet_remove(dp->factors, (void*) fctr);
    }

    free(fctr);
}

Factor* get_base_factor(Factor* fctr) {
    while (fctr->factor_type != BASE) {
        fctr = fctr->parent;
        if (fctr == NULL) {
            break;
        }
    }
    return fctr;
}

double get_factor_data_pt(Factor* fctr) {
    if (fctr->factor_type != DATA_PT) {
        //fprintf(stderr, "Attempted to access data point from non-leaf factor.\n");
        exit(EXIT_FAILURE);
    }
    return *(fctr->factor_data);
}

void get_factor_sum_internal(Factor* fctr, double* sum, int64_t* num_data) {
    if (fctr->factor_type == DATA_PT) {
        //printf("attempts to get data from data point factor\n");
        *sum += get_factor_data_pt(fctr);
        // TODO: there should be a way to use the parent's counters instead of recounting the data pts
        (*num_data)++;
    }
    else {
        //printf("attempts to get iterator for children\n");
        stSetIterator* child_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            //printf("recurses down from depth %d\n", fctr->dp->depth);
            get_factor_sum_internal(child_fctr, sum, num_data);
            //printf("exits recursion back to depth %d\n", fctr->dp->depth);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        //printf("attempts to destruct iterator\n");
        stSet_destructIterator(child_iter);
    }
}

void get_factor_ssd_internal(Factor* fctr, double center, double* sum_sq_dev) {
    if (fctr->factor_type == DATA_PT) {
        double dev = get_factor_data_pt(fctr) - center;
        *sum_sq_dev += dev * dev;
    }
    else {
        stSetIterator* child_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            get_factor_ssd_internal(child_fctr, center, sum_sq_dev);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        stSet_destructIterator(child_iter);
    }
}

void get_factor_stats(Factor* fctr, double* mean_out, double* sum_sq_dev_out, int64_t* num_data_out) {
    *mean_out = 0.0;
    *sum_sq_dev_out = 0.0;
    *num_data_out = 0;
    //printf("attempts to get sum and count\n");
    get_factor_sum_internal(fctr, mean_out, num_data_out);
    *mean_out /= (double) *num_data_out;
    //printf("attempts to get sum of sq devs\n");
    get_factor_ssd_internal(fctr, *mean_out, sum_sq_dev_out);
}

void add_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    double* param_array = fctr->factor_data;

    double mu_prev = param_array[0];
    double nu_prev = param_array[1];
    double two_alpha_prev = param_array[2];
    double beta_prev = param_array[3];

    double nu_post = nu_prev + num_data;
    double mu_post = (mu_prev * nu_prev + mean * num_data) / nu_post;
    double two_alpha_post = two_alpha_prev + num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_post = beta_prev + .5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_post, two_alpha_post, beta_post,
                                                          fctr->dp->hdp->log_sum_memo);

    cache_base_factor_params(fctr, mu_post, nu_post, two_alpha_post, beta_post, log_post_term);
}

void remove_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    double* param_array = fctr->factor_data;

    double mu_post = param_array[0];
    double nu_post = param_array[1];
    double two_alpha_post = param_array[2];
    double beta_post = param_array[3];

    double nu_prev = nu_post - num_data;
    double mu_prev = (mu_post * nu_post - mean * num_data) / nu_prev;
    double two_alpha_prev = two_alpha_post - num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_prev = beta_post - 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_prev, two_alpha_prev, beta_prev,
                                                          fctr->dp->hdp->log_sum_memo);

    cache_base_factor_params(fctr, mu_prev, nu_prev, two_alpha_prev, beta_prev, log_post_term);
}

double factor_parent_joint_log_likelihood(Factor* fctr, Factor* parent) {
    Factor* base_fctr = get_base_factor(parent);
    DirichletProcess* dp = fctr->dp;

    double num_reassign = (double) dp->cached_factor_size;
    double mean_reassign = dp->cached_factor_mean;
    double sum_sq_devs = dp->cached_factor_sum_sq_dev;

    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double nu_numer = nu_denom + num_reassign;

    double mean_dev = mean_reassign - mu_denom;
    double sq_mean_dev = nu_denom * num_reassign * mean_dev * mean_dev / nu_numer;

    double two_alpha_numer = two_alpha_denom + num_reassign;
    double beta_numer = beta_denom + 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_denom = param_array[4];
    double log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                      dp->hdp->log_sum_memo);

    return -0.5 * num_reassign * log(2.0 * M_PI) + log_numer - log_denom;
}

double data_pt_factor_parent_likelihood(Factor* data_pt_fctr, Factor* parent) {
    if (data_pt_fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Can only access data point likelihood for data point factors.\n");
        exit(EXIT_FAILURE);
    }

    double data_pt = get_factor_data_pt(data_pt_fctr);
    Factor* base_fctr = get_base_factor(parent);
    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double nu_numer = nu_denom + 1.0;

    double mean_dev = data_pt - mu_denom;
    double sq_mean_dev = nu_denom * mean_dev * mean_dev / nu_numer;

    double two_alpha_numer = two_alpha_denom + 1.0;
    double beta_numer = beta_denom + 0.5 * sq_mean_dev;

    double log_denom = param_array[4];
    double log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                      base_fctr->dp->hdp->log_sum_memo);

    return (1.0 / sqrt(2.0 * M_PI)) * exp(log_numer - log_denom);
}

void evaluate_posterior_predictive(Factor* base_fctr, double* x, double* pdf_out, int64_t length,
							       SumOfLogsMemo* log_sum_memo) {
    if (base_fctr->factor_type != BASE) {
        fprintf(stderr, "Can only evaluate posterior predictive of base factors.\n");
        exit(EXIT_FAILURE);
    }

    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double log_denom = param_array[4];

    double nu_numer = nu_denom + 1.0;
    double two_alpha_numer = two_alpha_denom + 1.0;
    double nu_ratio = nu_denom / nu_numer;
    double pi_factor = (1.0 / sqrt(2.0 * M_PI));

    double mean_dev;
    double sq_mean_dev;
    double beta_numer;
    double log_numer;
    for (int64_t i = 0; i < length; i++) {
        mean_dev = x[i] - mu_denom;
        sq_mean_dev = nu_ratio * mean_dev * mean_dev;
        beta_numer = beta_denom + 0.5 * sq_mean_dev;

        log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                   log_sum_memo);

        pdf_out[i] = pi_factor * exp(log_numer - log_denom);
    }
}

void evaluate_prior_predictive(HierarchicalDirichletProcess* hdp,
                               double* x, double* pdf_out, int64_t length) {
    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    int64_t two_alpha = (int64_t) hdp->two_alpha;
    double beta = hdp->beta;

    double nu_factor = nu / (2.0 * (nu + 1.0) * beta);
    double alpha_term = exp(log_gamma_half(two_alpha + 1, hdp->log_sum_memo)
                            - log_gamma_half(two_alpha, hdp->log_sum_memo));
    double beta_term = sqrt(nu_factor / M_PI);
    double constant_term = alpha_term * beta_term;
    double power = -0.5 * (two_alpha + 1.0);

    for (int64_t i = 0; i < length; i++) {
        double dev = x[i] - mu;
        double var_term = pow(1.0 + nu_factor * dev * dev, power);

        pdf_out[i] = constant_term * var_term;
    }
}

double prior_likelihood(HierarchicalDirichletProcess* hdp, Factor* fctr) {
    // TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    double dbl_two_alpha = hdp->two_alpha;
    int64_t two_alpha = (int64_t) dbl_two_alpha;
    double beta = hdp->beta;

    if (fctr->factor_type == DATA_PT) {
        double data_pt = get_factor_data_pt(fctr);
        double dev = data_pt - mu;

        double alpha_term = exp(log_gamma_half(two_alpha + 1, hdp->log_sum_memo)
                                - log_gamma_half(two_alpha, hdp->log_sum_memo));

        double nu_term = nu / (2.0 * (nu + 1.0) * beta);
        double beta_term = pow(1.0 + nu_term * dev * dev, -0.5 * (dbl_two_alpha + 1.0));

        return alpha_term * sqrt(nu_term / M_PI) * beta_term;
    }
    else {
        DirichletProcess* dp = fctr->dp;
        int64_t num_reassign = dp->cached_factor_size;
        double dbl_reassign = (double) num_reassign;
        double mean_reassign = dp->cached_factor_mean;
        double sum_sq_devs = dp->cached_factor_sum_sq_dev;

        double mean_dev = mean_reassign - mu;
        double sq_mean_dev = nu * dbl_reassign * mean_dev * mean_dev / (nu + dbl_reassign);

        double log_alpha_term = log_gamma_half(two_alpha + num_reassign, hdp->log_sum_memo)
                                - log_gamma_half(two_alpha, hdp->log_sum_memo);
        double log_nu_term = 0.5 * log(nu / (nu + dbl_reassign));
        double log_beta_term_1 = -0.5 * dbl_reassign * log(2.0 * M_PI * beta);
        double log_beta_term_2 = -0.5 * (dbl_two_alpha + dbl_reassign)
                                 * log(1.0 + (0.5 / beta) * (sum_sq_devs + sq_mean_dev));
        return exp(log_alpha_term + log_nu_term + log_beta_term_1 + log_beta_term_2);
    }
}

double unobserved_factor_likelihood(Factor* fctr, DirichletProcess* dp) {
    //printf("attempts to get parent dp\n");
    DirichletProcess* parent_dp = dp->parent;
    if (parent_dp == NULL) {
        //printf("no parent, attempts to compute prior\n");
        return prior_likelihood(dp->hdp, fctr);
    }
    else {
        //printf("attempts to get parent gamma\n");
        double parent_gamma = *(parent_dp->gamma);
        double likelihood = 0.0;
        
        stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
        Factor* parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
        double fctr_size;
        if (fctr->factor_type == DATA_PT) {
            //printf("attempts to sum up likelihood for data pt factor\n");
            while (parent_fctr != NULL) {
                fctr_size = (double) stSet_size(parent_fctr->children);
                likelihood += fctr_size * data_pt_factor_parent_likelihood(fctr, parent_fctr);

                parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
            }
        }
        else {
            //printf("attempts to sum up likelihood for middle factor\n");
            while (parent_fctr != NULL) {
            	fctr_size = (double) stSet_size(parent_fctr->children);
                likelihood += fctr_size * exp(factor_parent_joint_log_likelihood(fctr, parent_fctr));

                parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
            }
        }
        stSet_destructIterator(parent_fctr_iter);
        
        //printf("recurses down\n");
        likelihood += parent_gamma * unobserved_factor_likelihood(fctr, parent_dp);
        
        //printf("exits recursion, attempts to normalize likelihood\n");
        likelihood /= (parent_gamma + (double) parent_dp->num_factor_children);

        return likelihood;
    }
}

DirichletProcess* new_dir_proc() {
    DirichletProcess* dp = (DirichletProcess*) malloc(sizeof(DirichletProcess));

    dp->gamma = NULL;
    dp->depth = -1;
    dp->parent = NULL;
    dp->children = stList_construct();
    dp->factors = stSet_construct();
    dp->num_factor_children = 0;

    dp->cached_factor_mean = 0.0;
    dp->cached_factor_sum_sq_dev = 0.0;
    dp->cached_factor_size = 0;

    dp->base_factor_wt = 0.0;
    dp->posterior_predictive = NULL;
    dp->spline_slopes = NULL;

    dp->observed = false;
    return dp;
}

void clear_factor_tree(Factor* fctr) {
    //printf("enters clear factor tree function\n");
    if (fctr->children != NULL) {
        stSetIterator* child_fctr_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            //printf("recurses down\n");
            clear_factor_tree(child_fctr);
            //printf("exits recursion\n");
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
    }
    else {
        destroy_factor(fctr);
    }
    
    //printf("attempts to destroy factor\n");
    
}

void destroy_dir_proc_factor_tree(DirichletProcess* dp) {
    if (stSet_size(dp->factors) == 0) {
        return;
    }
    //printf("enters destroy DP tree function\n");
    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    while (fctr != NULL) {
        clear_factor_tree(fctr);
        fctr = (Factor*) stSet_getNext(fctr_iter);
    }
    stSet_destructIterator(fctr_iter);
}

void destroy_dir_proc(DirichletProcess* dp) {
    destroy_dir_proc_factor_tree(dp);
    stSet_destruct(dp->factors);
    
    if (dp->children != NULL) {
        stListIterator* st_iterator = stList_getIterator(dp->children);
        DirichletProcess* dp_child = (DirichletProcess*) stList_getNext(st_iterator);
        while (dp_child != NULL) {
            destroy_dir_proc(dp_child);
            dp_child = (DirichletProcess*) stList_getNext(st_iterator);
        }
        stList_destructIterator(st_iterator);
        
        stList_destruct(dp->children);
    }

    if (dp->parent != NULL) {
        stList_removeItem(dp->parent->children, (void*) dp);
    }

    free(dp->posterior_predictive);
    free(dp->spline_slopes);
    free(dp);
}

// fixed concentration parameters
HierarchicalDirichletProcess* new_hier_dir_proc(int64_t num_dps, int64_t depth,
                                                double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop,
                                                int64_t sampling_grid_length,
                                                double mu, double nu, double alpha,
                                                double beta) {
    if (nu <= 0.0) {
        fprintf(stderr, "nu parameter of Normal-Inverse Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    //if (alpha <= 0.0) {
    //    fprintf(stderr, "alpha parameter of Normal-Inverse Gamma distribution must be positive.\n");
    //    exit(EXIT_FAILURE);
    //}
    if (beta <= 0.0) {
        fprintf(stderr, "beta parameter of Normal-Inverse Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (2 * alpha != (int) 2 * alpha || alpha <= 1.0) {
            fprintf(stderr, "Normal-Inverse Gamma parameter 'alpha' must be integer or half-integer > 1.0.\n");
            exit(EXIT_FAILURE);
        }
    if (gamma != NULL) {
        for (int i = 0; i < depth; i++) {
            if (gamma[i] <= 0) {
                fprintf(stderr, "Concentration parameter gamma must be postive.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    double *grid = linspace(sampling_grid_start, sampling_grid_stop, sampling_grid_length);

    HierarchicalDirichletProcess* hdp = (HierarchicalDirichletProcess*) malloc(sizeof(HierarchicalDirichletProcess));

    // normal-inverse gamma parameters

    hdp->mu = mu;
    hdp->nu = nu;
    hdp->two_alpha = 2.0 * alpha;
    hdp->beta = beta;

    hdp->gamma = gamma;
    hdp->depth = depth;

    hdp->finalized = false;

    hdp->num_dps = num_dps;
    DirichletProcess** dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
        DirichletProcess* dp = new_dir_proc();
        dp->id = i;
        dp->hdp = hdp;
        dps[i] = dp;
    }

    hdp->dps = dps;
    hdp->base_dp = NULL;

    hdp->sampling_grid = grid;
    hdp->grid_length = sampling_grid_length;
    hdp->samples_taken = 0;
    hdp->splines_finalized = false;

    hdp->data = NULL;
    hdp->data_pt_dp_id = NULL;
    hdp->data_length = -1;

    hdp->log_sum_memo = new_log_sum_memo();

    hdp->sample_gamma = false;
    hdp->s_aux_vector = NULL;
    hdp->w_aux_vector = NULL;

    return hdp;
}

// Gamma prior on concentration parameters
HierarchicalDirichletProcess *new_hier_dir_proc_2(int64_t num_dps, int64_t depth,
                                                  double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start,
                                                  double sampling_grid_stop,
                                                  int64_t sampling_grid_length,
                                                  double mu, double nu,
                                                  double alpha, double beta) {

    for (int64_t i = 0; i < depth; i++) {
        if (gamma_alpha[i] <= 0.0) {
            st_errAbort("alpha parameter of Gamma distribution must be positive.\n");
            //fprintf(stderr, "alpha parameter of Gamma distribution must be positive.\n");
            //exit(EXIT_FAILURE);
        }
        if (gamma_beta[i] <= 0.0) {
            st_errAbort("beta parameter of Gamma distribution must be positive.\n");
            //fprintf(stderr, "beta parameter of Gamma distribution must be positive.\n");
            //exit(EXIT_FAILURE);
        }
    }
    
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dps, depth, NULL,
                                                          sampling_grid_start, sampling_grid_stop,
                                                          sampling_grid_length,
                                                          mu, nu, alpha,
                                                          beta);

    hdp->sample_gamma = true;
    
    hdp->gamma_alpha = gamma_alpha;
    hdp->gamma_beta = gamma_beta;

    double* w = (double*) malloc(sizeof(double) * num_dps);
    hdp->w_aux_vector = w;
    bool* s = (bool*) malloc(sizeof(bool) * num_dps);
    hdp->s_aux_vector = s;

    for (int64_t i = 0; i < num_dps; i++) {
        w[i] = 1.0;
        s[i] = false;
    }
    
    // init to prior expected value
    double* gamma = (double*) malloc(sizeof(double) * depth);
    hdp->gamma = gamma;
    for (int64_t i = 0; i < depth; i++) {
        gamma[i] = gamma_alpha[i] / gamma_beta[i];
    }
    
    
    return hdp;
}

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp) {
    //printf("1\n");
    destroy_dir_proc(hdp->base_dp);
    //printf("2\n");
    free(hdp->gamma);
    //printf("3\n");
    free(hdp->data);
    //printf("4\n");
    free(hdp->data_pt_dp_id);
    //printf("5\n");
    free(hdp->dps);
    //printf("6\n");
    free(hdp->sampling_grid);
    //printf("7\n");
    destroy_log_sum_memo(hdp->log_sum_memo);
    //printf("8\n");
    free(hdp->gamma_alpha);
    //printf("9\n");
    free(hdp->gamma_beta);
    //printf("10\n");
    free(hdp->w_aux_vector);
    //printf("11\n");
    free(hdp->s_aux_vector);
    //printf("12\n");
    free(hdp);
}

void establish_base_dp(HierarchicalDirichletProcess* hdp) {
    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;

    DirichletProcess* dp;
    for (int64_t i = 0; i < num_dps; i++) {
        dp = dps[i];
        if (dp->parent == NULL) {
            if (hdp->base_dp == NULL) {
                hdp->base_dp = dp;
            }
            else {
                fprintf(stderr, "Hierarchical Dirichlet process contains orphaned Dirichlet process.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (hdp->base_dp == NULL) {
        fprintf(stderr, "Hierarchical Dirichlet process does not contain base Dirichlet process.\n");
        exit(EXIT_FAILURE);
    }
}

// DFS to verify that Dirichlet processes follow tree structure
void verify_dp_tree(HierarchicalDirichletProcess* hdp) {
    int64_t num_dps = hdp->num_dps;
    bool* visited = (bool*) malloc(sizeof(bool) * num_dps);
    for (int i = 0; i < num_dps; i++) {
        visited[i] = false;
    }

    DirichletProcess* base_dp = hdp->base_dp;
    stList* stck = stList_construct();
    stList_append(stck, (void*) base_dp);

    DirichletProcess* dp;
    while (stList_length(stck) > 0) {
        dp = (DirichletProcess*) stList_pop(stck);
        if (visited[dp->id]) {
            fprintf(stderr, "Hierarchical Dirichlet process does not have tree structure.\n");
            exit(EXIT_FAILURE);
        }
        visited[dp->id] = true;

        stListIterator* child_iter = stList_getIterator(dp->children);
        DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
        while (child != NULL) {
            stList_append(stck, (void*) child);
            child = (DirichletProcess*) stList_getNext(child_iter);
        }
        stList_destructIterator(child_iter);
    }
    stList_destruct(stck);
    free(visited);
}

void verify_tree_depth(HierarchicalDirichletProcess* hdp, DirichletProcess* dp, int64_t current_depth,
                       int64_t leaf_depth) {
    dp->gamma = &(hdp->gamma[current_depth]);
    dp->depth = current_depth;

    if (stList_length(dp->children) == 0) {
        if (current_depth != leaf_depth) {
            fprintf(stderr, "Hierarchical Dirichlet process has leaf Dirichlet process at incorrect depth.\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        stListIterator* st_iterator = stList_getIterator(dp->children);
        DirichletProcess* child = (DirichletProcess*) stList_getNext(st_iterator);
        while (child != NULL) {
            verify_tree_depth(hdp, child, current_depth + 1, leaf_depth);
            child = (DirichletProcess*) stList_getNext(st_iterator);
        }
        stList_destructIterator(st_iterator);
    }
}

void verify_valid_dp_assignments(HierarchicalDirichletProcess* hdp) {

    int64_t length = hdp->data_length;
    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    int64_t* dp_ids = hdp->data_pt_dp_id;

    int64_t id;
    DirichletProcess* dp;
    for (int64_t i = 0; i < length; i++) {
        id = dp_ids[i];
        if (id >= num_dps || id < 0) {
            fprintf(stderr, "Data point is assigned to non-existent Dirichlet process.\n");
            exit(EXIT_FAILURE);
        }

        dp = dps[id];
        if (stList_length(dp->children) > 0) {
            fprintf(stderr, "Data point cannot be assigned to non-leaf Dirichlet process.\n");
            exit(EXIT_FAILURE);
        }
    }
}

void mark_observed_dps(HierarchicalDirichletProcess* hdp) {
    // mark newly observed dps
    int64_t length = hdp->data_length;
    DirichletProcess** dps = hdp->dps;
    int64_t* dp_ids = hdp->data_pt_dp_id;
    int64_t grid_length = hdp->grid_length;

    DirichletProcess* dp;
    double* pdf;
    int64_t id;

    for (int i = 0; i < length; i++) {
        id = dp_ids[i];
        dp = dps[id];
        while (dp != NULL) {
            if (dp->observed) {
                break;
            }
            dp->observed = true;

            pdf = (double*) malloc(sizeof(double) * grid_length);
            dp->posterior_predictive = pdf;
            for (int j = 0; j < grid_length; j++) {
                pdf[j] = 0.0;
            }

            dp = dp->parent;
        }
    }
}

void init_factors_internal(DirichletProcess* dp, Factor* parent_fctr, stList** data_pt_fctr_lists) {
    if (!dp->observed) {
        return;
    }
    Factor* fctr = new_middle_factor(dp);
    fctr->parent = parent_fctr;
    stSet_insert(parent_fctr->children, (void*) fctr);
    
    if (stList_length(dp->children) == 0) {
        stSet* children = fctr->children;
        stListIterator* data_pt_fctr_iter = stList_getIterator(data_pt_fctr_lists[dp->id]);
        Factor* data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        while (data_pt_fctr != NULL) {
            data_pt_fctr->parent = fctr;
            stSet_insert(children, (void*) data_pt_fctr);
            data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        }
        stList_destructIterator(data_pt_fctr_iter);
    }
    else {
        stListIterator* child_dp_iter = stList_getIterator(dp->children);
        DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        while (child_dp != NULL) {
            init_factors_internal(child_dp, fctr, data_pt_fctr_lists);
            child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        }
        stList_destructIterator(child_dp_iter);
    }
}

void init_factors(HierarchicalDirichletProcess* hdp) {
    int64_t data_length = hdp->data_length;
    int64_t* data_pt_dp_id = hdp->data_pt_dp_id;
    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;
    
    stList** data_pt_fctr_lists = (stList**) malloc(sizeof(stList*) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
        data_pt_fctr_lists[i] = NULL;
    }
    
    Factor* data_pt_fctr;
    int64_t dp_id;
    stList* fctr_list;
    for (int64_t data_pt_idx = 0; data_pt_idx < data_length; data_pt_idx++) {
        dp_id = data_pt_dp_id[data_pt_idx];
        fctr_list = data_pt_fctr_lists[dp_id];
        if (fctr_list == NULL) {
            fctr_list = stList_construct();
            data_pt_fctr_lists[dp_id] = fctr_list;
        }
        data_pt_fctr = new_data_pt_factor(hdp, data_pt_idx);
        stList_append(fctr_list, (void*) data_pt_fctr);
    }
    
    DirichletProcess* base_dp = hdp->base_dp;
    Factor* root_factor = new_base_factor(hdp);
    
    stListIterator* child_dp_iter = stList_getIterator(base_dp->children);
    DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    while (child_dp != NULL) {
        init_factors_internal(child_dp, root_factor, data_pt_fctr_lists);
        child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    }
    stList_destructIterator(child_dp_iter);
    
    for (int i = 0; i < num_dps; i++) {
        if (data_pt_fctr_lists[i] != NULL) {
            stList_destruct(data_pt_fctr_lists[i]);
        }
    }
    free(data_pt_fctr_lists);
    
    double mean, sum_sq_devs;
    int64_t num_data;
    get_factor_stats(root_factor, &mean, &sum_sq_devs, &num_data);
    add_update_base_factor_params(root_factor, mean, sum_sq_devs, (double) num_data);

    int64_t fctr_child_count;
    DirichletProcess* dp;
    Factor* fctr;
    stSetIterator* fctr_iter;
    for (int64_t i = 0; i < num_dps; i++) {
        dp = dps[i];
        
        fctr_child_count = 0;

        fctr_iter = stSet_getIterator(dp->factors);
        fctr = (Factor*) stSet_getNext(fctr_iter);
        while (fctr != NULL) {
            fctr_child_count += stSet_size(fctr->children);
            fctr = (Factor*) stSet_getNext(fctr_iter);
        }
        stSet_destructIterator(fctr_iter);

        dp->num_factor_children = fctr_child_count;
    }
    //printf("base dp starts with count %d\n", hdp->base_dp->num_factor_children);
}

void finalize_data(HierarchicalDirichletProcess* hdp) {
    verify_valid_dp_assignments(hdp);
    mark_observed_dps(hdp);
    init_factors(hdp);
}

void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int64_t child_id, int64_t parent_id) {
    if (hdp->finalized) {
        fprintf(stderr, "Hierarchical Dirichlet process structure has been finalized. Cannot set new parent.\n");
        exit(EXIT_FAILURE);
    }

    if (child_id >= hdp->num_dps || parent_id >= hdp->num_dps || child_id < 0 || parent_id < 0) {
        fprintf(stderr, "Dirichlet process ID does not exist.\n");
        exit(EXIT_FAILURE);
    }

    DirichletProcess* child_dp = hdp->dps[child_id];
    DirichletProcess* parent_dp = hdp->dps[parent_id];

    if (child_dp->parent != NULL) {
        fprintf(stderr, "Dirichlet process already has parent.\n");
        exit(EXIT_FAILURE);
    }

    child_dp->parent = parent_dp;
    stList_append(parent_dp->children, (void*) child_dp);
}

void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int64_t* dp_ids, int64_t length) {
    if (hdp->data != NULL) {
        fprintf(stderr, "Hierarchical Dirichlet process must be reset before passing new data.\n");
        exit(EXIT_FAILURE);
    }

    hdp->data = data;
    hdp->data_pt_dp_id = dp_ids;
    hdp->data_length = length;

    if (hdp->finalized) {
        finalize_data(hdp);
    }
}

void finalize_hdp_structure(HierarchicalDirichletProcess* hdp) {
    establish_base_dp(hdp);
    verify_dp_tree(hdp);
    verify_tree_depth(hdp, hdp->base_dp, 0, hdp->depth - 1);

    if (hdp->data != NULL) {
        finalize_data(hdp);
    }
    hdp->finalized = true;
}

void reset_hdp_data(HierarchicalDirichletProcess* hdp) {
    if (hdp->data == NULL && hdp->data_pt_dp_id == NULL) {
        return;
    }
    
    //printf("attempt to free data and dp ids\n");
    free(hdp->data);
    hdp->data = NULL;

    free(hdp->data_pt_dp_id);
    hdp->data_pt_dp_id = NULL;

    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;
    
    //printf("attempt to clear factor tree\n");
    
    destroy_dir_proc_factor_tree(hdp->base_dp);
    
    //printf("attempt to clear sampled distributions\n");
    DirichletProcess* dp;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];

        free(dp->posterior_predictive);
        dp->posterior_predictive = NULL;

        free(dp->spline_slopes);
        dp->spline_slopes = NULL;

        dp->observed = false;
    }
    
    hdp->splines_finalized = false;
    
    hdp->samples_taken = 0;
    
    if (hdp->sample_gamma) {
        //printf("attempt to reset gamma and auxilliary variables\n");
        double* gamma = hdp->gamma;
        double* gamma_alpha = hdp->gamma_alpha;
        double* gamma_beta = hdp->gamma_beta;
        
        for (int depth = 0; depth < hdp->depth; depth++) {
            gamma[depth] = gamma_alpha[depth] / gamma_beta[depth];
        }
        
        double* w = hdp->w_aux_vector;
        bool* s = hdp->s_aux_vector;
        
        for (int i = 0; i < num_dps; i++) {
            w[i] = 1.0;
            s[i] = false;
        }
    }
}

void unassign_from_parent(Factor* fctr) {
    if (fctr->factor_type == BASE) {
        fprintf(stderr, "Cannot unassign base factor's parent.\n");
        exit(EXIT_FAILURE);
    }
    //printf("attempts to get parent\n");
    Factor* parent = fctr->parent;
    //printf("attempts to get base factor\n");
    Factor* base_fctr = get_base_factor(parent);
    //printf("attempts to get base dp\n");
    DirichletProcess* base_dp = base_fctr->dp;
    
    //printf("attempts to remove from parent's children\n");
    stSet_remove(parent->children, (void*) fctr);
    fctr->parent = NULL;
    //printf("attempts to decrement count\n");
    (parent->dp->num_factor_children)--;
    
    if (stSet_size(parent->children) == 0) {
        //printf("attempts to destroy parent\n");
        destroy_factor(parent);
    }

    int64_t num_reassign;
    double mean_reassign;
    double sum_sq_devs;
    
    //printf("attempts to get factor stats\n");
    get_factor_stats(fctr, &mean_reassign, &sum_sq_devs, &num_reassign);
    
    // check to see if base factor has been destroyed
    if (stSet_search(base_dp->factors, (void*) base_fctr) != NULL) {
        //printf("attempts to update base factor params\n");
        remove_update_base_factor_params(base_fctr, mean_reassign, sum_sq_devs, (double) num_reassign);
    }
    
    //printf("attempts to cache stats in dp\n");
    DirichletProcess* dp = fctr->dp;
    if (dp != NULL) {
        dp->cached_factor_mean = mean_reassign;
        dp->cached_factor_size = num_reassign;
        dp->cached_factor_sum_sq_dev = sum_sq_devs;
    }
}

void assign_to_parent(Factor* fctr, Factor* parent) {
    if (fctr->factor_type == BASE) {
        fprintf(stderr, "Cannot assign base factor to a parent.\n");
        exit(EXIT_FAILURE);
    }
    
    if (parent->factor_type == DATA_PT) {
        fprintf(stderr, "Cannot assign data point factor to be parent.\n");
        exit(EXIT_FAILURE);
    }
    
    fctr->parent = parent;
    stSet_insert(parent->children, (void*) fctr);
    (parent->dp->num_factor_children)++;

    Factor* base_fctr = get_base_factor(parent);
    if (base_fctr == NULL) {
        return;
    }

    if (fctr->factor_type == DATA_PT) {
        double data_pt = get_factor_data_pt(fctr);
        add_update_base_factor_params(base_fctr, data_pt, 0.0, 1.0);
    }
    else {
        DirichletProcess* dp = fctr->dp;
        add_update_base_factor_params(base_fctr, dp->cached_factor_mean, dp->cached_factor_sum_sq_dev,
                                      (double) dp->cached_factor_size);
    }
}

Factor* sample_from_data_pt_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Attempted a data point factor sample from non-data point factor.\n");
        exit(EXIT_FAILURE);
    }
    
    //printf("attempt to get pool of parents to sample\n");
    stSet* pool = dp->factors;
    int64_t num_fctrs = stSet_size(pool);
    
    //printf("attempt to malloc arrays\n");
    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);

    double* cdf = (double*) malloc(sizeof(double) * (num_fctrs + 1));
    double cumul = 0.0;

    //printf("attempt to get iterator\n");
    stSetIterator* pool_iter = stSet_getIterator(pool);
    Factor* fctr_option;
    double fctr_size;
    //printf("attempt to compute cdf and get factor order\n");
    for (int64_t i = 0; i < num_fctrs; i++) {
        fctr_option = (Factor*) stSet_getNext(pool_iter);
        fctr_order[i] = fctr_option;
        
        fctr_size = (double) stSet_size(fctr_option->children);
        cumul += fctr_size * data_pt_factor_parent_likelihood(fctr, fctr_option);
        cdf[i] = cumul;
    }
    stSet_destructIterator(pool_iter);
    
    //printf("attempt to get gamma\n");
    double gamma = *(dp->gamma);
    //printf("attempt to compute unobserved factor likelihood\n");
    cumul += gamma * unobserved_factor_likelihood(fctr, dp);
    cdf[num_fctrs] = cumul;
    
    //printf("attempt to choose a factor\n");
    int64_t choice_idx = bisect_left(rand_uniform(cumul), cdf, num_fctrs + 1);
    free(cdf);
    
    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            //printf("attempt to make new base factor\n");
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            //printf("attempt to make a new middle factor\n");
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_data_pt_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
    }

    free(fctr_order);

    return fctr_choice;
}

Factor* sample_from_middle_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != MIDDLE) {
        fprintf(stderr, "Attempted a middle factor sample from non-middle factor.\n");
        exit(EXIT_FAILURE);
    }
    
    stSet* pool = dp->factors;
    int64_t num_fctrs = stSet_size(pool);

    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);
    double* log_probs = (double*) malloc(sizeof(double) * num_fctrs);
    
    stSetIterator* pool_iter = stSet_getIterator(pool);
    Factor* fctr_option;
    for (int64_t i = 0; i < num_fctrs; i++) {
        fctr_option = (Factor*) stSet_getNext(pool_iter);
        fctr_order[i] = fctr_option;
        log_probs[i] = factor_parent_joint_log_likelihood(fctr, fctr_option);
    }
    stSet_destructIterator(pool_iter);

    double* cdf = (double*) malloc(sizeof(double) * (num_fctrs + 1));
    double cumul = 0.0;
    double normalizing_const = median(log_probs, num_fctrs);
    
    double fctr_size;
    for (int64_t i = 0; i < num_fctrs; i++) {
        fctr_size = (double) stSet_size(fctr_option->children);
        cumul += fctr_size * exp(log_probs[i] - normalizing_const);
        cdf[i] = cumul;
    }
    free(log_probs);

    double unseen_fctr_prob = unobserved_factor_likelihood(fctr, dp);
    double gamma = *(dp->gamma);

    cumul += gamma * unseen_fctr_prob * exp(-normalizing_const);
    cdf[num_fctrs] = cumul;

    int64_t choice_idx = bisect_left(rand_uniform(cumul), cdf, num_fctrs + 1);
    free(cdf);

    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_middle_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
    }
    free(fctr_order);

    return fctr_choice;
}

Factor* sample_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type == DATA_PT) {
        //printf("attempt to sample parent of data pt factor\n");
        return sample_from_data_pt_factor(fctr, dp);
    }
    else if (fctr->factor_type == MIDDLE) {
        //printf("attempt to sample parent of middle factor\n");
        return sample_from_middle_factor(fctr, dp);
    } else {
        fprintf(stderr, "Cannot sample base factor parent assignments.\n");
        exit(EXIT_FAILURE);
    }
}

void gibbs_factor_iteration(Factor* fctr) {
    //printf("begins gibbs iteration\n");
    DirichletProcess* parent_dp = fctr->parent->dp;
    //printf("going to sample from factors in dp with id = %d\n", parent_dp->id);
    unassign_from_parent(fctr);
    //printf("finishes unassignment\n");
    Factor* new_parent = sample_factor(fctr, parent_dp);
    //printf("finishes sampling parents\n");
    assign_to_parent(fctr, new_parent);
    //printf("finishes gibbs iteration by assigning to parent\n");
}

void cache_prior_contribution(DirichletProcess* dp, double parent_prior_prod) {
    if (!(dp->observed)) {
        return;
    }
    double gamma = *(dp->gamma);
    double total_children = (double) dp->num_factor_children;
    double prior_prod = (gamma / (gamma + total_children)) * parent_prior_prod;
    dp->base_factor_wt += prior_prod;

    stListIterator* child_iter = stList_getIterator(dp->children);
    DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
    while (child != NULL) {
        cache_prior_contribution(child, prior_prod);
        child = (DirichletProcess*) stList_getNext(child_iter);
    }
    stList_destructIterator(child_iter);
}

void cache_base_factor_weight(Factor* fctr) {
    DirichletProcess* dp = fctr->dp;
    
    //printf("attempt to get gamma at depth %d\n", dp->depth);
    double gamma = *(dp->gamma);
    //printf("attempt to get dp size\n");
    double total_children = (double) dp->num_factor_children;
    //printf("attempt to calc weight\n");;
    double wt = ((double) stSet_size(fctr->children)) / (gamma + total_children);
    //printf("attempt to store weight\n");
    dp->base_factor_wt += wt;
    
    // stop pushing the weight cache before reaching data point factors
    if (stList_length(dp->children) > 0) {
        //printf("attempt to store in cache\n");
        stSetIterator* child_fctr_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            cache_base_factor_weight(child_fctr);
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        
        //printf("attempt cache prior portion\n");
        stListIterator* child_dp_iter = stList_getIterator(dp->children);
        DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        while (child_dp != NULL) {
            cache_prior_contribution(child_dp, wt);
            child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        }
        stList_destructIterator(child_dp_iter);
    }
}

void push_factor_distr(DirichletProcess* dp, double* distr, int64_t length) {
    double* sample_collector = dp->posterior_predictive;
    double wt = dp->base_factor_wt;

    for (int64_t i = 0; i < length; i++) {
        sample_collector[i] += wt * distr[i];
    }

    dp->base_factor_wt = 0.0;

    stListIterator* child_iter = stList_getIterator(dp->children);
    DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
    while (child != NULL) {
        if (child->observed) {
            push_factor_distr(child, distr, length);
        }
        child = (DirichletProcess*) stList_getNext(child_iter);
    }
    stList_destructIterator(child_iter);
}

void take_distr_sample(HierarchicalDirichletProcess* hdp) {
    DirichletProcess* base_dp = hdp->base_dp;

    //printf("OOO attempt to allocate vectors\n");
    double* grid = hdp->sampling_grid;
    int64_t length = hdp->grid_length;
    double* pdf = (double*) malloc(sizeof(double) * length);
    
    SumOfLogsMemo* log_sum_memo = hdp->log_sum_memo;

    stSetIterator* base_fctr_iter = stSet_getIterator(base_dp->factors);
    Factor* base_fctr = (Factor*) stSet_getNext(base_fctr_iter);
    while (base_fctr != NULL) {
        //printf("OOO attempt to cache weight\n");
        cache_base_factor_weight(base_fctr);
        //printf("OOO attempt to eval post pred\n");
        evaluate_posterior_predictive(base_fctr, grid, pdf, length, log_sum_memo);
        //printf("OOO attempt to push distr\n");
        push_factor_distr(base_dp, pdf, length);

        base_fctr = (Factor*) stSet_getNext(base_fctr_iter);
    }
    stSet_destructIterator(base_fctr_iter);
    
    //printf("OOO attempt to cache prior weight\n");
    cache_prior_contribution(base_dp, 1.0);
    //printf("OOO attempt to eval prior pred\n");
    evaluate_prior_predictive(hdp, grid, pdf, length);
    //printf("OOO attempt to push prior distr\n");
    push_factor_distr(base_dp, pdf, length);

    (hdp->samples_taken)++;
    
    free(pdf);
}

// Knuth shuffle algorithm
DirichletProcess** get_shuffled_dps(HierarchicalDirichletProcess* hdp) {
    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess** shuffled_dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    int64_t pos;
    for (int64_t i = 0; i < num_dps; i++) {
        pos = rand() % (i + 1);
        shuffled_dps[i] = shuffled_dps[pos];
        shuffled_dps[pos] = dps[i];
    }
    return shuffled_dps;
}

void sample_dp_factors(DirichletProcess* dp, int64_t* iter_counter, int64_t burn_in, int64_t thinning,
                       int64_t* sample_counter, int64_t num_samples) {

    if (!dp->observed) {
        return;
    }

    int64_t iter = *iter_counter;
    int64_t samples_taken = *sample_counter;
    
    // have to pre-allocate the array of sampling factors in case reassignment triggers
    // destruction of the set the iterator is iterating through
    int64_t num_factor_children = dp->num_factor_children;
    Factor** sampling_fctrs = (Factor**) malloc(sizeof(Factor*) * num_factor_children);
    int64_t i = 0;

    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);

    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    stSetIterator* child_fctr_iter;
    Factor* child_fctr;
    while (fctr != NULL) {
        child_fctr_iter = stSet_getIterator(fctr->children);
        child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            sampling_fctrs[i] = child_fctr;
            i++;
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        fctr = (Factor*) stSet_getNext(fctr_iter);
    }
    //printf("### attempts to destruct iterator \n");
    stSet_destructIterator(fctr_iter);
    if (i != num_factor_children) {
        exit(EXIT_FAILURE);
    }
    
    for (int64_t j = 0; j < num_factor_children; j++) {
        //printf("attempt to sample parents of factor %d out of %d\n", j + 1, num_factor_children);
        gibbs_factor_iteration(sampling_fctrs[j]);
        iter++;
        
        //printf("attempt to check iteration count \n");
        if (iter > burn_in && iter % thinning == 0) {
            //printf("attempt to take sample\n");
            take_distr_sample(dp->hdp);
            samples_taken++;
            
            if (samples_taken >= num_samples) {
                break;
            }
        }
    }
    free(sampling_fctrs);
    
    //printf("### attempts to update counters\n");
    *sample_counter = samples_taken;
    *iter_counter = iter;
    //printf("### finished dp factors sampling function\n");
}

double sample_auxilliary_w(DirichletProcess* dp) {
    return (double) genbet((float) *(dp->gamma) + 1.0, (float) dp->num_factor_children);
}

bool sample_auxilliary_s(DirichletProcess* dp) {
    double num_children = (double) dp->num_factor_children;
    return rand_bernoulli(num_children / (num_children + *(dp->gamma)));
}

void sample_gamma_aux_vars(HierarchicalDirichletProcess* hdp) {
    double* w = hdp->w_aux_vector;
    bool* s = hdp->s_aux_vector;

    DirichletProcess** dps = hdp->dps;
    int num_dps = hdp->num_dps;
    DirichletProcess* dp;
    for (int id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        w[id] = sample_auxilliary_w(dp);
        s[id] = sample_auxilliary_s(dp);
    }
}

void sample_base_gamma_internal(HierarchicalDirichletProcess* hdp, double log_w, int64_t num_factors) {
    // Escobar and West's (1995) algorithm
    DirichletProcess* base_dp = hdp->base_dp;
    double gamma_alpha = hdp->gamma_alpha[0];
    double gamma_beta = hdp->gamma_beta[0];

    double num_children = (double) base_dp->num_factor_children;

    double gamma_beta_post = gamma_beta - log_w;
    double gamma_alpha_post = gamma_alpha + (double) num_factors;

    double frac = (gamma_alpha_post - 1.0)
                  / (num_children * gamma_beta_post);

    double wt = frac / (1.0 + frac);
    // note: different parameterization switches alpha and beta
    float sample_gamma = wt * gengam(gamma_beta_post, gamma_alpha_post)
                         + (1 - wt) * gengam(gamma_beta_post, gamma_alpha_post - 1.0);

    hdp->gamma[0] = (double) sample_gamma;
}

void sample_middle_gammas_internal(HierarchicalDirichletProcess* hdp, int64_t depth,
                         double sum_log_w, int64_t sum_s, int64_t num_depth_fctrs) {
    double gamma_alpha = hdp->gamma_alpha[depth];
    double gamma_beta = hdp->gamma_beta[depth];

    float gamma_alpha_post = (float) (gamma_alpha + (double) (num_depth_fctrs - sum_s));
    float gamma_beta_post = (float) (gamma_beta - sum_log_w);
    // note: different parameterization switches alpha and beta
    hdp->gamma[depth] = (double) gengam(gamma_beta_post, gamma_alpha_post);
}

void sample_gammas(HierarchicalDirichletProcess* hdp, int64_t* iter_counter, int64_t burn_in,
                   int64_t thinning, int64_t* sample_counter, int64_t num_samples) {
    int64_t iter = *iter_counter;
    int64_t samples_taken = *sample_counter;
    
    //printf("start sampling gammas\n");

    int64_t tree_depth = hdp->depth;
    double* w = hdp->w_aux_vector;
    bool* s = hdp->s_aux_vector;

    int64_t* num_depth_fctrs = (int64_t*) malloc(sizeof(int64_t) * tree_depth);
    double* sum_log_w = (double*) malloc(sizeof(double) * tree_depth);
    int64_t* sum_s = (int64_t*) malloc(sizeof(int64_t) * tree_depth);

    for (int64_t depth = 0; depth < tree_depth; depth++) {
        num_depth_fctrs[depth] = 0;
        sum_log_w[depth] = 0.0;
        sum_s[depth] = 0;
    }

    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    int64_t dp_depth;
    for (int64_t id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        dp_depth = dp->depth;
        num_depth_fctrs[dp_depth] += stSet_size(dp->factors);
        sum_log_w[dp_depth] += log(w[id]);
        if (s[id]) sum_s[dp_depth]++;
    }

    for (int64_t depth = 0; depth < tree_depth; depth++) {
        if (depth == 0) {
            sample_base_gamma_internal(hdp, sum_log_w[depth], num_depth_fctrs[depth]);
        }
        else {
            sample_middle_gammas_internal(hdp, depth, sum_log_w[depth],
                                          sum_s[depth], num_depth_fctrs[depth]);
        }
        iter++;

        if (iter > burn_in && iter % thinning == 0) {
            take_distr_sample(dp->hdp);
            samples_taken++;

            if (samples_taken >= num_samples) {
                break;
            }
        }
    }
    free(sum_log_w);
    free(sum_s);
    free(num_depth_fctrs);

    *iter_counter = iter;
    *sample_counter = samples_taken;
    
    //printf("finish sampling gammas\n");
}

void sample_gamma_params(HierarchicalDirichletProcess* hdp, int64_t* iter_counter, int64_t burn_in,
                         int64_t thinning, int64_t* sample_counter, int64_t num_samples) {
    sample_gamma_aux_vars(hdp);
    sample_gammas(hdp, iter_counter, burn_in, thinning, sample_counter, num_samples);
}

void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp,
                            int64_t num_samples,
                            int64_t burn_in,
                            int64_t thinning) {
    if (hdp->data == NULL || hdp->data_pt_dp_id == NULL) {
        fprintf(stderr, "Cannot perform Gibbs sampling before passing data to HDP.\n");
        exit(EXIT_FAILURE);
    }

    int64_t iter_counter = 0;
    int64_t sample_counter = 0;
    int64_t num_dps = hdp->num_dps;

    DirichletProcess** sampling_dps;
    while (sample_counter < num_samples) {
        
        sampling_dps = get_shuffled_dps(hdp);

        for (int64_t i = 0; i < num_dps; i++) {
            sample_dp_factors(sampling_dps[i], &iter_counter, burn_in, thinning,
                              &sample_counter, num_samples);
            if (sample_counter >= num_samples) {
                break;
            }
        }
        
        //printf("finishes looping through dps on one iteration\n");

        free(sampling_dps);

        if (hdp->sample_gamma && sample_counter < num_samples) {
            sample_gamma_params(hdp, &iter_counter, burn_in, thinning, &sample_counter,
                                num_samples);
        }
        
        //printf("finishes sampling gamma params on one iteration\n");
    }
    //printf("finishes sampling\n");
}

void finalize_distributions(HierarchicalDirichletProcess* hdp) {
    if (hdp->samples_taken <= 0) {
        fprintf(stderr, "Must perform Gibbs sampling before finalizing sampled distributions.\n");
        exit(EXIT_FAILURE);
    }
    
    if (hdp->splines_finalized) {
        fprintf(stderr, "Distributions have already been finalized.\n");
        exit(EXIT_FAILURE);
    }

    double inv_sample_size = 1.0 / ((double) hdp->samples_taken);
    int64_t grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;

    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    double* distr;
    for (int64_t id = 0; id < num_dps; id++){
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }

        distr = dp->posterior_predictive;

        for (int i = 0; i < grid_length; i++) {
            distr[i] = distr[i] * inv_sample_size;
        }

        dp->spline_slopes = spline_knot_slopes(grid, distr, grid_length);
    }
    
    hdp->splines_finalized = true;
}

double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int64_t dp_id) {
    if (!hdp->splines_finalized) {
        fprintf(stderr, "Must finalize distributions before querying densities.\n");
        exit(EXIT_FAILURE);
    }

    if (dp_id < 0 || dp_id >= hdp->num_dps) {
        fprintf(stderr, "Hierarchical Dirichlet process has no Dirichlet process with this ID.\n");
        exit(EXIT_FAILURE);
    }

    DirichletProcess* dp = hdp->dps[dp_id];
    while (!dp->observed) {
        dp = dp->parent;
    }

    double interp =  grid_spline_interp(x, hdp->sampling_grid, dp->posterior_predictive,
                                        dp->spline_slopes, hdp->grid_length);
    if (interp > 0.0) {
        return interp;
    }
    else {
        return 0.0;
    }
}
