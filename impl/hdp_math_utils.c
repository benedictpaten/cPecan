#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "hdp_math_utils.h"

#define LOG_ROOT_PI 0.572364942924700087071713
#define LOG_4 1.386294361119890618834464

typedef struct LogGammaHalfMemo LogGammaHalfMemo;

struct LogGammaHalfMemo {
    double alpha;
    double* zero_offset_memo;
    int zero_offset_final_entry;
    int zero_offset_length;
    double* half_offset_memo;
    int half_offset_final_entry;
    int half_offset_length;
};

LogGammaHalfMemo* new_log_gamma_memo(double alpha) {
    LogGammaHalfMemo* memo = (LogGammaHalfMemo*) malloc(sizeof(LogGammaHalfMemo));
    memo->alpha = alpha;
    double* zero_base_case = (double*) malloc(sizeof(double));
    zero_base_case[0] = lgamma(alpha);
    memo->zero_offset_final_entry = 0;
    memo->zero_offset_memo = zero_base_case;
    memo->zero_offset_length = 1;

    double* half_base_case = (double*) malloc(sizeof(double));
    half_base_case[0] = lgamma(alpha + .5);
    memo->half_offset_final_entry = 0;
    memo->half_offset_memo = half_base_case;
    memo->half_offset_length = 1;

    return memo;
}

void destroy_log_gamma_memo(LogGammaHalfMemo* memo) {
    free(memo->half_offset_memo);
    free(memo->zero_offset_memo);
    free(memo);
}

void extend_gamma_zero_offset_memo(LogGammaHalfMemo* memo) {
    int final_entry = memo->half_offset_final_entry + 1;
    memo->zero_offset_final_entry = final_entry;
    double* current_array = memo->zero_offset_memo;

    int current_length = memo->zero_offset_length;
    if (current_length == final_entry) {

        int new_array_length = current_length * 2;
        double* new_array = (double*) malloc(sizeof(double) * new_array_length);

        for (int i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }

        memo->zero_offset_length = new_array_length;
        memo->zero_offset_memo = new_array;
        free(current_array);
        current_array = new_array;
    }
    double log_term = log(memo->alpha - 1.0 + (double) final_entry);
    current_array[final_entry] = current_array[final_entry - 1] + log_term;
}

void extend_gamma_half_offset_memo(LogGammaHalfMemo* memo) {
    int final_entry = memo->half_offset_final_entry + 1;
    memo->half_offset_final_entry = final_entry;
    double* current_array = memo->half_offset_memo;

    int current_length = memo->half_offset_length;
    if (current_length == final_entry) {

        int new_array_length = current_length * 2;
        double* new_array = (double*) malloc(sizeof(double) * new_array_length);

        for (int i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }

        memo->half_offset_length = new_array_length;
        memo->half_offset_memo = new_array;
        free(current_array);
        current_array = new_array;
    }
    double log_term = log(memo->alpha -.5 + (double) final_entry);
    current_array[final_entry] = current_array[final_entry - 1] + log_term;
}

// returns log(Gamma(memo->alpha + n / 2))
double offset_log_gamma_half(int n, LogGammaHalfMemo* memo) {
    int idx = n / 2;
    if (n % 2 == 0) {
        while (memo->zero_offset_final_entry < idx) {
            extend_gamma_zero_offset_memo(memo);
        }
        return memo->zero_offset_memo[idx];
    }
    else {
        while (memo->half_offset_final_entry < idx) {
            extend_gamma_half_offset_memo(memo);
        }
        return memo->half_offset_memo[idx];
    }
}

struct SumOfLogsMemo {
    double* memo_array;
    int64_t final_entry;
    int64_t array_length;
};

SumOfLogsMemo* new_log_sum_memo() {
    SumOfLogsMemo* memo = (SumOfLogsMemo*) malloc(sizeof(SumOfLogsMemo));
    double* base_case = (double*) malloc(sizeof(double));
    base_case[0] = 0.0;
    memo->memo_array = base_case;
    memo->final_entry = 1;
    memo->array_length = 1;
    return memo;
}

void destroy_log_sum_memo(SumOfLogsMemo* memo) {
    free(memo->memo_array);
    free(memo);
}

void extend_log_sum_memo(SumOfLogsMemo* memo) {
    int64_t final_entry = memo->final_entry;
    int64_t current_length = memo->array_length;
    if (current_length == final_entry) {
        double *current_array = memo->memo_array;

        int64_t new_array_length = current_length * 2;
        double *new_array = (double*) malloc(sizeof(double) * new_array_length);

        for (int64_t i = 0; i < current_length; i++) {
            new_array[i] = current_array[i];
        }

        memo->array_length = new_array_length;
        memo->memo_array = new_array;
        free(current_array);
    }
    double log_term = log((double) final_entry + 1);
    memo->memo_array[final_entry] = memo->memo_array[final_entry - 1] + log_term;
    (memo->final_entry)++;
}

double sum_of_logs(SumOfLogsMemo* memo, int64_t n) {
    while (n > memo->final_entry) {
        extend_log_sum_memo(memo);
    }
    return memo->memo_array[n - 1];
}

// returns log(Gamma(n / 2)) in amortized constant time with low risk of overflow
double log_gamma_half(int64_t n, SumOfLogsMemo* sum_of_logs_memo) {
    if (n <= 2) {
        fprintf(stderr, "log_gamma_half only supports n > 2\n");
        exit(EXIT_FAILURE);
    }
    if (n % 2 == 0) {
        return sum_of_logs(sum_of_logs_memo, n / 2 - 1);
    }
    else {
        return LOG_ROOT_PI - (n / 2) * LOG_4 + sum_of_logs(sum_of_logs_memo, n - 1)
               - sum_of_logs(sum_of_logs_memo, n / 2);
    }
}

// quick-select algorithm on array copy (does not alter original array)
double quickselect(double* arr, int length, int target_idx) {
    if (target_idx < 0 || target_idx >= length) {
        fprintf(stderr, "Order statistic outside of array bounds\n");
        exit(EXIT_FAILURE);
    }

    double* arr_copy = (double*) malloc(sizeof(double) * length);
    for (int i = 0; i < length; i++ ) {
        arr_copy[i] = arr[i];
    }

    int low = 0;
    int hi = length - 1;
    int mid;
    int median;
    double temp;

    while (true) {
        // median of three technique
        mid = (hi + low) / 2;
        if (arr_copy[hi] > arr_copy[mid]) {
            if (arr_copy[hi] > arr_copy[low]) {
                if (arr_copy[mid] > arr_copy[low]) {
                    median = mid;
                }
                else {
                    median = low;
                }
            }
            else {
                median = hi;
            }
        }
        else {
            if (arr_copy[hi] > arr_copy[low]) {
                median = hi;
            }
            else {
                if (arr_copy[mid] > arr_copy[low]) {
                    median = low;
                }
                else {
                    median = mid;
                }
            }
        }

        // remove pivot
        temp = arr_copy[median];
        arr_copy[median] = arr_copy[hi];
        arr_copy[hi] = temp;

        // partition array
        int pivot = low;
        for (int i = low; i < hi; i++) {
            if (arr_copy[i] < arr_copy[hi]) {
                temp = arr_copy[i];
                arr_copy[i] = arr_copy[pivot];
                arr_copy[pivot] = temp;
                pivot++;
            }
        }

        temp = arr_copy[pivot];
        arr_copy[pivot] = arr_copy[hi];
        arr_copy[hi] = temp;

        if (pivot == target_idx) {
            return arr_copy[pivot];
        }
        else if (pivot < target_idx) {
            low = pivot + 1;
        }
        else {
            hi = pivot - 1;
        }
    }
}

double median(double* arr, int64_t length) {
    return quickselect(arr, length, length / 2);
}


// returns the index of the first element of arr greater or equal to x, assuming arr is sorted
// returns final index if x is greater than all elements of arr
int64_t bisect_left(double x, double* arr, int64_t length) {
    if (x <= arr[0]) {
        return 0;
    }
    int64_t low = 0;
    int64_t hi = length - 1;
    int64_t mid;
    double arr_mid;
    while (hi > low + 1) {
        mid = (hi + low) / 2;
        
        arr_mid = arr[mid];
        if (x <= arr_mid) {
            hi = mid;
        }
        else {
            low = mid;
        }
    }
    return hi;
}

void spline_knot_slopes_internal(double* x, double* y, double* k, int64_t idx, double center_coef_prev,
                                 double right_coef_prev, double rhs_prev, int64_t final_idx) {

    if (idx == final_idx) {
        double left_coef = 1.0 / (x[idx] - x[idx - 1]);
        double center_coef = 2.0 * left_coef;
        double rhs = 3.0 * (y[idx] - y[idx - 1]) * left_coef * left_coef;
        // Cramer's rule
        k[idx] = (rhs * center_coef_prev - rhs_prev * left_coef) /
                 (center_coef * center_coef_prev - right_coef_prev * left_coef);
        return;
    }

    double left_coef = 1.0 / (x[idx] - x[idx - 1]);
    double right_coef = 1.0 / (x[idx + 1] - x[idx]);
    double center_coef = 2.0 * (left_coef + right_coef);

    double rhs = 3.0 * ((y[idx] - y[idx - 1]) * left_coef * left_coef +
                        (y[idx + 1] - y[idx]) * right_coef * right_coef);

    center_coef -= left_coef * right_coef_prev / center_coef_prev;
    rhs -= left_coef * rhs_prev / center_coef_prev;

    spline_knot_slopes_internal(x, y, k, idx + 1, center_coef, right_coef, rhs, final_idx);

    k[idx] = (rhs - right_coef * k[idx + 1]) / center_coef;
}

double* spline_knot_slopes(double* x, double* y, int64_t length) {
    double* k = (double*) malloc(sizeof(double) * length);

    double right_coef = 1.0 / (x[1] - x[0]);
    double center_coef = 2.0 * right_coef;
    double rhs = 3.0 * (y[1] - y[0]) * right_coef * right_coef;

    spline_knot_slopes_internal(x, y, k, 1, center_coef, right_coef, rhs, length - 1);

    k[0] = (rhs - right_coef * k[1]) / center_coef;

    return k;
}

double spline_interp(double query_x, double* x, double* y, double* slope, int64_t length) {
    if (query_x <= x[0]) {
        return y[0] - slope[0] * (x[0] - query_x);
    }
    else if (query_x >= x[length - 1]) {
        int64_t n = length - 1;
        return y[n] + slope[n] * (query_x - x[n]);
    }
    else {
        int64_t idx_right = bisect_left(query_x, x, length);
        int64_t idx_left = idx_right - 1;

        double dx = x[idx_right] - x[idx_left];
        double dy = y[idx_right] - y[idx_left];

        double a = slope[idx_left] * dx - dy;
        double b = dy - slope[idx_right] * dx;

        double t_left = (query_x - x[idx_left]) / dx;
        double t_right = 1.0 - t_left;

        return t_right * y[idx_left] + t_left * y[idx_right] +
               t_left * t_right * (a * t_right + b * t_left);
    }
}

// assumes even spacing of x points
double grid_spline_interp(double query_x, double* x, double* y, double* slope, int64_t length) {
    if (query_x <= x[0]) {
        return y[0] - slope[0] * (x[0] - query_x);
    }
    else if (query_x >= x[length - 1]) {
        int64_t n = length - 1;
        return y[n] + slope[n] * (query_x - x[n]);
    }
    else {
        double dx = x[1] - x[0];
        int64_t idx_left = (int64_t) ((query_x - x[0]) / dx);
        int64_t idx_right = idx_left + 1;
        
        double dy = y[idx_right] - y[idx_left];
        
        double a = slope[idx_left] * dx - dy;
        double b = dy - slope[idx_right] * dx;
        
        double t_left = (query_x - x[idx_left]) / dx;
        double t_right = 1.0 - t_left;
        
        return t_right * y[idx_left] + t_left * y[idx_right]
               + t_left * t_right * (a * t_right + b * t_left);
    }
}

double double_max(double a, double b) {
    if (a > b) {
        return a;
    }
    else {
        return b;
    }
}

double* linspace(double start, double stop, int64_t length) {
    if (start >= stop) {
        fprintf(stderr, "linspace requires stop > start\n");
        exit(EXIT_FAILURE);
    }
    double* lin = (double*) malloc(sizeof(double) * length);
    int64_t n = length - 1;
    double dx = (stop - start) / ((double) n);
    for (int i = 0; i < n; i++) {
        lin[i] = start +  i * dx;
    }
    lin[n] = stop;
    return lin;
}

double rand_standard_uniform() {
    return ((double) rand()) / ((double) RAND_MAX);
}

double rand_uniform(double a) {
    return ((double) rand()) / ((double) RAND_MAX / a);
}

bool rand_bernoulli(double p) {
    return (rand_standard_uniform() < p);
}

double rand_exponential(double lambda) {
    double draw;
    do {
        draw = rand_standard_uniform();
    } while (draw == 1.0);
    return -log(1.0 - draw) / lambda;
}

double log_posterior_conditional_term(double nu_post, double two_alpha_post,
                                      double beta_post, SumOfLogsMemo* memo) {

    return log_gamma_half((int64_t) two_alpha_post, memo)
           - .5 * (log(nu_post) + two_alpha_post * log(beta_post));
}

void normal_inverse_gamma_params(double* x, int64_t length, double* mu_out, double* nu_out,
                                 double* alpha_out, double* beta_out) {
    double mean = 0.0;
    for (int64_t i = 0; i < length; i++) {
        mean += x[i];
    }
    mean /= (double) length;

    double dev;
    double sum_sq_devs = 0.0;
    for (int64_t i = 0; i < length; i++) {
        dev = x[i] - mean;
        sum_sq_devs += dev * dev;
    }

    *mu_out = mean;
    *nu_out = (double) length;
    *alpha_out = ((double) length - 1.0) / 2.0;
    *beta_out = .5 * sum_sq_devs;
}
