// -----------------------------------------------------------------//
// header for math_functions.c                                      //			                                      //
//------------------------------------------------------------------//
#include <math.h>
#define TIME_DELAY 12
// empirical mean of vector entries
double mean(double * input, int nb_samples);
// center vector by substracting mean
void center(double * input, int nb_samples);
// dot product of two vectors
double dot_product(double * input_1, double * input_2, int nb_samples);
// correlation of two vectors
double correlation(double * input_1, double * input_2, int nb_samples);
// vector norm
double norm(double * input, int nb_samples);
// vector standard deviation
double std(double * input, int nb_samples);
// Matrix Trace
double tr(double (*M)[TIME_DELAY], int dim_M);
// Off function
double off(double (*M)[TIME_DELAY], int dim_M);
// Correlation Matrix
void compute_correlation(double (*M)[TIME_DELAY], double * input1, double * input2, int n_samples);
