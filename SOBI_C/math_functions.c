// -----------------------------------------------------------------//
// Code property:                                                   //
// ==============                                                   //
// TELECOM BRETAGNE						    //
// Dpt. Signal et Communications				    //
// Technopole Brest-Iroise					    //
// CS 83818 - 29238 Brest Cedex 3 - France			    //
//                                                                  //
// General info:                                                    //
// =============                                                    //
// program: math_functions.c                                        //
// last upate: 14/09/2016                                           //
// info: thierry.chonavel@telecom-bretagne.eu                       //
//       thierry.legall1@telecom-bretagne.eu                        //
//                                                                  //
// Code info:  	                                                    //
// ==========                                                       //
// A bunch of math functions for program Sobi.c                     //
//------------------------------------------------------------------//
#include <math.h>
#include <stdlib.h>
#define TIME_DELAY 12
// empirical mean of vector entries
double mean(double * input, int nb_samples) {
	double sum = 0.0;
	for (int i=0; i<nb_samples; i++) sum += input[i];
  return sum/nb_samples;
}
// center vector by substracting mean
void center(double * input, int nb_samples) {
	double m = mean(input, nb_samples);
	for (int i=0; i<nb_samples; i++) input[i] -= m;
}
// dot product of two vectors
double dot_product(double * input_1, double * input_2, int nb_samples) {
	double dot_prod = 0.0;
	for (int i=0; i<nb_samples; i++) dot_prod += input_1[i]*input_2[i];
	return dot_prod;
}
// correlation of two vectors
double correlation(double * input_1, double * input_2, int nb_samples) {
	double s = 0;
	double m1 = mean(input_1, nb_samples);
	double m2 = mean(input_2, nb_samples);
	for (int i=0; i<nb_samples; i++) s+= (input_1[i]-m1)*(input_2[i]-m2);
	return s/nb_samples;
}
// squared vector norm
double squared_norm(double * input, int nb_samples) {
	return dot_product(input, input, nb_samples);
}
// vector norm
double norm(double * input, int nb_samples) {
	return sqrt(squared_norm(input, nb_samples));
}
// vector standard deviation
double std(double * input, int nb_samples) {
	double s = 0;
	double m = mean(input,nb_samples);
	for (int i=0; i<nb_samples; i++) s+= (input[i]-m)*(input[i]-m);
	return sqrt(s/nb_samples);
}
// Matrix Trace
double tr(double (*M)[TIME_DELAY], int dim_M) {
    double trace = 0;
    for (int i=0; i<dim_M; i++) {
        trace += M[i][i];
    }
    trace /= dim_M;
    return trace;
}
// Off function
double off(double (*M)[TIME_DELAY], int dim_M) {
    double F = 0;
    for (int i=0; i<dim_M; i++) {
        for (int j=0; j<dim_M; j++) {
            if (i != j) {
                F += M[i][j];
            }
        }
    }
    F /= dim_M * (dim_M - 1);
    return F;
}
// Correlation Matrix
void compute_correlation(double (*M)[TIME_DELAY], double * input1, double * input2, int n_samples) {
	int K = n_samples / TIME_DELAY; // Nombre de fenetres
	double aux_vector1[TIME_DELAY]; // vecteur pour socker les elements du premier vecteur
	double aux_vector2[TIME_DELAY]; // vecteur pour socker les elements du second vecteur

	for(int i=0; i < K; i++) {
    int index_count = 0;

		for (int index=i*TIME_DELAY; index < (i+1)*TIME_DELAY; index++ ){
      aux_vector1[index_count] = input1[index];
        aux_vector2[index_count] = input2[index];
        index_count++;
	    }

		for(int j=0; j < TIME_DELAY; j++ ){
      for (int k = 0; k< TIME_DELAY;k++ ){
        M[j][k] = M[j][k] + (aux_vector1[j] * aux_vector2[k]);
      }
    }
	}

	for (int line=0; line< TIME_DELAY ; line++){
    for(int col=0; col<TIME_DELAY; col++){
      M[line][col] =  M[line][col]/K;
    }
	}
}
