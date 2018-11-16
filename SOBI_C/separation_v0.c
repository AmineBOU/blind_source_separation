///////////////////////////////////////
// separation_v0.c                   //
// Separation of sources with SOBI   //
// Test function: matrix A is known  //
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "math_functions.h"
#define TIME_DELAY 12

void separation_v0(double * input_1, double * input_2, double * output_1, double * output_2,
	                     int nb, int time_delay){ //nb stands for nb_samples, the number of samples
  double normalisation;
  double A11, A12, A21, A22, D;
  int i;

  // Signals normalisation
  normalisation = sqrt(2/(std(input_1,nb)*std(input_1,nb)+std(input_2,nb)*std(input_2,nb)));
  for (int i=0; i<nb; i++){
    input_1[i] *= normalisation;
    input_2[i] *= normalisation;
  }

  // Mixing matrix
  //A11 =  0.25508; A12 = 0.83463;
  //A21 = 0.45562; A22 = -0.74592;
  //D   = A11*A22-A12*A21;

  int window_size = time_delay;

	double R12[TIME_DELAY][TIME_DELAY] = {0};
	double R22[TIME_DELAY][TIME_DELAY] = {0};

	compute_correlation(R11, input_1, input_1,nb);
	compute_correlation(R12, input_1, input_2,nb);
	compute_correlation(R22, input_2, input_2,nb);


  double T1 = tr(R11, window_size);
  double T2 = tr(R22, window_size);
  double T12 = tr(R12, window_size);

  double F1 = off(R11, window_size);
  double F2 = off(R22, window_size);
  double F12 = off(R12, window_size);

  double alpha = 2*F12*T12-F1*T2-F2*T1;
  double beta = 2*(T12*T12-T1*T2);
  double gamma = sqrt((F1*T2-F2*T1)*(F1*T2-F2*T1) +4*(F12*T2-T12*F2)*(F12*T1-T12*F1));


  double d1 = alpha - gamma;
  double d2 = alpha + gamma;

//calcul des coefficients de la matrice Â
  A11 = beta * F1 - T1 * d1;
  A12 = beta * F12 - T12 * d2;
  A21 = beta*F12 - T12 * d1;
  A22 = beta * F2 -T2 * d2;

// compute D
  D   = A11*A22-A12*A21;
  printf("\nA  = [[%f, %f]\n      [%f, %f]]\n\n", A11, A12, A21, A22);

  // Sources separation
  for (i = 0; i<nb; i++) {
    output_1[i] =  (A22*input_1[i] - A12*input_2[i])/D;
    output_2[i] =  (-A21*input_1[i] + A11*input_2[i])/D;
  }

}
