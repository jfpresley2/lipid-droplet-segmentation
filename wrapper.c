#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>

 
/*#define SEED 2343*/

#define BUF_SIZE 10000


 int poisson(double mean) {
  static int initialize = 0;
  static gsl_rng *rand;
  if (!initialize) {
     rand = gsl_rng_alloc(gsl_rng_mt19937);
     gsl_rng_set(rand, time(NULL));
     initialize = 1;
  }
  return (int)gsl_ran_poisson(rand, mean);
}

double gaussian(double mean, double sigma) {
  static int initialize = 0;
  static gsl_rng *rand;
  if (!initialize) {
    rand = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rand, time(NULL));
    initialize = 1;
  }
  return gsl_ran_gaussian(rand, sigma) + mean;
}

/* wrapper for gnu matrix solver for D.
   Will output a solution vector with cols elements.
   This memory will not be known to the D garbage collector!!
   Free it manually once finished with it!!

   For fitting circles (the first use of this), this function will be 
   called only once per circle, to determine the initial values.  This
   will be an insignificant part of the total run time, so whether calling
   the library function or writing something specialized for a 3x3 matrix
   that might be faster is of little consequence.  The library routine should
   be safer to use esp. with special cases.
 */
void solveMatrix(double* output_x, double* m, int cols, int rows, double* b) {
  /*printf("C--solveMatrix entered.  %d cols;  %d rows.\n", cols, rows);*/
  gsl_matrix_view m2 = gsl_matrix_view_array(m,cols,rows);
  gsl_vector_view b2 = gsl_vector_view_array(b,cols);
  gsl_vector *x      = gsl_vector_alloc(cols);
  /*double* output_x  = calloc(cols, sizeof(double));*/

  int s, i;
  
  gsl_permutation *p = gsl_permutation_alloc(cols);
  gsl_linalg_LU_decomp(&m2.matrix, p, &s);
  gsl_linalg_LU_solve(&m2.matrix, p, &b2.vector, x);
  
  /*printf("x= \n");
  gsl_vector_fprintf(stdout, x, "%g");
  printf("C--matrix solved, values printed.\n"); */
  for (i=0; i<cols; i++) {
    output_x[i] = gsl_vector_get(x, i);
  }
  gsl_permutation_free(p);
  gsl_vector_free(x);
  /*printf("C--exiting back to D\n");*/
  return;
}

/*
int main (void) {
  double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

  double b_data[]  = { 1.0, 2.0, 3.0, 4.0 };

  gsl_matrix_view m = gsl_matrix_view_array(a_data, 4, 4);
  gsl_vector_view b = gsl_vector_view_array(b_data, 4);
  gsl_vector *x     = gsl_vector_alloc(4);
  
  int s;

  gsl_permutation *p = gsl_permutation_alloc(4);
  gsl_linalg_LU_decomp(&m.matrix, p, &s);
  gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
  
  printf("x= \n");
  gsl_vector_fprintf(stdout, x, "%g");
  gsl_permutation_free(p);
  gsl_vector_free(x);
  return 0;
}
*/
