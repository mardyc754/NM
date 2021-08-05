#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/usr/local/include/gsl/gsl_math.h"
#include "/usr/local/include/gsl/gsl_linalg.h"
#include "/usr/local/include/gsl/gsl_eigen.h" 
#include "/usr/local/include/gsl/gsl_complex.h" 
#include "/usr/local/include/gsl/gsl_complex_math.h"


double x_i(int i, double L, double dx)
{
  return -L/2 + dx * (i+1);
}

double rho_i(int alpha, double x)
{
  return 1. + 4 * alpha * x * x;
}

void print_matrix(gsl_matrix *M, int size)
{
  for(int i=0; i<size; i++)
  {
    for(int j=0; j<size; j++)
    {
      printf("%20g ", gsl_matrix_get(M, i,j));
    }
      printf("\n");
  }
}

 
double get_real_i(int i, gsl_vector_complex *eval)
{
  return GSL_REAL(gsl_vector_complex_get(eval, i));
}

void write_to_file(const char *filename, gsl_matrix_complex *evec, int n)
{
  double L = 10.;
  double dx = L/(n+1);

  FILE *fp = fopen(filename, "w");
  
  for(int i=0; i<n; i++)
  {
    fprintf(fp, "%lf", x_i(i, L, dx));
    for(int num=0; num<6; num++)
    {
      fprintf(fp, "\t%lf", GSL_REAL(gsl_matrix_complex_get(evec, i, num)));
    }
    fprintf(fp, "\n");
  }
  
  fclose(fp);
}

void eigen_values(gsl_matrix *A, int alpha, int n, gsl_vector_complex *eval,
	  gsl_matrix_complex *evec, gsl_eigen_nonsymmv_workspace *w, FILE *fp)
{
  double L = 10.;
  int N = 1;
  double dx = L/(n+1);
    
    
  gsl_matrix_set_zero(A);
    
  for(int i=0; i<n; i++)
  {
    gsl_matrix_set(A, i, i, 2*N / (rho_i(alpha, x_i(i, L, dx))*dx*dx));
    if(i-1 >= 0)
    {
      gsl_matrix_set(A, i, i-1, (- N) / (rho_i(alpha, x_i(i, L, dx))*dx*dx));
    }
    if(i+1 < n)
    {
      gsl_matrix_set(A, i, i+1, (- N) / (rho_i(alpha, x_i(i, L, dx))*dx*dx));
    }
  }
  
  gsl_eigen_nonsymmv(A, eval, evec, w); 
  gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    
  fprintf(fp, "%d", alpha);
  for(int i=0; i<6; i++)
  {
    fprintf(fp, "\t%lf", sqrt(get_real_i(i, eval))); 
  }
  fprintf(fp, "\n");
    
  if(fabs(alpha)<1e-6)
  {
    write_to_file("out2.txt", evec, n);
  }
  else if(fabs(alpha-100.0)<1e-6)
  {
    write_to_file("out3.txt", evec, n);
  }
    
}


int main(void)
{
  
  int n = 200;
    
  FILE *fp = fopen("out1.txt", "w");
  gsl_matrix *A = gsl_matrix_calloc(n,n);
  gsl_vector_complex *eval = gsl_vector_complex_calloc(n);
	gsl_matrix_complex *evec = gsl_matrix_complex_calloc(n, n);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    
  for(int alpha = 0; alpha<=100; alpha+=2)
  {
    eigen_values(A, alpha, n, eval, evec, w, fp);
  }
  gsl_matrix_free(A);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(eval);
  gsl_eigen_nonsymmv_free(w);
  fclose(fp);
    
  return 0;
} 
