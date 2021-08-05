#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#include "nrutil.c"
#include "gammp.c"
#include "gcf.c"
#include "gammln.c"
#include "gser.c"

#define N 10000
#define PI 3.141527


double average(double *x)
{
    double res = 0.0;
    for(int i=0; i<N; i++)
    {
        res += x[i];
    }
    return res / N;
}

double standard_deviation(double *x)
{
    double res = 0.0;
    double avg = average(x);
    for(int i=0; i<N; i++)
    {
        res += pow(x[i]-avg,2);
    }
    return sqrt(res / N);
}

void set_n(double x_min, double x_max, double *x, int *n, const int k, double delta)
{
    for(int j=0; j<k; j++)
    {
        n[j] = 0;
    }
    for(int i=0; i<N; i++)
    {
        int j = (x[i] - x_min)/delta;
        n[j]++;
    }   
}

void uniform(long a, long c, long m, FILE *fp, FILE *fp1, FILE *fp2)
{
    
    double *x = (double *)malloc(N*sizeof(double));
    long X_0 = 10;
    
    for(int i=0; i<N; i++)
    {
        X_0 = (a*X_0+c) % m;
        x[i] = X_0/(m + 1.0);
    } 

    
    for(int i=0; i<N-1; i++)
    {
        fprintf(fp1, "%15g %15g\n", x[i], x[i+1]);
    }
    fprintf(fp1, "\n\n");


    double mu = average(x);
    double sigma = standard_deviation(x);
    fprintf(fp, "Srednia: %lf\n", mu);
    fprintf(fp, "Odchylenie standardowe: %lf\n", sigma);

    const int k = 12;

    double x_min = 0.0;
    double x_max = 1.0;
    int *n = (int *)malloc(k*sizeof(int));

    double delta = (x_max-x_min)/k;
    set_n(x_min, x_max, x, n, k, delta);

    for(int j=0; j<k; j++)
    {
        double x_j_min = x_min + delta * j;
        double x_j_max = x_min + (j+1) * delta;
        fprintf(fp2, "%15g %15g\n", (x_j_min + x_j_max)/2, n[j]/(1.0*N));
    }
    fprintf(fp2, "\n\n");

    free(x);
    free(n);
}

double f(double x, double mu, double sigma)
{
    return 1.0/(sigma*sqrt(2*PI))*exp(-pow(x-mu,2)/(2*sigma*sigma));
}

double F(double x, double mu, double sigma)
{
    return (1.0 + erf((x-mu)/(sqrt(2.0)*sigma)))/2.0;
}

int main()
{
    FILE *fp, *U, *U_hist;
    fp = fopen("U_N_Info.dat", "w");
    U = fopen("U.dat", "w");
    U_hist = fopen("U_hist.dat", "w");
    
    fprintf(fp, "#### Rozklad jednorodny ####\n");
    fprintf(fp, "\nPrzypadek 1: a=123, c=1, m=2^15\n");
    uniform(123, 1, pow(2,15), fp, U, U_hist);
    
    fprintf(fp, "\nPrzypadek 2: a=69069, c=1, m=2^32\n");
    uniform(69069, 1, pow(2,32),fp, U, U_hist);
    
    fclose(U);
    fclose(U_hist);
    
    fprintf(fp, "\n\n#### Rozklad normalny ####\n");
    fprintf(fp, "a=69069, c=1, m=2^32\n\n");

    long a = 69069;
    long c = 1;
    long m = pow(2,32);
    long X_0 = 10;
    long X_i;
    
    double mu = 0.2;
    double sigma = 0.5;

    double x_min = mu - 3.0*sigma;
    double x_max = mu + 3.0*sigma;

    double *x = (double *)malloc(N*sizeof(double));

    double u_1;
    double u_2;

    for(int i=0; i<N; i++)
    {
        do
        {
            X_i = (a*X_0+c) % m;
            X_0 = (a*X_i+c) % m;

            u_1 = (X_i/(m + 1.0)) * (x_max-x_min) + x_min;
            u_2 = X_0/(m + 1.0);
        } while (u_2 > f(u_1, mu, sigma));
        x[i] = u_1;
    } 

    fprintf(fp, "Srednia: %g\nOdchylenie standardowe: %g\n", 
                average(x), standard_deviation(x));
    fprintf(fp, "Wariancja: %g\n", pow(standard_deviation(x),2));
    
    FILE *fp2 = fopen("N_hist.dat", "w");
    const int k = 12;

    int *n = (int *)malloc(k*sizeof(int));

    double delta = (x_max-x_min)/k;
    set_n(x_min, x_max, x, n, k, delta);

    for(int j=0; j<k; j++)
    {
        double x_j_min = x_min + delta * j;
        double x_j_max = x_min + (j+1) * delta;
        fprintf(fp2, "%15g %15g\n", (x_j_min + x_j_max)/2, n[j]/(1.0*N));
    }
    fclose(fp2);

    double chi_2 = 0.0;
    
    fprintf(fp, "\n\n%2s %15s %15s\n", "j", "p_j", "N * p_j");
    for(int j=0; j<k; j++)
    {
        double x_j_min = x_min + delta * j;
        double x_j_max = x_min + (j+1) * delta;
        double p_j = F(x_j_max, 0.2, 0.5) - F(x_j_min, 0.2, 0.5);
        chi_2 += pow(n[j] - N*p_j, 2)/(N*p_j);
        fprintf(fp, "%2d %15g %15g\n", j, p_j, N*p_j);
    }

    fprintf(fp, "\nStatystyka testowa: %g", chi_2);

    if(chi_2 < 16.91)
    {
        fprintf(fp, " -  Hipoteza nie zostala odrzucona dla alpha = 0.05\n");
    }
    else
    {
        fprintf(fp, " -  Hipoteza zostala odrzucona dla alpha = 0.05\n");
    }

    double nu = 9.0;
    double confidence_level = gammp(nu/2.0, chi_2/2.0);
    double alpha_pr = 1.0 - confidence_level;

    fprintf(fp, "Poziom ufnosci: %g\nPoziom istotnosci: %g\n", confidence_level, alpha_pr);

    free(x);
    free(n);
    fclose(fp);
    return 0;
}