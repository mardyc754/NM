#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nrutil.h"


int main(void)
{
    srand(time(NULL));
    float T = 1.0;
    float t_max = 3.0*T;
    int N[3] = {8,10,12};

    for(int k=0; k<3; k++)
    {
        /*FILE *fp = N[k] == 8 ? fopen("k8.dat", "w") : (N[k] == 10 ? fopen("k10.dat", "w") : fopen("k12.dat", "w"));*/ 
        FILE *fp;
        if(N[k]==8)
        {
            fp = fopen("k8.dat", "w");
        }
        else if(N[k]==10)
        {
            fp = fopen("k10.dat", "w");
        }
        else
        {
            fp = fopen("k12.dat", "w");
        }
        int N_k = (int)pow(2,N[k]);
        float *f0 = vector(1,2*N_k); // sygnał bez szumu (5)
        float *f = vector(1,2*N_k); // sygnał z szumem (4)
        float *g1 = vector(1,2*N_k); // tablice 
        float *g2 = vector(1,2*N_k); // funkcji wagowej (6)
        float *g = vector(1, 2*N_k);
        float dt = t_max/N_k;
        float sigma = T/20.0;
        printf("k = %d, dt = %g\n", N[k], dt);

        for(int i=1; i<=N_k; i++)
        {
            float omega = 2*M_PI/T;
            float t = dt*(i-1);
            f0[2*i-1] = sin(omega*t) + sin(2*omega*t) + sin(3*omega*t);
            f[2*i-1] = f0[2*i-1] + (float)rand() / RAND_MAX - 0.5;
            g1[2*i-1] = g2[2*i-1] = 1.0/(sigma*sqrt(2*M_PI)) * exp(-(t*t)/(2*sigma*sigma));
            fprintf(fp, "%15g %15g\n", t, f[2*i-1]); 
            
            g2[2*i] = g1[2*i] = f[2*i] = f0[2*i] = 0.0;
        }

        four1(f, N_k, 1);
        four1(g1, N_k, 1);
        four1(g2, N_k, -1);

        for(int i=1; i<=2*N_k; i++)
        {
            g[i] = g1[i] + g2[i];
        }

        for(int i=1; i<=N_k; i++)
        {
            float a_1 = f[2*i-1];
            float b_1 = f[2*i];
            float a_2 = g[2*i-1];
            float b_2 = g[2*i];
            f[2*i-1] = a_1 * a_2 - b_1 * b_2;
            f[2*i] = a_1 * b_2 + a_2 * b_1;
        }

        four1(f, N_k, -1);

        float f_max = 0.0;
        for(int i=1; i<=N_k; i++)
        {
            f_max = fabs(f[2*i-1]) > f_max ? fabs(f[2*i-1]) : f_max;
        }

        printf("k = %d, f_max = %g\n\n", N[k], f_max);
        fprintf(fp, "\n\n");
        for(int i=1; i<=N_k; i++)
        {
            fprintf(fp, "%15g %15g\n", dt*(i-1), f[2*i-1]/(f_max/2.5));
        }

        free_vector(f0, 1, 2*N_k);
        free_vector(f, 1, 2*N_k);
        free_vector(g1, 1, 2*N_k);
        free_vector(g2, 1, 2*N_k);
        free_vector(g, 1, 2*N_k);
        fclose(fp);
    }
}