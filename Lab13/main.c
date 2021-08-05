#include <stdio.h> 
#include <math.h> 

#include "nrutil.h"
#include "nrutil.c"
#include "gammln.c"
#include "gauher.c"
#include "gaulag.c"
#include "gauleg.c"


float f1(float x)
{
    return 1.0 / (x*sqrt(x*x-1.0));
}

float f2a(float x)
{
    return 0.5*log(fabs(x));
}

float f2b(float x)
{
    return log(x)*exp(-x*x);
}

float f3(float x)
{
    return sin(2*x)*exp(-2.0*x);
}

int main()
{

    FILE *fp;
    fp = fopen("out.dat", "w");
    float c1a = 3.1415927 / 3.0;
    
    for (int n = 2; n <= 100; n++)
    {
        float *x = vector(1, n);
        float *w = vector(1, n);
        gauleg(1.0, 2.0, x, w, n);
        float c1 = 0.0;
        for(int i=1; i <= n; i++)
        {
            c1 += w[i]*f1(x[i]);
        }
        fprintf(fp, "%4d %15g\n", n, fabs(c1a-c1));

        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }
    fprintf(fp, "\n\n");
    
    float c2a = -0.8700577; 
    for (int n = 2; n <= 100; n += 2)
    {
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauher(x, w, n);
        float c2 = 0.0;
        for(int i=1; i <= n; i++)
        {
            c2 += w[i]*f2a(x[i]);
        }
        fprintf(fp, "%4d %15g\n", n, fabs(c2a-c2));

        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }
    fprintf(fp, "\n\n");


    for (int n = 2; n <= 100; n++)
    {
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauleg(0.f, 5.f, x, w, n);
        float c2 = 0.0;
        for(int i=1; i <= n; i++)
        {
            c2 += w[i]*f2b(x[i]);
        }
        fprintf(fp, "%4d %15g\n", n, fabs(c2a-c2));

        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }
    fprintf(fp, "\n\n");

    float c3a = 2.0/13.0;
    for (int n = 2; n <= 20; n++)
    {
        float* x = vector(1, n);
        float* w = vector(1, n);
        
        for(int i=1; i<=n; i++)
        {
            x[i] = 0.0;
            w[i] = 0.0;
        }

        gaulag(x, w, n, 0.0);
        float c3 = 0.0;
        for(int i=1; i <= n; i++)
        {
            c3 += w[i]*f3(x[i]);
        }
        fprintf(fp, "%4d %15f\n", n, fabs(c3a-c3));

        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }

    fclose(fp);
}