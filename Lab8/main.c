#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"

typedef float (*fun)(float);

float f_1(float x)
{
    return 1/(1+x*x);
}

float f_2(float x)
{
    return cos(2*x);
}

void print_vector(float *v, int n)
{
    for(int i=1; i<=n; i++)
    {
        printf("%10g\n", v[i]);
    }
    printf("\n");
}

void print_matrix(float **M, int n, int m)
{
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=m; j++)
        {
            printf("%15g", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}


void wyzM(float *xw,float *yw, float *w, int n, float alfa, float beta)
{
    float **A = matrix(1,n,1,n);
    float **d = matrix(1,n,1,1);
    d[1][1] = alfa;
    d[n][1] = beta; 

    float h = 10.0/(n-1); // x_max = 5.0, x_min = -5.0

    for(int i=1; i<=n; i++){
        d[i][1] = 0.0;
        for(int j=1; j<=n; j++){
            A[i][j] = 0.0;
        }
    }
    A[1][1] = A[n][n] = 1.0;

    for(int i=2; i<=n-1; i++)
    {
            A[i][i] = 2.0;
            A[i][i+1] = h/(h+h); // lambda_i
            A[i][i-1] = 1.0 - A[i][i+1]; // mu_i
    }

    for(int i=2; i<=n-1; i++)
    {
        d[i][1] = 6.0/(2*h)*( (yw[i+1]-yw[i])/h - (yw[i]-yw[i-1])/h );
    }
    
    gaussj(A, n, d, 1);
    for(int i=1; i<=n; i++)
    {
        w[i] = d[i][1];
    }
    
    free_matrix(A,1,n,1,n);
    free_matrix(d,1,n,1,1); 
}

float wyzSx(float *xw,float *yw, float *m, int n, float x)
{
    int i = 0;
    float h = 10.0/(n-1.0);
    for(int j=2; j<=n; j++)
    {
        if(x >= xw[j-1] && x <= xw[j])
        {
            i = j;
            break;
        }
    }
    float A_i = ((yw[i] - yw[i-1])/h) - ((h/6.0) * (m[i]-m[i-1]));
    float B_i = yw[i-1] - (m[i-1] * ((h*h) / 6.0));
     
    return (m[i-1]*pow(xw[i] - x, 3)/(6.0*h)) + (m[i]*pow(x-xw[i-1], 3)/(6.0*h)) + A_i*(x- xw[i-1]) + B_i;
}

float d2x(float x, float dx, fun f)
{
    return (f(x-dx)-2.0*f(x)+f(x+dx))/(dx*dx);
}

void interpolation(int n, fun f, FILE *fp)
{
    float *xw = vector(1, n);
    float *yw = vector(1,n);
    float *m = vector(1, n);

    xw[1] = -5.0;
    float h = 10.0/(n-1);
    for(int i=2; i<=n; i++)
    {
        xw[i] = xw[i-1] + h; 
    }

    for(int i=1; i<=n; i++)
    {
        yw[i] = f(xw[i]); 
    }
    
    for(int i=1; i<=n; i++)
    {
        m[i] = 0.0; 
    } 

    wyzM(xw, yw, m, n, 0.0, 0.0);
    if(n==10)
    {
        for(int i=1; i<=n; i++)
        {
            fprintf(fp, "%7.2g %15g %15g\n", 
                xw[i], m[i], d2x(xw[i], 0.01, f));
        }
    } 
    else 
    {
        for(float x=-5.0; x<=5.0; x+= 0.01)
        {
            fprintf(fp, "%7.3g %15g\n", x, wyzSx(xw, yw, m, n, x));
        }
    }
    fprintf(fp, "\n\n");
    free_vector(xw, 1, n);
    free_vector(yw, 1, n);
    free_vector(m, 1, n);
}

int main(void)
{
    FILE *fp1 = fopen("f1.dat", "w");
    interpolation(5, &f_1, fp1);
    interpolation(8, &f_1, fp1);
    interpolation(21, &f_1, fp1);
    fclose(fp1);

    FILE *fp2 = fopen("f2.dat", "w");
    interpolation(5, &f_2, fp2);
    interpolation(8, &f_2, fp2);
    interpolation(21, &f_2, fp2);
    fclose(fp2);

    FILE *fp3 = fopen("pochodne.dat", "w");
    interpolation(10, &f_1, fp3);
    fclose(fp3);
}
