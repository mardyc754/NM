#include <iostream>
#include <cmath>
#include <stdio.h>

double f(double x)
{
    return log(x*x*x + 3*x*x + x + 0.1)*sin(18*x);
}

double simpson(double x[], double h, int N)
{
    double result = 0.0;
    for(int i=0; i<=(N/2-1); i++)
    {
        result += h/3 * (f(x[2*i]) + 4*f(x[2*i+1]) + f(x[2*i+2]));
    }
    return result;
}

double milne(double x[], double h, int N)
{
    double result = 0.0;
    for(int i=0; i<=(N/4-1); i++)
    {
        result += (4.0*h)/90.0 * (7*f(x[4*i]) + 32*f(x[4*i+1]) + 12*f(x[4*i+2]) + 32*f(x[4*i+3]) + 7*f(x[4*i+4]));
    }
    return result;
}

int main()
{
    const int n = 8;
    double **D = new double*[n+1];
    double a = 0.0, b = 1.0;
    FILE *out = fopen("out.txt", "w"); 
    for(int i=0; i<=n; i++)
    {
        D[i] = new double[n+1];
    }

    for(int w=0; w<=n; w++)
    {
        int N = pow(2, w+1);
        double h_w = (b-a)/N;
        double *x = new double[N+1];
        for(int j=0; j<=N; j++)
        {
            x[j] = a + j*h_w;
        }

        D[w][0] = simpson(x, h_w, N);
        for(int k=1; k<=w; k++)
        {
            D[w][k] = (pow(4, k)*D[w][k-1] - D[w-1][k-1])/(pow(4,k)-1);    
        }
        delete [] x;
        fprintf(out, "%d %15.10f %15.10f\n", w, D[w][0], D[w][w]);
    }

    fprintf(out, "\n\n");

    for(int w=0; w<=n; w++)
    {
        int N = pow(2, w+2);
        double h_w = (b-a)/N;
        double *x = new double[N+1];
        for(int j=0; j<=N; j++)
        {
            x[j] = a + j*h_w;
        }

        D[w][0] = milne(x, h_w, N);
        for(int k=1; k<=w; k++)
        {
            D[w][k] = (pow(4, k)*D[w][k-1] - D[w-1][k-1])/(pow(4,k)-1);    
        }
        delete [] x;
        fprintf(out, "%d %15.10f %15.10f\n", w, D[w][0], D[w][w]);
    }
    fprintf(out, "\n\n");

    for(int i=0; i<=n; i++)
    {
        delete [] D[i];
    }
    delete [] D;
    fclose(out); 
}
