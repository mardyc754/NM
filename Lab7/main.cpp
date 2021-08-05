#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <functional>

double f(double x)
{
    return 1./(1.+x*x);
}

// równoodległe położenie węzłów
void create_nodes_1(double *x_m, double x_min, double x_max, int n)
{
    double h = (x_max - x_min)/n;
    x_m[0] = x_min;
    for(int i=1; i<=n; i++)
    {
        x_m[i] = x_m[i-1] + h;
    }
}

// miejsca zerowe wielomianu Czebyszewa
void create_nodes_2(double *x_m, double x_min, double x_max, int n)
{
    for(int i=0; i<=n; i++)
    {
        x_m[i]=0;
    }
    for(int i=0; i<=n; i++)
    {
        x_m[i] = 0.5 * ((x_min-x_max)*cos(M_PI*(2*i+1)/(2*n+2))+ (x_min + x_max));
    }
}

double W_x(double x, int m, double *x_m, double **f_m)
{
    double result = 0;
    for(int j=0; j<=m; j++)
    {
        double temp = f_m[j][j];
        for(int i=0; i<=j-1; i++)
        {
            temp *= (x - x_m[i]);
        }
        result += temp;
    }
    return result;
}

void print_vals(double *x, double *y, int n)
{
    for(int i=0; i<=n; i++)
    {
        printf("%2d %15.1f %15g\n", i, x[i], y[i]);
    }
    printf("\n");
}

void print_matrix(double **M, int n)
{
    for(int i=0; i<=n; i++)
    {
        for(int j=0; j<=n; j++)
            printf("%15g", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

void interpolation(int n, std::function<void(double *, double, double, int)> nodes, FILE *fp)
{
    double *xm = new double[n+1];

    nodes(xm, -5.0, 5.0, n);
    double *ym = new double[n+1];
    for(int i=0; i<=n; i++)
    {
        ym[i] = f(xm[i]);
    }

    double **fm = new double *[n+1];
    for(int i=0; i<=n; i++)
    {
        fm[i] = new double [n+1];
    }
    for(int i=0; i<=n; i++)
    {
        for(int j=0; j<=n; j++)
            fm[i][j] = j==0 ? ym[i] : 0;
    }

    for(int j=1; j<=n; j++)
    {
        for(int i=j; i<=n; i++)
        {
            fm[i][j] = (fm[i][j-1]-fm[i-1][j-1])/(xm[i]-xm[i-j]);
        }
    }
    
    for(double x=-5.0; x<=5.0; x+= 0.01)
    {
        fprintf(fp, "%6.2f %15g\n", x, W_x(x, n, xm, fm));
    }
    fprintf(fp, "\n\n");
    delete [] xm;
    delete [] ym;
    for(int i=0; i<=n; i++)
    {
        delete [] fm[i];
    }
    delete [] fm;
}

int main(void)
{
    FILE *fp1 = fopen("zad_1.dat", "w");
    for(int i=5; i<=20; i+=5)
    {
        interpolation(i, create_nodes_1,fp1);
    }
    fclose(fp1);

    FILE *fp2 = fopen("zad_2.dat", "w");
    for(int i=5; i<=20; i+=5)
    {
        interpolation(i,  create_nodes_2, fp2);
    }
    fclose(fp2);
    return 0;
}

