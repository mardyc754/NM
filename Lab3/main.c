#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 2000


double sum_of_squares(double *a)
{
    double sum = 0.;
    for(int i=0; i<=N; i++)
    {
        sum += a[i]*a[i];
    }
    return sum;
}

void iterate_jacobi(double beta, double F_0, FILE *fp)
{
	printf("beta=%.1lf, F_0=%.1lf\n", beta, F_0);

	double b[N+1] = {0.0};

    double v_0 = 0., w = 1.;
    double h = 0.02;
    double Omega = 0.8;

    b[0] = 1.;
    b[1] = 0.;
    double f = 1.0;
    for(int i=2; i<=N; i++)
    {

        b[i] = F_0 * sin(Omega*h*f)*h*h;
        f += 1.0;
    }

    double a_1 = 1.;
    double a_2 = w*w*h*h - 2. - beta*h;
    double a_3 = 1. + beta*h;

    double d_0[N+1] = {0.0};
    double d_1[N+1] = {0.0};
    double d_2[N+1] = {0.0};

    printf("a_1=%f a_2=%f a_3=%f\n", a_1, a_2, a_3);
    d_0[0] = d_0[1] = 1.0;
    d_1[1] = -1.0;
    for(int i=2; i<=N; i++)
    {
        d_0[i] = a_3;
        d_1[i] = a_2;
        d_2[i] = a_1;
    }

    double x_n[N+1] = {1.0};
    double x_s[N+1] = {1.0};
    x_n[0] = x_s[0] = 1.0;
    x_n[1] = x_n[0] - v_0*h;
    x_s[1] = x_s[0] - v_0*h;
    
    int it = 0;
    while (it<3000)
    {

        for(int i=2; i<=N; i++)
        {
            x_n[i] = (1/d_0[i]) * (b[i] - d_1[i]*x_s[i-1] - d_2[i]*x_s[i-2]);
        }
        
        it++;
        
        if(fabs(sum_of_squares(x_n) - sum_of_squares(x_s)) < 1e-6)
        {
            break;
        }

        for(int i=0; i<=N; i++)
        {   
            x_s[i] = x_n[i];
    }
    }

    printf("Liczba iteracji: %d\n ", it);
 

    for(int i=0; i<=N; i++)
    {
        fprintf(fp, "%5g %15lf\n", h, x_n[i]);
        h += 0.02;
    }
    fprintf(fp, "\n\n");

    for(int i=0; i<=5; i++)
    {
        printf("%f %g\n", b[i], b[N-i]);
    }
    printf("\n");

}

int main(void)
{

    FILE *fp = fopen("out.dat", "w");
    iterate_jacobi(0.0, 0.0, fp);	// beta = 0.4, F_0 = 0.0
    iterate_jacobi(0.4, 0.0, fp);	// beta = 0.4, F_0 = 0.0
    iterate_jacobi(0.4, 0.1, fp);	// beta = 0.4, F_0 = 0.1
    fclose(fp);

	return 0;
}
