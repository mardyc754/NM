#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 5

double licz_r(double a[N+1], double b[N+1], int n, double x_j){
    b[n] = 0;
    for(int k=n-1; k>=0; k--){
        b[k] = a[k+1] + x_j*b[k+1];
    }
    return a[0] + x_j*b[0];
}

int main(void)
{
    double a[N+1] = {240, -196, -92, 33, 14, 1};
    double b[N+1] = {0};
    int n;
    double c[N] = {0};
    double x0, x1, Rj, Rj_pr;
    FILE *fp = fopen("out.txt", "w"); 
    for(int L=1; L<=N; L++){
        n = N - L + 1;
        x0 = 0;
        for(int it=1; it<=30; it++){
            Rj = licz_r(a, b, n, x0);
            Rj_pr = licz_r(b, c, n-1, x0);
            x1 = x0 - Rj/Rj_pr;

            fprintf(fp, "%5d %5d %15g %15g %15g\n", L, it, x1, Rj, Rj_pr);
            if(fabs(x1-x0)<1.0e-7) break;
            
            x0 = x1;
            
        }
        fprintf(fp, "\n");
        for(int i=0; i<=(n-1); i++) a[i] = b[i];
    } 
    fclose(fp);
	return 0;
}

