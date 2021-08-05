#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "nrutil.h"
//#include "nrutil.c"
//#include "ludcmp.c" 
//#include "lubksb.c"

// #include "/opt/NR/numerical_recipes.c/nrutil.h"
// #include "/opt/NR/numerical_recipes.c/nrutil.c"
// #include "/opt/NR/numerical_recipes.c/gaussj.c"


#define N 3 

void print_matrix(float **M, int size)
{
    for(int i=1; i<=size; ++i)
    {
        for(int j=1; j<=size; ++j)
        {
            printf("%10g ", M[i][j]);
        }
        printf("\n");
    }  
    printf("\n");
}

void print_matrix_to_file(float **M, int size, FILE *pt)
{
  for(int i=1; i<=size; ++i)
  {
        for(int j=1; j<=size; ++j)
        {
            fprintf(pt, "%15g ", M[i][j]);
        }
        fprintf(pt, "\n");
    }  
    fprintf(pt, "\n");
}

void copy_matrix(float **M, float **copy, int size)
{
    for(int i=1; i<=size; ++i)
    {
        for(int j=1; j<=size; ++j)
        {
            copy[i][j] = M[i][j];
        }
    }  
}

void product(float **X, float **Y, float **Z, int size)
{
    for (int i=1; i<=size; ++i)
    {
        for (int j=1; j<=size; ++j)
        {
            Z[i][j]= 0.;
            for (int k=1; k<=size; ++k)
                Z[i][j]+= X[i][k] * Y[k][j]; 
        }
    }
}  

//norma macierzy
float norm(float **M, int size)
{
    float max_m = M[1][1] >= 0 ? M[1][1] : -M[1][1];
    for(int i=1; i<=N; ++i)
    {
        for(int j=1; j<=N; ++j)
        {
            float temp_m = M[i][j]  >= 0 ? M[i][j] : -M[i][j];
            max_m = temp_m > max_m ? temp_m : max_m;
        }
    }
    return max_m;
}

void transpose(float **M)
{
    for (int i=1; i<=N-1; ++i)
    {
        for(int j=i+1; j<=N; ++j)
        {
            float temp = M[i][j];
            M[i][j] = M[j][i];
            M[j][i] = temp; 
        }
    }
}


int main(void)
{
	float **A, **B, **A_copy, **B_copy, **product_A, **product_B;

	//	Alokacja macierzy
	A = matrix(1, N, 1, N);
	B = matrix(1, N, 1, N);
    
    A_copy = matrix(1, N, 1, N);
    B_copy = matrix(1, N, 1, N);

    product_A = matrix(1, N, 1, N);
    product_B = matrix(1, N, 1, N);

    int *indx;
    float d;
    float f = 1.; // zmienna pomocnicza służąca do wypełnienia macierzy A i B
    for(int i=1; i<=N; ++i)
    {
        for(int j=1; j<=N; ++j)
        {
            A[i][j] = f;
            B[i][j] = f;
            f += 1.;
        }
    }
    B[1][1] = 1.1;
    
    FILE *fp = fopen("out.txt", "w");
    copy_matrix(A, A_copy, N);
    copy_matrix(B, B_copy, N);

    printf("A:\n");
    print_matrix(A, N);
    
    printf("B:\n");
    print_matrix(B, N);

    indx = ivector(1,N);


    ludcmp(A_copy, N, indx, &d);
    ludcmp(B_copy, N, indx, &d);

    printf("Rozkład LU macierzy A:\n");
    print_matrix(A_copy, N);

    printf("Rozkład LU macierzy B:\n");
    print_matrix(B_copy, N);

    float **a = matrix(1, N, 1, N); //wektor wyrazów wolnych dla A
    float **b = matrix(1, N, 1, N); //wektor wyrazów wolnych dla B

    for(int i=1; i<=N; ++i)
    {
        for(int j=1; j<=N; ++j)
        {
            a[i][j] = i==j ? 1. : 0.;
            b[i][j] = i==j ? 1. : 0.;
        }
    }

    
    for(int i=1; i<=N; i++){
        lubksb(A_copy, N, indx, a[i]);  
        lubksb(B_copy, N, indx, b[i]);
    }

    /*  Ponieważ lubksb zwraca kolumny odwróconej macierzy a nie wiersze,
        to, aby otrzymać macierze a i b, które są odwrotne do macierzy odpowiednio A i B,
        należy transponować otrzymane macierze a i b
    */ 
    transpose(a);
    transpose(b);

    printf("A^(-1):\n");
    print_matrix(a, N);

    printf("B^(-1):\n");
    print_matrix(b, N);

    float max_A = norm(A, N);
    float max_a = norm(a, N);
    
    float max_B = norm(B, N);
    float max_b = norm(b, N);
    

    fprintf(fp, "Norma macierzy A: %g \n", max_A);
    fprintf(fp, "Norma macierzy A^(-1): %g \n", max_a);
    fprintf(fp, "Wskaźnik uwarunkowania macierzy A (norma A * norma A^(-1)): %g\n\n", max_A * max_a);

    fprintf(fp, "Norma macierzy B: %g \n", max_B);
    fprintf(fp, "Norma macierzy B^(-1): %f \n", max_b);
    fprintf(fp, "Wskaźnik uwarunkowania macierzy B (norma B * norma B^(-1)): %f \n\n", max_B * max_b);

    fprintf(fp, "A*A^(-1):\n");
    product(A, a, product_A, N);
    print_matrix_to_file(product_A, N, fp);

    fprintf(fp, "B*B^(-1):\n");
    product(B, b, product_B, N);
    print_matrix_to_file(product_B, N, fp);
    
    //	Zwolnienie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);

    free_matrix(a, 1, N, 1, N);
    free_matrix(b, 1, N, 1, N);

    free_ivector(indx, 1, N);
    free_matrix(A_copy, 1, N, 1, N);
	free_matrix(B_copy, 1, N, 1, N);
    free_matrix(product_A, 1, N, 1, N);
	free_matrix(product_B, 1, N, 1, N);
    fclose(fp);
	return 0;
}
