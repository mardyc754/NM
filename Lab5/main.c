#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7

void print_matrix(double (*A)[N])
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%15g", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(double *v)
{
    for(int i=0; i<N; i++)
    {
        printf("%15g", v[i]);
    }
    printf("\n\n");
}

void copy_matrix(double (*A)[N], double (*W)[N])
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            W[i][j] = A[i][j];
        }
    }
}

void copy_column_to_vector(double (*M)[N], double *V, int col)
{
    for(int i=0; i<N; i++)
    {
        V[i] = M[i][col];
    }
}

void mm_product(double (*X)[N], double (*Y)[N], double (*Z)[N])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            Z[i][j]= 0.;
            for (int k=0; k<N; k++)
                Z[i][j]+= X[i][k] * Y[k][j]; 
        }
    }
}  


void mv_product(double (*M)[N], double (*V)[N], int col, double *R)
{
    for (int i=0; i<N; i++)
    {
        R[i] = 0;
        for (int j=0; j<N; j++)
        {
            R[i] += M[i][j] * V[j][col]; 
        }
    }
}

double vv_scalar_product(double *A, double *B)
{
    double product = 0;
    for(int i=0; i<N; i++)
    {
        product += A[i] * B[i]; 
    }
    return product;
}

void vv_tensor_product(double *A, double *B, double (*M)[N])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            M[i][j] = A[i] * B[j];
        }
    }
}

void product_by_scalar(double (*M)[N], double num)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            M[i][j] *= num;
        }
    }
}

void initialize_vector(double (*M)[N], int col)
{
    for(int i=0; i<N; i++)
    {
        M[i][col] = 1.;
    }
}

double norm(double *X)
{
    double sum_of_squares = 0.;
    for(int i=0; i<N; i++)
    {
        sum_of_squares += X[i] * X[i];
    }
    return sqrt(sum_of_squares);
}

void transpose(double (*M)[N], double (*M_tr)[N])
{
    for (int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            M_tr[i][j] = M[j][i];
        }
    }
}

void matrix_substract(double (*A)[N], double (*B)[N])
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            A[i][j] -= B[i][j];
        }
    }
}

int main(void)
{
    double A[N][N] = {0};
    double W[N][N] = {0};
    double X[N][N] = {0};
    double X_col[N] = {0};
    double X_tr[N][N] = {0};
    double tensor_pr[N][N] = {0};
    double lambda[N] = {0};
    double temp[N] = {0};
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            A[i][j] = 1./sqrt(2.+fabs(i-j));
        }
    }
    FILE *fp = fopen("lambda.dat", "w");

    copy_matrix(A, W);

    for(int k=0; k<N; k++)
    {
        initialize_vector(X, k); 
        for(int i=1; i<=12; i++)
        {    
            mv_product(W, X, k, temp);
            copy_column_to_vector(X, X_col, k);
            
            lambda[k] = vv_scalar_product(temp, X_col)/vv_scalar_product(X_col, X_col);
            
            fprintf(fp, "%d %g\n", i, lambda[k]);

            for(int j=0; j<N; j++)
            {
                X[j][k] = temp[j] / norm(temp);
            }
            
        }
        copy_column_to_vector(X, X_col, k);
        vv_tensor_product(X_col, X_col, tensor_pr);
            
        product_by_scalar(tensor_pr, lambda[k]);
            
        matrix_substract(W, tensor_pr);
            
        fprintf(fp, "\n\n");
    }
	
	printf("Macierz wektorów własnych X:\n");
    print_matrix(X);
	
    transpose(X, X_tr);
    
    double D[N][N] = {0};
    double temp_matrix[N][N] = {0};

    mm_product(X_tr, A, temp_matrix);
    mm_product(temp_matrix, X, D);

    fclose(fp);

    FILE *fp1 = fopen("macierzD.dat", "w");
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            fprintf(fp1, "%15g", D[i][j]);
        }
        fprintf(fp1, "\n");
    }
    fclose(fp1);
	return 0;
}

