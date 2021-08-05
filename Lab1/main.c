#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Dyrektywy zakladajace, ze te trzy pliki sa skopiowane do aktualnego katalogu. */
#include "nrutil.h"
//#include "nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
//#include "gaussj.c" // To tez mozna usunac, jesli plik jest dodany w poleceniu kompilacji.

/* Dyrektywy dla Taurusa (nie wymagaja kopiowania plikow, ale Taurus musi dzialac...) */
// #include "/opt/NR/numerical_recipes.c/nrutil.h"
// #include "/opt/NR/numerical_recipes.c/nrutil.c"
// #include "/opt/NR/numerical_recipes.c/gaussj.c"

#define N 200 // rozmiar macierzy M: NxN

int main(void)
{
    float A = 1.;
    float h = 0.1;
    float w = 1.;

	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);

	// 	Wypelnienie macierzy M i wektora b
	for (int i = 1; i <= N; ++i)
	{
		b[i][1] = 0.0;
		for (int j = 1; j <= N; ++j)
			M[i][j] = 0.0;
	}

    for (int i=1; i<=N; ++i)
    {
        for(int j=1; j<=N; ++j)
        {
            if(i==j || i==(j+2))
            {
                M[i][j] = 1.0;
            }
            else if(i==(j+1))
            {
                M[i][j] = w*w*h*h-2.;
            }
            else 
            {
                M[i][j] = 0.0;
            }
        }  
    }
    M[2][1] = -1; // bo x_1 jest rowne x_0
    b[1][1] = A;
	
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	for (int i = 1; i <= N; ++i)
	{
		printf("%g %g\n", h, b[i][1]);
        h += 0.1;
    }
	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}
