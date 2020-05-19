#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/gaussj.c"

#define N 400 // rozmiar macierzy M: NxN

int main(void)
{
    double omega = 1;
    double h = 0.1;
    double vpocz = 0;
    double A = 1;
    double dt;

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

    M[2][1] = -1;

    for(int i = 1; i <= N; i++){
        M[i][i] = 1;
    	dt = h;
        
        if(i > 2){
            M[i][i-1] = ((omega * omega) * (dt * dt) - 2);
            M[i][i-2] = 1;
        }
        
        b[1][1] = A;
        b[2][1] = vpocz * dt;
    }

	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	for (int i = 1; i <= N; ++i){
		printf("%f ", dt*(i-1));
		printf("%g\n", b[i][1]);
	}
	
	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}

