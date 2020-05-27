#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/pythag.c"


float scalar_product(float *x1) //funkcja pozwalająca na otrzymanie normy euklidesowej
{
    float norm = 0;

    for (int j = 1; j <= n; j++)
    {
        norm += x1[j] * x1[j];
    }

    norm = sqrt(norm);

    return norm;
}

int main()
{
    float **A = matrix(1, n, 1, n); //MACIERZ A
    float **W = matrix(1, n, 1, n); //MACIERZ W (zad. 4)

    float *d = vector(1, n); //WEKTOR d
    float *e = vector(1, n); //WEKTOR e

    float *x0 = vector(1, n); //WEKTOR STARTOWY
    float *x1 = vector(1, n);

    float *lmbd = vector(1, n); //WEKTOR WARTOŚCI WŁASNYCH
    float *A_wartosci_wlasne = vector(1, n);

    int it = 8;        //8 iteracji
    float nomin = 0;   //licznik lambdy
    float denomin = 0; //mianownik lambdy
    float norm;        //norma euklidesowa

    int i, j;

    //ZADANIE 1
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
        {
            A[i][j] = sqrt(i + j);
            W[i][j] = sqrt(i + j); //TWORZYMY MACIERZ W JAKO NIENARUSZONĄ KOPIĘ A
        }

    printf("\n Macierz wejsciowa A:\n");

    //SPRAWDZAMY ZGODNOŚĆ MACIERZY WEJŚCIOWEJ Z ODPOWIEDZIAMI
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            printf("%8.5f", A[i][j]);
        }

        printf("\n");
    }

    //ZADANIE 2
    tred2(A, n, d, e);

    //ZADANIE 3
    tqli(d, e, n, A);

    printf("\nWartosci wlasne macierzy A, uzyskane przy pomocy biblioteki Numerical Recipes:\n");

    for (i = 1; i <= n; i++)
    {
        printf("%13.5e\n", d[i]);
    }

    printf("\n");

    //ZADANIE 4
    printf("\nKolejne wartosci lambda_1 we wszystkich osmiu iteracjach\n");

    for (int k = 1; k <= n; k++)
    {

        for (int i = 1; i <= n; i++)
        {
            x0[i] = 1;
        }

        //POCZATEK ITEROWANIA
        for (int i = 1; i <= it; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                x1[j] = 0;
                for (int k = 1; k <= n; k++)
                {
                    x1[j] = x1[j] + W[j][k] * x0[k];
                }
            }

            nomin = 0;
            denomin = 0;

            for (int j = 1; j <= n; j++)
            {
                nomin = nomin + x1[j] * x0[j];
                denomin = denomin + x0[j] * x0[j];
            }

            lmbd[k] = nomin / denomin;

            norm = 0;
            norm = scalar_product(x1);

            for (int j = 1; j <= n; j++)
                x1[j] = x1[j] / norm;

            for (int j = 1; j <= n; j++)
                x0[j] = x1[j];

            //test wartości lambda_1 we wszystkich ośmiu iteracjach
            if (k <= 1)
            {
                printf("lambda_1%d = %8.4f\n", i, lmbd[2]);
            }
        }

        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                W[i][j] = W[i][j] - x0[i] * lmbd[k] * x0[j];
            }
        }
    }

    printf("\nWartosci wlasne macierzy A uzyskane metoda potegowa:\n");

    for (i = 1; i <= n; i++)
    {
        printf("%13.5e\n", lmbd[i]);
    }

    printf("\n");

    return 0;
}