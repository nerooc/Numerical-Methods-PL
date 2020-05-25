#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/tqli.c"
#include "numerical_recipes.c/tred2.c"
#include "numerical_recipes.c/pythag.c"

#define nx 20
#define ny 20
#define n nx *ny

int main()
{

    int m = 10;
    double t = -0.021;

    float **H, **Y, **X;
    float *d, *e;

    H = matrix(1, n, 1, n);
    Y = matrix(1, n, 1, n);
    X = matrix(1, n, 1, n);

    d = vector(1, n);
    e = vector(1, n);

    //zmienne potrzebne do wypełnienia macierzy H
    int i = 0;
    int j = 0;
    int l = 0;
    int k = 0;

    //WYPEŁNIAMY Z
    for (i = 1; i <= nx; i++)
    {
        for (j = 1; j <= ny; j++)
        {
            l = j + (i - 1) * ny;
            for (k = 1; k <= n; k++)
                H[l][k] = 0.;
            if (i > 1)
                H[l][l - ny] = t; //dla i=1 nie ma sasiada z lewej strony
            if (i < nx)
                H[l][l + ny] = t; //dla i=nx nie ma sasiada z prawej strony
            H[l][l] = -4 * t;
            if (j > 1)
                H[l][l - 1] = t; //dla j=1 nie ma sasiada ponizej siatki
            if (j < ny)
                H[l][l + 1] = t; //dla j=ny nie ma sasiada powyzej siatki
        }
    }

    //WYPEŁNIAMY Y JAKO MACIERZ JEDNOSTKOWĄ
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
            Y[i][j] = (float)(i == j);
    }

    tred2(H, n, d, e);
    tqli(d, e, n, Y);

    //MNOŻENIE MACIERZY
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            float temp = 0.;
            for (k = 1; k <= n; k++)
            {
                temp += H[i][k] * Y[k][j];
            }
            X[i][j] = temp;
        }
    }

    int *indx = ivector(1, n);
    float e1 = 0.;
    float e2 = 0.;
    float l1 = 0.;
    float l2 = 0.;

    for (l = 1; l <= n; l++)
        indx[l] = l; // inicjalizacja
    for (l = 1; l <= n - 1; l++)
    {
        for (k = n; k >= l + 1; k--)
        {
            e1 = d[k - 1];
            e2 = d[k];
            l1 = indx[k - 1];
            l2 = indx[k];
            if (e2 < e1)
            { //wymieniamy energie i indeksy wektorów miejscami
                d[k] = e1;
                d[k - 1] = e2;
                indx[k] = l1;
                indx[k - 1] = l2;
            }
        }
    }

    FILE *fp;
    fp = fopen("dane.dat", "w");
    for (i = 1; i <= nx; i++)
    {
        for (j = 1; j <= ny; j++)
        {
            l = j + (i - 1) * ny;
            fprintf(fp, "%6d %6d ", i, j);
            for (k = 1; k <= m; k++)
                fprintf(fp, " %12.6f ", X[l][indx[k]]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    FILE *fc;
    fc = fopen("danewarwl.dat", "w");
    for (i = 1; i <= nx; i++)
    {
        fprintf(fc, " %12.6f ", d[i]);
    }
    fclose(fc);
}