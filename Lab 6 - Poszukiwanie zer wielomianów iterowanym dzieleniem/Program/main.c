#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"

#define N 5


double licz_r(double *, double *, unsigned, double);

double licz_r(double *a, double *b, unsigned n, double xj)
{
    //na podstawie wzorow (3), (4) i (5)
    b[n] = 0;

    for (int k = n - 1; k >= 0; k--)
    {
        b[k] = a[k + 1] + xj * b[k + 1];
    }
    return a[0] + xj * b[0];
}

int main(void)
{

    int IT_MAX = 30; //maksymalna ilosc iteracji
    int n;           //aktualny stopien wielomianu

    double x0;
    double x1;

    double Rj;
    double RjPrim;

    double *a = dvector(0, N); //6-elementowe wektory
    double *b = dvector(0, N);
    double *c = dvector(0, N);

    a[0] = 240.;
    a[1] = -196.;
    a[2] = -92.;
    a[3] = 33.;
    a[4] = 14.;
    a[5] = 1.;

    FILE *file; //plik do ktorego bedziemy zapisywac

    file = fopen("results.dat", "w"); //otwieramy plik do zapisu

    for (int L = 1; L <= N; L++)
    {
        x0 = 0.;
        n = N - L + 1;
        fprintf(file, "%5s  %5s  %5s  %2s  %15s\n", "L", "it", "x0", "Rj", "Rj'");

        for (int it = 1; it <= IT_MAX; it++)
        {
            Rj = licz_r(a, b, n, x0);
            RjPrim = licz_r(b, c, n - 1, x0);

            x1 = x0 - Rj / RjPrim;

            fprintf(file, "%5.d  %5.d  %5.lf  %5.5e  %5.lf\n", L, it, x0, Rj, RjPrim);

            if (fabs(x1 - x0) < 1.e-7)
            {
                break;
            }

            x0 = x1;
        }

        for (int i = 0; i <= (n - 1); i++)
            a[i] = b[i];

        fprintf(file, "\n%d miejsce zerowe: x = %.0f\n\n", L, x1);
    }

    fclose(file); //zamykamy plik

    //zwalniamy pamiec wektorow
    free_dvector(a, 0, N);
    free_dvector(b, 0, N);
    free_dvector(c, 0, N);

    return 0;
}