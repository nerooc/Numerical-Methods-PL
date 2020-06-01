#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"

#define x0 2
#define frand() ((double)rand()) / (RAND_MAX + 1.0)
#define x_min -4
#define x_max 4
#define Crand (frand() - 0.5) / 5.;
#define sigma (x_max - x_min) / 16.
#define PI 3.14159265359

double function(double x)
{
    double r = 1;

    r *= sin((14 * PI * x) / (x_max - x_min));
    r *= (exp(-1 * ((x - x0) * (x - x0)) / (2 * sigma * sigma)) + exp(-1 * ((x + x0) * (x + x0)) / (2 * sigma * sigma)));

    return r;
}

int main(void)
{

    int N = 201;
    int m_test = 50;

    float alpha = 0.;
    float beta = 0.;

    double alpha_licznik, alpha_mianownik, beta_licznik, beta_mianownik;

    FILE *file1;
    FILE *file2;

    float *x = vector(1, N);
    float *y = vector(1, N);
    float *F = vector(1, N);

    for (int i = 1; i <= N; i++)
    {
        x[i] = x_min + ((x_max - x_min) / (N - 1.)) * (i - 1.);
    }

    for (int i = 1; i <= N; i++)
    {
        y[i] = function(x[i]);
    }

    //for (int i = 1; i <= N; ++i)
    //    printf("%d %f %f\n", i, x[i], y[i]);

    float **phi = matrix(1, N, -1, m_test);

    for (int i = 1; i <= N; i++)
    {
        phi[i][-1] = 0.;
        phi[i][0] = 1.;
    }

    for (int k = 1; k <= m_test; k++)
    {
        alpha_licznik = 0.;
        alpha_mianownik = 0.;

        beta_licznik = 0.;
        beta_mianownik = 0.;

        for (int i = 1; i <= N; i++)
        {

            alpha_licznik += x[i] * phi[i][k - 1] * phi[i][k - 1];
            alpha_mianownik += phi[i][k - 1] * phi[i][k - 1];

            beta_licznik += x[i] * phi[i][k - 1] * phi[i][k - 2];
            beta_mianownik += phi[i][k - 2] * phi[i][k - 2];

            //printf("al= %f am = %f\n", alpha_licznik, alpha_mianownik);
            //printf("bl= %f bm = %f\n", beta_licznik, beta_mianownik);
        }

        alpha = alpha_licznik / alpha_mianownik;

        if (k == 1)
        {
            beta = 0;
        }
        else
        {
            beta = beta_licznik / beta_mianownik;
        }

        //printf("alfa = %f beta = %f\n", alpha, beta);

        for (int i = 1; i <= N; i++)
        {
            phi[i][k] = (x[i] - alpha) * (phi[i][k - 1]) - beta * phi[i][k - 2];
        }
    }

    file1 = fopen("Gram.dat", "w");

    for (int i = 1; i <= N; i++)
        fprintf(file1, " %f %f %f %f %f %f %f %f\n", x[i], phi[i][0] / phi[1][0], phi[i][1] / phi[1][1], phi[i][2] / phi[1][2], phi[i][3] / phi[1][3], phi[i][4] / phi[1][4], phi[i][5] / phi[1][5], phi[i][6] / phi[1][6]);

    fclose(file1);

    file1 = fopen("pkt.dat", "w");
    file2 = fopen("approx.dat", "w");

    float cj = 0.;
    float sj = 0.;

    for (int m = 10; m <= 50; m += 20)
    {

        for (int k = 1; k <= N; k++)
        {

            F[k] = 0.;

            for (int j = 0; j <= m; j++)
            {
                cj = 0.;
                sj = 0.;

                for (int i = 1; i <= N; i++)
                {
                    cj += y[i] * phi[i][j];
                    sj += phi[i][j] * phi[i][j];

                    //printf("sj = %f\n", sj);
                    // printf("cj= %f\n", cj);
                }

                F[k] += ((cj / sj) * phi[k][j]);
            }

            //printf("%f %f\n", x[k], F[k]);

            fprintf(file1, "%f %f\n", x[k], y[k]);
            fprintf(file2, "%f %f\n", x[k], F[k]);
        }

        fprintf(file1, "\n\n");
        fprintf(file2, "\n\n");
    }

    fclose(file1);
    fclose(file2);

    free_vector(x, 1, N);
    free_vector(y, 1, N);
    free_vector(F, 1, N);
    free_matrix(phi, 1, N, -1, m_test);

    return 0;
}
