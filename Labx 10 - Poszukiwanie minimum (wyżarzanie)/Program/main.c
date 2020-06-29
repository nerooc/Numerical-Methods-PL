#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"

#define PI 3.14159265359
#define N 200


double d_rand(const double min, const double max)
{
    double r = (double)rand() / RAND_MAX; // Przedzial [0 , 1]
    r = r * (max - min) + min;            // Przeskalowanie do [min , max]
    return r;
}

double f(double x, double y)
{
    double r = 1;

    r *= sin(x) * sin(y);
    r -= (exp(-1 * ((x + (PI / 2)) * (x + (PI / 2))) - ((y - (PI / 2)) * (y - (PI / 2)))));

    return r;
}

int main(void)
{
    int iT = 0;
    int k = 0;
    int i = 0;
    float T = 0.;

    float MIN = 1.;
    float x_dla_MIN = 0.;
    float y_dla_MIN = 0.;
    int index = 0;

    float x_min = -10.;
    float x_max = 10.;

    float y_min = -10.;
    float y_max = 10.;

    float *x = vector(1, N);
    float *y = vector(1, N);
    float *deltaX = vector(1, N);
    float *deltaY = vector(1, N);

    FILE *file1;
    FILE *file2;

    for (int j = 1; j <= N; j = j + 1)
    {
        x[j] = 5;
        y[j] = 5;
    }

    file1 = fopen("T.dat", "w");
    file2 = fopen("w0.dat", "w");

    for (iT = 0; iT <= 20; iT = iT + 1)
    {
        T = 1;

        //printf("iteracja numer %d, T = %f\n\n", iT, T);

        for (k = 1; k < 100; k = k + 1)
        {
            fprintf(file2, "%f\n", f(x[1], y[1]));

            for (i = 1; i <= N; i = i + 1)
            {

                deltaX[i] = d_rand(-1, 1);
                deltaY[i] = d_rand(-1, 1);

                if (x[i] + deltaX[i] > x_max)
                {
                    deltaX[i] = x_max - x[i];
                }
                else if (x[i] + deltaX[i] < x_min)
                {
                    deltaX[i] = x_min - x[i];
                }
                else
                {
                    deltaX[i] = d_rand(-1, 1);
                }

                if (y[i] + deltaY[i] > y_max)
                {
                    deltaY[i] = y_max - y[i];
                }
                else if (y[i] + deltaY[i] < y_min)
                {
                    deltaY[i] = y_min - y[i];
                }
                else
                {
                    deltaY[i] = d_rand(-1, 1);
                }

                if (f(x[i] + deltaX[i], y[i] + deltaY[i]) < f(x[i], y[i]))
                {
                    x[i] = x[i] + deltaX[i];
                    y[i] = y[i] + deltaY[i];
                }
                else if (d_rand(0, 1) < exp(-1 * ((f(x[i] + deltaX[i], y[i] + deltaY[i]) - f(x[i], y[i])) / T)))
                {
                    x[i] = x[i] + deltaX[i];
                    y[i] = y[i] + deltaY[i];
                }
            }
        }

        if (iT == 0)
        {
            for (int i = 1; i <= N; i++)
            {
                fprintf(file1, "%f %f\n", x[i], y[i]);
            }
            fprintf(file1, "\n\n");
        }

        if (iT == 7)
        {
            for (int i = 1; i <= N; i++)
            {
                fprintf(file1, "%f %f\n", x[i], y[i]);
            }
            fprintf(file1, "\n\n");
        }

        if (iT == 20)
        {
            for (int i = 1; i <= N; i++)
            {
                fprintf(file1, "%f %f\n", x[i], y[i]);
            }
            fprintf(file1, "\n\n");
        }
    }

    for (int i = 1; i <= N; i++)
    {
        if (f(x[i], y[i]) < MIN)
        {
            MIN = f(x[i], y[i]);
            x_dla_MIN = x[i];
            y_dla_MIN = y[i];
            index = i;
        }
    }

    printf("\nMinimum funkcji wynosi f(x[%d], y[%d]) = %f, dla x[%d] = %f, y[%d] = %f\n\n", index, index, MIN, index, x_dla_MIN, index, y_dla_MIN);

    free_vector(x, 1, N);
    free_vector(y, 1, N);
    free_vector(deltaX, 1, N);
    free_vector(deltaY, 1, N);
}
