#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "numerical_recipes.c/nrutil.h"
#include "numerical_recipes.c/nrutil.c"
#include "numerical_recipes.c/gaussj.c"

float f1(float x);
float f2(float x);
double deriv(double x);
void secDeriv();
void wyzM(float *xw, float *yw, float *m, int n, float alfa, float beta);
float wyzSx(float *xw, float *yw, float *m, int n, float x);

//funkcja f1
float f1(float x)
{
    return 1. / (1 + x * x);
}

//funkcja f2
float f2(float x)
{
    return cos(2 * x);
}

//funkcja zwracajaca "dokladniejsza" wartosc drugiej pochodnej
double deriv(double x)
{
    double value = 0.0;
    value = (f1(x - 0.01) - (2 * f1(x)) + (f1(x + 0.01))) / pow(0.01, 2);
    return value;
}

//funkcja powtarzajaca nasze dzialania dla n = 10
void secDeriv()
{
    const int n = 10;

    const float x_max = 5;
    const float x_min = -5;

    float alfa = 0.;
    float beta = 0.;

    //otwieramy plik do ktorego bedziemy zapisywac
    FILE *file = fopen("pochodne.dat", "w");

    //wektor xw i jego wypelnienie
    float *xw = vector(1, n);

    for (int i = 1; i <= n; i++)
    {
        float h = (x_max - x_min) / ((float)(n)-1);
        xw[i] = x_min + h * (i - 1);
    }

    //wektor yw i jego wypelnienie
    float *yw = vector(1, n);

    for (int i = 1; i <= n; i++)
    {
        yw[i] = f1(xw[i]);
    }

    //wektor m
    float *m = vector(1, n);

    //wywolanie funkcji wyzM
    wyzM(xw, yw, m, n, alfa, beta);

    //wypisanie wartosci w odpowiedni sposob
    for (int i = 1; i <= n; i++)
    {
        fprintf(file, "%lf %lf %lf\n", xw[i], m[i], deriv(xw[i]));
    }

    //zamykamy plik
    fclose(file);

    free_vector(xw, 1, n);
    free_vector(yw, 1, n);
    free_vector(m, 1, n);
}

//funkcja wyzM wyznaczajaca wartosci drugich pochodnych w wezlach
void wyzM(float *xw, float *yw, float *m, int n, float alfa, float beta)
{

    float **A = matrix(1, n, 1, n);

    //wypelniamy macierz zerami i dwojkami na przekatnej
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
        {
            if (i != j)
                A[i][j] = 0.;
            else
                A[i][j] = 2.;
        }

    //pierwszy i ostatni element rowne 1
    A[1][1] = 1.;
    A[n][n] = 1.;

    //wektor wyrazow wolnych jako macierz
    float **d = matrix(1, n, 1, n);

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            d[i][j] = 0;
        }
    }

    //ustawiamy poczatek i koniec jako alfa i beta
    d[1][1] = alfa;
    d[n][1] = beta;

    //tworzymy wektory h, lambda i mikro (troche za pozno zorientowalem sie ze w sumie to nie trzeba az tak)
    float *h = vector(1, n);

    for (int i = 1; i <= n; i++)
    {
        h[i] = 0;
    }

    float *lambda = vector(1, n);

    for (int i = 1; i <= n; i++)
    {
        lambda[i] = 0;
    }

    float *mikro = vector(1, n);

    for (int i = 1; i <= n; i++)
    {
        mikro[i] = 0;
    }

    //testujemy wartosci
    for (int i = 1; i <= n; i++)
    {
        //printf("xw[%d] = %f\n", i, xw[i]);
    }

    //printf("\n\n");

    for (int i = 2; i <= n; i++)
    {
        h[i] = xw[i] - xw[i - 1];
        //printf(" h[%d] = %f\n", i, h[i]);
    }
    //printf("\n\n");

    for (int i = 2; i < n; i++)
    {
        lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
        //printf(" lmbd[%d] = %f\n", i, lambda[i]);
    }
    //printf("\n\n");

    for (int i = 2; i < n; i++)
    {
        mikro[i] = 1. - lambda[i];
        //printf(" mikro[%d] = %f\n", i, mikro[i]);
    }
    //printf("\n\n");

    //wypelniamy wektor d
    for (int i = 2; i < n; ++i)
        d[i][1] = (6. / (h[i] + h[i + 1])) * (((yw[i + 1] - yw[i]) / (h[i + 1])) - (((yw[i]) - yw[i - 1]) / (h[i])));

    //sprawdzenie
    /*
    for (int i = 1; i <= n; ++i)
    {
        printf("%g\n", d[i][1]);
    }
    printf("\n\n");
    */

    //wypelniamy macierz A
    for (int i = 2; i < n; i++)
    {
        A[i][i + 1] = lambda[i];
        A[i][i - 1] = mikro[i];
    }

    //sprawdzamy macierz A
    /*
    for (int i = 1; i <= n; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            printf("%f ", A[i][j]);
        }

        printf("\n");
    }

    printf("\n\n");
    */

    //korzystamy z funkcji gaussj do obliczenia rownania
    gaussj(A, n, d, 1);

    //wypisujemy wektor d i przepisujemy
    for (int i = 1; i <= n; ++i)
    {
        //printf("%g ", d[i][1]);
        m[i] = d[i][1];
    }

    free_matrix(A, 1, n, 1, n);
    free_matrix(d, 1, n, 1, n);
    free_vector(h, 1, n);
    free_vector(lambda, 1, n);
    free_vector(mikro, 1, n);
}

//funkcja sluzaca do wyznaczania wartosci funkcji w polozeniu miedzywezlowym
float wyzSx(float *xw, float *yw, float *m, int n, float x)
{
    int i;
    float sx;

    for (int c = 1; c < n; c++)
    {
        if (x >= xw[c] && x <= xw[c + 1])
            i = c;
    }

    float Ai = ((yw[i + 1] - yw[i]) / (xw[i + 1] - xw[i])) - ((xw[i + 1] - xw[i]) / 6) * (m[i + 1] - m[i]);
    float Bi = yw[i] - (m[i] * (((xw[i + 1] - xw[i]) * (xw[i + 1] - xw[i])) / 6));

    sx += m[i] * (((xw[i + 1] - x) * (xw[i + 1] - x) * (xw[i + 1] - x)) / (6 * (xw[i + 1] - xw[i]))) + m[i + 1] * (((x - xw[i]) * (x - xw[i]) * (x - xw[i])) / (6 * (xw[i + 1] - xw[i]))) + Ai * (x - xw[i]) + Bi;

    return sx;
}

int main(void)
{

    //robimy to dla f1 i f2
    for (int funkcja = 1; funkcja <= 2; funkcja++)
    {
        //granice
        const float x_max = 5;
        const float x_min = -5;

        float alfa = 0.;
        float beta = 0.;

        FILE *file;

        //otwieramy plik do zapisu
        if (funkcja == 1)
        {
            file = fopen("f1.dat", "w");
        }
        else
        {
            file = fopen("f2.dat", "w");
        }

        //powtarzamy 3 razy dla roznych n
        for (int p = 1; p <= 3; p++)
        {
            int n;

            if (p == 1)
            {
                n = 5;
            }
            else if (p == 2)
            {
                n = 8;
            }
            else if (p == 3)
            {
                n = 21;
            }

            //*to samo co opisane w funkcji secDeriv*//
            float *xw = vector(1, n);

            for (int i = 1; i <= n; i++)
            {
                float h = (x_max - x_min) / ((float)(n)-1);
                xw[i] = x_min + h * (i - 1);
            }

            float *yw = vector(1, n);

            for (int i = 1; i <= n; i++)
            {
                if (funkcja == 1)
                {
                    yw[i] = f1(xw[i]);
                }
                else
                {
                    yw[i] = f2(xw[i]);
                }
            }

            float *m = vector(1, n);

            wyzM(xw, yw, m, n, alfa, beta);

            float arg;

            for (float i = -5; i < 5; i += 0.01)
            {
                arg = i;
                fprintf(file, "%f  %f\n", arg, wyzSx(xw, yw, m, n, arg));
            }

            fprintf(file, "\n\n");

            free_vector(xw, 1, n);
            free_vector(yw, 1, n);
            free_vector(m, 1, n);
        }

        fclose(file);

        if (funkcja == 1)
            secDeriv();
    }
    return 0;
}
