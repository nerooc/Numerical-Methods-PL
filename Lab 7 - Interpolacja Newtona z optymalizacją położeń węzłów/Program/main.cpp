#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#define PI 3.14159265359

//FUNKCJA WYKONUJACA INTERPOLACJE
double interpol(double *p, double *xm, int n, double x)
{
    double w;
    double suma = 0;
    int j, i;
    for (j = 0; j <= n; j++)
    {
        w = 1;
        w *= p[j];

        for (i = 0; i <= j - 1; i++)
        {
            w *= (x - xm[i]);
        }

        suma += w;
    }

    return suma;
}

int main()
{

    for (int wezly = 1; wezly <= 2; wezly++)
    {
        std::ofstream data;
        if (wezly == 1)
        {
            data.open("zad_1.dat");
        }
        else if (wezly == 2)
        {
            data.open("zad_2.dat");
        }

        for (int n = 5; n <= 20; n += 5)
        {

            double h;

            //TABLICA
            double *xm = new double[n + 1];
            double *ym = new double[n + 1];

            //MIEJSCE NA TABLICE DWUWYMIAROWA
            double **fm = new double *[n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                fm[i] = new double[n + 1];
            }

            double x_min = -5.0; //X MINIMALNE
            double x_max = 5.0;  //X MAKSYMALNE

            std::cout << "\n\n TABLICA XM dla n = " << n << "\n";

            for (int l = 0; l < n + 1; l++)
            {
                if (wezly == 1)
                {
                    h = (x_max - x_min) / n;
                    xm[l] = x_min + h * l;
                }
                else if (wezly == 2)
                {
                    xm[l] = 0.5 * ((x_min - x_max) * cos(PI * (2 * l + 1) / (2 * n + 2)) + (x_min + x_max));
                }

                ym[l] = 1 / (1 + xm[l] * xm[l]);
                std::cout << "xm[" << l << "] = " << xm[l] << std::endl;
                fm[l][0] = ym[l];
            }

            std::cout << "\n\n TABLICA YM dla n = " << n << "\n";

            for (int l = 0; l < n + 1; l++)
            {
                std::cout << "ym[" << l << "] = " << ym[l] << std::endl;
            }

            //inicjalizujemy zerami
            for (int j = 1; j <= n; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    fm[i][j] = 0;
                }
            }

            for (int j = 1; j <= n; j++)
            {
                for (int i = j; i <= n; i++)
                {
                    fm[i][j] = (fm[i][j - 1] - fm[i - 1][j - 1]) / (xm[i] - xm[i - j]);
                }
            }

            //TABLICA ZAWIERAJACA JEDYNIE ELEMENTY PRZEKATNEJ
            double *tab = new double[n + 1];
            for (int j = 0; j <= n; j++)
            {
                tab[j] = fm[j][j];
            }

            double arg = 0;

            //WYPISUJEMY ELEMENTY WRAZ Z WYNIKAMI INTERPOLACJI
            for (double i = -5; i < 5; i += 0.01)
            {
                arg = i;
                data << arg << " " << interpol(tab, xm, n, arg) << std::endl;
            }

            data << "\n\n";

            //ZWALNIAMY PAMIEC
            for (int i = 0; i < n + 1; i++)
            {
                delete[](*fm);
            }

            delete fm;
        }

        data.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
