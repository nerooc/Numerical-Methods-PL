#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include </usr/include/gsl/gsl_math.h>
#include </usr/include/gsl/gsl_linalg.h>
#include </usr/include/gsl/gsl_errno.h>
#include </usr/include/gsl/gsl_fft_complex.h>

#define PI 3.14159265359

int main()
{
	srand(time(NULL)); //ziarno

	//deklarujemy parametry

	//okres
	double T = 1.0;

	//odchylenie
	double sigma = T / 20;

	//maksymalny okres czasu rejestracji sygnalu
	double tmax = 3 * T;

	//krok
	double d_t;

	//pulsacja
	double omega = 2 * PI / T;

	for (int k = 8; k <= 12; k += 2)
	{
		const int N = pow(2, k);
		d_t = tmax / N;
		FILE *file;

		double f[2 * N];
		double f0[2 * N];
		double g1[2 * N];
		double g2[2 * N];

		if (k == 8)
		{
			file = fopen("k8.dat", "w");
		}
		else if (k == 10)
		{
			file = fopen("k10.dat", "w");
		}
		else
		{
			file = fopen("k12.dat", "w");
		}

		for (int l = 0; l < N; l++)
		{
			double t = d_t * l;

			f[2 * l] = sin(omega * t) + sin(2 * omega * t) + sin(3 * omega * t) + rand() / (RAND_MAX + 1.0) - 0.5;

			f0[2 * l] = sin(omega * t) + sin(2 * omega * t) + sin(3 * omega * t);

			f[2 * l + 1] = 0;
			f0[2 * l + 1] = 0;

			g1[2 * l] = 1.0 / (sigma * sqrt(2 * PI)) * exp(-(t * t) / (2 * sigma * sigma));
			g1[2 * l + 1] = 0;

			g2[2 * l] = g1[2 * l];
			g2[2 * l + 1] = g1[2 * l + 1];

			fprintf(file, "%f %f\n", t, f[2 * l]);
		}

		gsl_fft_complex_radix2_forward(g1, 1, N);
		gsl_fft_complex_radix2_backward(g2, 1, N);
		gsl_fft_complex_radix2_forward(f, 1, N);

		for (int m = 0; m < N; m++)
		{
			double a1 = f[2 * m];
			double b1 = f[2 * m + 1];
			double a2 = g1[2 * m] + g2[2 * m];
			double b2 = g1[2 * m + 1] + g2[2 * m + 1];

			f[2 * m] = a1 * a2 - b1 * b2;
			f[2 * m + 1] = a1 * b2 + a2 * b1;
		}

		gsl_fft_complex_radix2_backward(f, 1, N);

		double f_max = pow(f[0], 2) + pow(f[1], 2);

		for (int n = 1; n < N; n++)
		{
			f_max < pow(f[2 * n], 2) + pow(f[2 * n + 1], 2) ? f_max = pow(f[2 * n], 2) + pow(f[2 * n + 1], 2) : 0;
		}

		f_max = sqrt(f_max);

		//robimy przerwy miedzy danymi
		fprintf(file, "\n\n");

		for (int o = 0; o < N; o++)
		{
			double t = o * d_t;
			fprintf(file, "%f %f\n", t, f[2 * o] * 2.5 / f_max);
		}

		fclose(file);
	}

	return 0;
}
