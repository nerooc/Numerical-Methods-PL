#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

void licz(double bet, double F0, double Omega, std::string x)
{

	//STATYCZNA ZMIENNA POZWALAJĄCA NA ZMIANĘ NAZW PLIKÓW POMIĘDZY PRZYPADKAMI
	static int nrIteracji = 0;

	//LICZBA KROKÓW CZASOWYCH
	int n = 2000;

	//ELEMENTY PRZETĄTNYCH MACIERZY TRÓJPRZEKĄTNIOWEJ
	double a1;
	double a2;
	double a3;

	//KROK CZASOWY
	double h = 0.02;

	//CZĘSTOŚĆ KĄTOWA
	double om = 1; // om - mala omega, Omega - duza omega

	//TRZY N-ELEMENTOWE WEKTORY PRZECHOWUJĄCE TRZY PRZEKĄTNE MACIERZY RZADKIEJ
	std::vector<double> d0;
	std::vector<double> d1;
	std::vector<double> d2;

	//WEKTORY PRZYBLIŻENIA
	std::vector<double> Xn;
	std::vector<double> Xs;

	//WEKTOR ROZWIĄZAŃ
	std::vector<double> b;

	//NADAJEMY ELEMENTOM WARTOŚCI
	a1 = 1;

	a2 = (om * om) * (h * h) - 2 - bet * h;

	a3 = 1 + bet * h;

	//TESTOWE WYPISANIE ELEMENTÓW:
	std::cout << "\n"
			  << x << std::endl;
	std::cout << "a1 = " << a1 << " "
			  << "a2 = " << a2 << " "
			  << "a3 = " << a3 << std::endl;

	//DWA PIERWSZE ELEMENTY TO 1 i 0
	for (int i = 1; i >= 0; i--)
	{
		b.push_back(i);
	}

	//RESZTA ELEMENTÓW OPISANA WZOREM
	for (int i = 1; i <= n; i++)
	{
		b.push_back(F0 * sin(Omega * h * i) * h * h);
	}

	//WSTĘPNIE WYPEŁNIAMY WEKTORY ZERAMI
	for (int i = 0; i <= n; i++)
	{
		Xn.push_back(0);
		Xs.push_back(0);
	}

	//TEST PIERWSZYCH ELEMENTÓW WEKTORA B
	for (int i = 0; i < 4; i++)
	{
		std::cout << "b[" << i << "] = " << b[i] << " " << std::endl;
	}

	//TEST OSTATNICH ELEMENTÓW WEKTORA B
	for (int i = n - 2; i <= n; i++)
	{
		std::cout << "b[" << i << "] = " << b[i] << " " << std::endl;
	}

	std::cout << "\n";

	//WYPEŁNIAMY WEKTORY PRZEKĄTNYCH
	d0.push_back(1); //d0 = [1, 1, a3, a3, ..., a3]
	d0.push_back(1);

	d1.push_back(0); //d1 = [0, -1, a2, a2, ..., a2]
	d1.push_back(-1);

	d2.push_back(0); //d2 = [0, 0, a1, a1, ..., a1]
	d2.push_back(0);

	for (int i = 3; i <= n; i++)
	{
		d0.push_back(a3);

		d1.push_back(a2);

		d2.push_back(a1);
	}

	//LICZNIK ITERACJI
	int counter = 0;

	//SUMA KWADRATÓW ELEMENTÓW WEKTORA Xs
	double Ss = 0;

	//SUMA KWADRATÓW ELEMENTÓW WEKTORA Xn
	double Sn = 0;

	//PĘTLA
	while (counter < 100000)
	{

		//WYPEŁNIENIE DWÓCH PIERWSZYCH ELEMENTÓW WEKTORA Xn
		Xn[0] = b[0] / d0[0]; //Xs[-2] i Xs[-1] MOGĄ BYĆ DOWOLNE, DLATEGO TUTAJ = 0
		Xn[1] = (1 / d0[1]) * (b[1] - d1[1] * Xs[0]);

		//WYPEŁNIENIE RESZTY NOWEGO WEKTORA Xn
		for (int i = 2; i <= n; i++)
		{
			Xn[i] = (1 / d0[i]) * (b[i] - d1[i] * Xs[i - 1] - d2[i] * Xs[i - 2]);
		}

		//WYZEROWANIE WARTOŚCI SUM KWADRATÓW Z POPRZEDNIEJ ITERACJI
		Ss = 0;
		Sn = 0;

		//OBLICZAMY SUMĘ KWADRATÓW ELEMENTÓW WEKTORÓW BY POTEM WYKORZYSTAĆ ICH RÓŻNICĘ W WARUNKU PRZERWANIA PĘTLI
		for (int i = 0; i < n; i++)
		{
			Ss += (Xs[i] * Xs[i]);
			Sn += (Xn[i] * Xn[i]);

			//PRZENOSIMY WEKTOR Xn Z TEJ ITERACJI DO WEKTORA Xs BY WYKORZYSTAĆ GO W NASTĘPNEJ
			Xs[i] = Xn[i];
		}

		//ZWIĘKSZAMY ILOŚĆ ITERACJI O JEDNĄ
		counter++;

		//BADANIE ZBIEŻNOŚCI - WARUNEK PRZERWANIA PĘTLI
		if (fabs(Sn - Ss) < 10e-6)
		{
			break;
		}
	}

	//COUT WYPISUJĄCY ILOŚĆ ITERACJI POLICZONĄ PRZEZ LICZNIK COUNTER
	std::cout << "Ile iteracji: " << counter << "\n"
			  << std::endl;

	//DEFINIUJEMY PLIK I JEGO NAZWĘ
	std::ofstream FILE;
	std::string fileName = "wyniki .txt";

	//NAZWA ZMIENIA NUMER W ZALEŻNOŚCI OD PRZYPADKU: wynik1.txt/wynik2.txt/wynik3.txt
	fileName[6] = '1' + nrIteracji;

	//OTWIERAMY PLIK
	FILE.open(fileName);

	//ZAPISUJEMY DO PLIKU W FORMACIE: 	ti   Xn[i]   ->   i*0.02   Xn[i]
	for (int i = 0; i < n; i++)
	{
		double ti = i * 0.02;
		FILE << ti << " " << Xn[i] << std::endl;
	}

	//ZAMYKAMY PLIK
	FILE.close();

	//INKREMENTUJEMY, ŻEBY NASTĘPNE WYWOŁANIE FUNKCJI ZAPISAŁO WYNIK W PLIKU Z CYFRĄ O 1 WYŻSZĄ
	nrIteracji++;
}

int main()
{
	licz(0.0, 0.0, 0.8, "/////PROBA PIERWSZA/////");
	licz(0.4, 0.0, 0.8, "/////PROBA DRUGA/////");
	licz(0.4, 0.1, 0.8, "/////PROBA TRZECIA/////");
}

/* Komentarz prowadzącej */
/*
Wyniki faktycznie są prawidłowe, więc nie musi Pan już niczego zmieniać, chociaż mam trochę uwag w sprawie zakończenia różnych pętli, o czym wspominał Pan w komentarzu. W każdym razie nie musi Pan już przesyłać programu ze sprawozdaniem -- ta wersja wystarczy.
Wszystkie wektory miały mieć n+1 elementów. Tymczasem w programie wektor b ma n+2 elementy (instrukcja b.push_back(F0 * sin(Omega * h * i) * h * h); powinna być wykonana n−1 razy, a nie n razy). Niczego to oczywiście nie psuje, ponieważ po prostu nie korzysta Pan z tego dodatkowego elementu na końcu.
Z kolei wszystkie trzy wektory d w Pańskim programie mają po n elementów, a powinny mieć o 1 więcej: pętla
	for (int i = 3; i <= n; i++)
	{
		d0.push_back(a3);
		d1.push_back(a2);
		d2.push_back(a1);
	}
powinna wykonywać o jeden obieg więcej. To jest już duży problem, ponieważ później w obliczeniach pętla metody używa tych nieistniejących elementów. Myślę, że dlatego miał Pan problem ze zbieżnością (i pewnie dlatego kończy Pan obliczenia sum Sn i Ss o jeden obieg za wcześnie). Jeśli dodamy w wektorach d0, d1, d2 ten brakujący element, to pętla metody (obliczanie Xn) działa prawidłowo, a pętlę obliczającą Ss oraz Sn można wydłużyć o 1 obieg, podobnie wypisywanie wyników.
Za to wektory Xn oraz Xs mają tyle elementów, ile powinny :)
Wektor Xn jest nowym oszacowaniem rozwiązania. Nowe oszacowanie bazuje tylko na starym Xs. Gdyby przyjrzeć się wzorowi metody Jacobiego, to można zauważyć, że w każdej iteracji metody tylko jeden z elementów ma szansę być obliczony prawidłowo -- najpierw pierwszy, potem drugi, itp., ponieważ wszystkie następne elementy nowego oszacowania są obliczane jeszcze na podstawie niedoszacowanych wartości z Xs (one są w danej iteracji stałe).
Jak Pan zauważył, metoda Gaussa-Seidla zbiega się dużo szybciej: nie, nie jest to błąd :) Powodem jest fakt, że w tej metodzie w każdej iteracji obliczenia dla danego elementu bazują na poprzednich już poprawionych elementach, skoro jest jedna tablica. Po wykonaniu obliczeń dla xi przechodzimy do kolejnego elementu, który już może bazować na poprawionym poprzednim (chociaż jesteśmy jeszcze w ten samej iteracji). W metoda Jacobiego było to możliwe dopiero w następnej iteracji.
/*