#include "gsl/gsl_linalg.h"
#define N 4


void main(){

    int signum;
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    
	for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
           gsl_matrix_set(A, i, j, 1./(i+j+2)); //dodajemy 2 poniewaz korzystamy z biblioteki gsl
        }
    }
	
	printf("Macierz A:\n");
	
    for (int i = 0 ; i < N ; i++){
        for (int j = 0; j < N; j++){
            printf("%13g ", gsl_matrix_get(A, i, j)); 
        }
        
        printf("\n");
    }

    printf("\n");    
    
    //ZADANIE 1
    
    gsl_permutation *p = gsl_permutation_calloc(N);

    gsl_linalg_LU_decomp(A, p, &signum);

    printf("//ZADANIE 1//Macierz LU: \n\n");
	
        for (int i = 0 ; i < N ; i++){
            for (int j = 0; j < N; j++){
                printf("%13g ", gsl_matrix_get(A, i, j)); //najkrÃ³tsza moÅ¼liwa reprezentacja wartoÅ›ci float/double
            }
            printf("\n");
        }

    printf("\n");   

    //ZADANIE 2
    float detA = 1;
	
    printf("//ZADANIE 2//Elementy diagonali macierzy LU:\n\n");
    
	  for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i == j){
				printf("%13g ", gsl_matrix_get(A, i, j));
				detA *= gsl_matrix_get(A, i, j);
			}
		}
	}
	printf("\n");
   
    printf("\n\ndet(A)=%g\n\n", (detA*signum));
    

    //ZADANIE 3
    gsl_matrix *Matrix = gsl_matrix_calloc(N, N);
    
    gsl_vector *b = gsl_vector_calloc(N);
    gsl_vector *x = gsl_vector_calloc(N);
    
    printf("//ZADANIE 3// Macierz odwrotna A^-1: \n\n");

    gsl_matrix *Inversed= gsl_matrix_calloc(N, N); 
    
	for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
           gsl_matrix_set(Inversed, i, j, 0);
        }
    }
    
    for (int i = 0; i < N; i++)
    {    

        for(int j = 0; j < N; j++){
            gsl_vector_set(b, j, 0);
        }

        gsl_vector_set(b, i, 1);
        
    

        gsl_linalg_LU_solve(A, p, b, x);

        for (int k = 0; k < N; k++){

            gsl_matrix_set(Inversed, k, i, gsl_vector_get(x, k));
            printf("%13g ", gsl_vector_get(x, k));
        }
        
        printf("\n");
    }   

    printf("\n");
    
/*
    printf("Macierz odwrotna 2:\n");
	
    for (int i = 0 ; i < N ; i++){
        for (int j = 0; j < N; j++){
            printf("%13g ", gsl_matrix_get(Inversed, i, j)); 
        }
        
        printf("\n");
    }

    printf("\n");   

*/
    //ZADANIE 4
    gsl_matrix *Product= gsl_matrix_calloc(N, N); 
    
	for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
           gsl_matrix_set(Product, i, j, 0);
        }
    }

    gsl_matrix *Ax= gsl_matrix_calloc(N, N); //moglbym skopiowac wczesniej, ale zrobilem tak :)
    
	for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
           gsl_matrix_set(Ax, i, j, 1./(i+j+2));
        }
    }
    
       
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Ax, Inversed, 0, Product); 


    printf("//ZADANIE 4// Wynik mnozenia macierzy:\n\n");
	
    for (int i = 0 ; i < N ; i++){
        for (int j = 0; j < N; j++){
            printf("%13g ", gsl_matrix_get(Product, i, j)); 
        }
        
        printf("\n");
    }

    printf("\n");


    //ZADANIE 5
    

    
    double max = gsl_matrix_get(Ax, 0, 0);
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N; j++){
            if (gsl_matrix_get(Ax, i, j) > max){
                max = gsl_matrix_get(Ax, i, j);
            }
        }
    }
    
    printf("//ZADANIE 5// \n\nMaksymalny element macierzy A: %g\n\n", max);
    
    double max2 = gsl_matrix_get(Inversed, 0, 0);
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N; j++){
            if (gsl_matrix_get(Inversed, i, j) > max2){
                max2 = gsl_matrix_get(Inversed, i, j);
            }
        }
    }
    
    printf("Maksymalny element macierzy odwrotnej: %g\n\n", max2);
    
    printf("Wskaznik uwarunkowania: %g\n\n", max*max2);

    

    //ZWALNIANIE PAMIECI
    
    gsl_matrix_free(A);
    gsl_matrix_free(Inversed);
    gsl_matrix_free(Ax);
    gsl_matrix_free(Product);
    gsl_vector_free(b);
    gsl_vector_free(x);
}

//INFORMACJA ZWROTNA
/*
Przy obliczaniu normy (czyli u nas: maksimum) brakuje wartoœci bezwzglêdnej, która jest wymagana wed³ug wzoru (5) oraz wed³ug definicji normy w ogólnoœci (norma jest uogólnieniem d³ugoœci, a d³ugoœæ nie mo¿e byæ ujemna :) ). U nas nie robi to ró¿nicy w wynikach, poniewa¿ najwiêkszy element z dok³adnoœci¹ do wartoœci bezwzglêdnej i tak jest oryginalnie dodatni, ale nie dla ka¿dej macierzy by tak by³o.
Wszystkie zaalokowane struktury wymagaj¹ zwolnienia -- przy wektorze permutacji oraz nieu¿ywanej macierzy Matrix tego brakuje. 
*/

