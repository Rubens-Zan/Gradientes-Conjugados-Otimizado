#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
// #include "utils.h"
#include <string.h>

#define ABS(num)  ((num) < 0.0 ? -(num) : (num))
#define ALIGNMENT 16
// Alocaçao de matriz em memória.
SistLinear_t *alocaSisLin(unsigned int n, unsigned int k)
{
  SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));

  if (SL)
  {

    SL->n = n;
	  SL->A = aligned_alloc(ALIGNMENT, ((n * k * sizeof(double)) + (ALIGNMENT - ((n * k * sizeof(double)) % ALIGNMENT))));

    SL->b = (double *)malloc(n * sizeof(double)+1);
    SL->nDiagonais = k; 

    if (!(SL->A) || !(SL->b))
    {
      liberaSisLin(SL);
      return NULL;
    }
  }

  return (SL);
}
/**
 * @brief      Função que retorna o índice no array de diagonais de um elemento aij
 * @param      i   coordenada linha do elemento
 * @param      j   coordenada coluna do elemento
 * @param      k   quantidade de diagonais da matriz A
 * @return 	   indice de aij no array de diagonais
 */
unsigned int index(unsigned int i, unsigned int j, unsigned int k)
{
  if ((j > i ? (j - i) : (i - j)) > (k / 2))
  {
    return (0);
  }

	return ((k / 2) + ((i - 1) * k) + j - i);
}

// Liberacao de memória
void liberaSisLin(SistLinear_t *SL)
{
  if (SL)
  {
    if (SL->A)
    {
      free(SL->A);
    }

    // if (SL->b)
      // free(SL->b);
    
    free(SL);
  }
}

// /***********************
//  * Função que gera os coeficientes de um sistema linear k-diagonal
//  * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
//  * k: numero de diagonais da matriz A
//  ***********************/
static inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((i == j) ? (double)(k << 1) : 1.0) * (double)rand() * invRandMax;
}

// /***********************
//  * Função que gera os termos independentes de um sistema linear k-diagonal
//  * k: numero de diagonais da matriz A
//  ***********************/
static inline double generateRandomB(unsigned int k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k << 2) * (double)rand() * invRandMax;
}

/**
 * @brief - Inicia o sistema linear e aplica as propriedades para garantir a convergência
 *
 * @param SL - Sistema linear a ser iniciado
 * @param nDiagonais - Numero de diagonais
 */
void iniSisLin(SistLinear_t *SL, unsigned int nDiagonais, double *matT)
{
  double termo;
  SL->A[0] = 0.0;
  matT[0] = 0.0;
  SL->b[0] = 0.0;

	for (unsigned int i = 1; i <= SL->n; ++i)
	{
		SL->b[i] = generateRandomB(nDiagonais);

	 	for (unsigned int j = 1; j <= SL->n; ++j)
	 	{
	 		if ((j > i ? (j - i) : (i - j)) <= (nDiagonais / 2))
	 		{
	 			termo = generateRandomA(i, j, nDiagonais);
				
	 			SL->A[index(i, j, nDiagonais)] = termo;
	 			matT[index(j, i, nDiagonais)] = termo;	 			
	 		}
	 	}
	}
}

void calcMatrizTransp(SistLinear_t *SL, SistLinear_t *SLtransp){
  SLtransp->A[0] = 0.0;
  for (unsigned int i=1;i <= SL->n;++i){
    for (unsigned int j=1;j <= SL->n;++j){
      SLtransp->A[index(j,i, SL->nDiagonais)] = SL->A[index(i,j, SL->nDiagonais)]; 
    }
  }
}

void calcularAtxB(unsigned int k, unsigned int n, double *matrizTransp, double *vb, double *atvb)
{
  double at;
	double soma;

	for (unsigned int i = 1; i <= n ; ++i)
	{
		soma = 0.0;
	  for (unsigned int j = 1; j <= n; ++j)
		{
			at = matrizTransp[index(i,j,k)];
			soma += at * vb[j-1];
		}
		atvb[i] = soma;
	}

}

// o calculo de A*At vai gerar uma matriz com 2x mais diagonais nao nulas
void calcularMatrizAtxA( SistLinear_t   *restrict SL ,double *  matT ,SistLinear_t   *restrict novoSisLin )
{ 
  double a;
	double at;
	double soma;
	unsigned int simk = (2*SL->nDiagonais) - 1;

	for (unsigned int i = 1; i < (SL->n + 1); ++i)
	{
		for (unsigned int j = 1; j < (SL->n + 1); ++j)
		{
			if ((j > i ? (j - i) : (i - j)) <= (simk / 2))
	 		{
	 			soma = 0.0;
	 			for (unsigned int k = 1; k < (SL->n + 1); ++k)
	 			{
	 				at = matT[index(i,k,SL->nDiagonais)]; 
	 				a = SL->A[index(k,j, SL->nDiagonais)];
	 				soma += at * a;
	 			}

	 			novoSisLin->A[index(i, j, simk)] = soma;
	 		} 
		}
	}
}

/**
 * @brief - Printa o vetor V no arquivo arqSaida
 *
 * @param v - vetor
 * @param n - tamanho do vetor
 * @param arqSaida - arquivo de saida
 */
void prnVetorArq(double *v, unsigned int n, FILE *arqSaida)
{
  int i;
  fprintf(arqSaida, "\n");
  for (i = 0; i < n; ++i)
  {
    fprintf(arqSaida, "%.15g ", v[i]);
  }
  fprintf(arqSaida, "\n");
}


/***********************************************************************/

void prnVetor(double *v, unsigned int n)
{
  int i;

  printf("\n");
  for (i = 0; i < n; ++i)
    printf("%10g ", v[i]);
  printf("\n\n");
}

/**
 * @brief - Calcula o produto entre dois vetores que geram um valor
 * produto = vetA + vetB
 * @param vetA - Vetor A
 * @param vetB - Vetor B
 * @param n - Tamanho do vetor
 * @return double - Produto dos vetores
 */
double multiplicaVetores( double  *restrict vetA,  double  *restrict vetB, unsigned int n)
{
  double produto = 0;

  for (int i = 0; i < n; ++i)
  {
    produto = produto + vetA[i] * vetB[i];
  }
  return produto;
}

/**
 * @brief Copia vetor
 *
 * @param a vetor origem
 * @param b vetor destino
 * @param n tamanho do vetor
 */
void copiaVetor( double  *restrict a,  double  *restrict b, unsigned int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    b[i] = a[i];
  }
}
