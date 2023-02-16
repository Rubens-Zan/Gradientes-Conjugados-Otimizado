#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
// #include "utils.h"
#include <string.h>

#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

// Alocaçao de matriz em memória.
SistLinear_t *alocaSisLin(unsigned int n, unsigned int k)
{
  SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));

  if (SL)
  {

    SL->n = n;
    SL->A = (double *)malloc((n +1)+ (2 * calcElementoDiagonaisInf(n,k) + 100)  * sizeof(double));
    SL->b = (double *)malloc(n * sizeof(double));

    if (!(SL->A) || !(SL->b))
    {
      liberaSisLin(SL);
      return NULL;
    }
  }

  return (SL);
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

    if (SL->b)
      free(SL->b);

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
 * @brief - Gera uma diagonal da matriz 
 * 
 * @param A - Vetor de diagonais principais
 * @param nElementos - n elementos na diagonal
 * @param totalElementos - n elementos na diagonal principal
 * @param nDiagonais - n de diagonais
 * @param indexInicial - index inicial do inicio da diagonal no vetor
 */
void geraDiagonal(double *A,unsigned int nElementos,unsigned int totalElementos,unsigned int nDiagonais, unsigned int indexInicial){
  unsigned int idx = indexInicial;
  #ifdef DEBUG
  printf("\n");
  #endif 
  if (nElementos == totalElementos)
    for (unsigned int i = 0;i < nElementos;++i){
      A[idx]= generateRandomA(i,i, nDiagonais); 
      #ifdef DEBUG
      printf("%.5g ",A[idx]);
      #endif
      ++idx; 
    }
  else
    for (unsigned int i = 0;i < nElementos;++i){
      A[idx]= generateRandomA(i,i+1, nDiagonais); 
      #ifdef DEBUG
      printf("%.5g ",A[idx]);
      #endif
      ++idx; 
    } 
}

unsigned int calcElementoDiagonaisInf(unsigned int tamSl, unsigned int nDiagonais){
  unsigned int elementoesDiagsInfs=0; 
  unsigned int elemDiag = tamSl;
  for (unsigned int i=0; i < (nDiagonais /2);++i){
    --elemDiag;
    elementoesDiagsInfs+= elemDiag;
  }

  return elementoesDiagsInfs;
}

/**
 * @brief - Inicia o sistema linear e aplica as propriedades para garantir a convergência
 *
 * @param SL - Sistema linear a ser iniciado
 * @param nDiagonais - Numero de diagonais
 */
void iniSisLin(SistLinear_t *SL, unsigned int nDiagonais)
{
  unsigned int indMatriz = 1,indMatrizInf=SL->n+1,indMatrizSup=indMatrizInf+calcElementoDiagonaisInf(SL->n, nDiagonais);
  unsigned int elementosDiagonal = SL->n - 1;

  SL->nDiagonais = nDiagonais; 
  SL->A[0]=0.0; 
  // o residuo do unroll eh a diagonal principal, como tem sempre k-impar diagonais na matriz
  geraDiagonal(SL->A,SL->n,SL->n, nDiagonais,indMatriz);
  indMatriz+= SL->n + 1; 
  // unroll de 2 
  for (int i = 1; i < nDiagonais; i+=2) {
    geraDiagonal(SL->A, elementosDiagonal, SL->n, nDiagonais, indMatrizInf);
    geraDiagonal(SL->A, elementosDiagonal, SL->n, nDiagonais, indMatrizSup);
    indMatrizInf+= elementosDiagonal;
    indMatrizSup+=elementosDiagonal; 

    --elementosDiagonal; 
  }

  #ifdef DEBUG
  for (int lin=0;lin < SL->n; ++lin){
    printf("\n"); 
    for (int col=0;col < SL->n; ++col){
      printf("%.5g |",SL->A[indexa(lin,col,nDiagonais, SL->n)]);
    }
  }
  #endif
  
  for (unsigned int i=0;i < SL->n;++i)
    SL->b[i] = generateRandomB(nDiagonais);

}


//Função que dado o indice da matriz, retorna o indice do elemento no vetor

/**
 * @brief Funcao para indexar a matriz, dado que sao armazenadas apenas as k-diagonais
 * [0, DIAG_PRINC, DIAG_INFS, DIAG_SUPS]
 * @param i - Linha  
 * @param j - Coluna
 * @param nDiagonais - N totais de diagonais na matriz 
 * @param tamSL - Dimensao do SL
 * @return unsigned int 
 */
unsigned int indexa (int i, int j, int nDiagonais, unsigned int tamSL) {
  // como a primeira diagonal eh a diagonal principal
  unsigned int distDiagPrinc = ABS(i - j);
  unsigned int diagonalIdx = (i > j) ? 1 : 0;
  unsigned int nElemDiag = tamSL - distDiagPrinc; 
  unsigned int meuIdx=0;
  
  for (unsigned int dist = distDiagPrinc;dist > 1 ;--dist){
    meuIdx+= (tamSL - dist);
  }

  // posicao igual a 0 
  if (distDiagPrinc > (nDiagonais / 2)){
    return 0; 
  }
  if (j == i)
    return (i+1); //posicao 0 guardada pro 0 
  // uma das diagonais inferiores
  else if (i > j){
    unsigned int comecoBloco = tamSL;
    return(comecoBloco + meuIdx + i ); 
  // uma das diagonais superiores
  }else {
    unsigned int comecoBloco = calcElementoDiagonaisInf(tamSL, nDiagonais);
    comecoBloco += tamSL; 
    return ( comecoBloco + meuIdx + j);
  }
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
double multiplicaVetores(double *vetA, double *vetB, unsigned int n)
{
  double produto = 0;

  for (int i = 0; i < n; ++i)
  {
    produto = produto + vetA[i] * vetB[i];
    // Testa valores inválidos.
    if (isnan(produto) || isinf(produto))
    {
      fprintf(stderr, "Erro variavel invalida: produto(multiplicaVetores): %g é NaN ou +/-Infinito\n", produto);
      exit(1);
    }
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
void copiaVetor(double *a, double *b, unsigned int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    b[i] = a[i];
  }
}
