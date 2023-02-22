#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.


// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  double **A; // coeficientes
  double *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;


SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);
void iniSisLin (SistLinear_t *SL, unsigned int nDiagonais);
void prnMat (double **mat, unsigned int n, unsigned int m);
void prnVetorArq(double *v, unsigned int n, FILE *arqSaida);
void prnSisLin (SistLinear_t *SL);
void prnVetor (double *vet, unsigned int n);
double multiplicaVetores(double *vetA, double *vetB, unsigned int n);
void copiaVetor (double *a, double *b, unsigned int N);
void copiaMatriz(double **A, double **B, unsigned int n); 

void calcularAtxB(SistLinear_t *SL, double *b); 
void calcularMatrizAtxA(SistLinear_t *SL, double **A);

#endif // __SISLIN_H__

