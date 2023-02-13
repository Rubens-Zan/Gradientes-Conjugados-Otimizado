#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.


// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  double *A; // coeficientes nao nulos da matriz, seguindo a ordem [diagonal_superior, diagonais_inferiores, diagonais_superiores]
  double *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;


SistLinear_t* alocaSisLin (unsigned int n, unsigned int k);
unsigned int indexa (int i, int j, int nDiagonais, unsigned int tamSL); 
void iniSisLin (SistLinear_t *SL, unsigned int nDiagonais);
void liberaSisLin (SistLinear_t *SL);

void prnVetor (double *vet, unsigned int n);
unsigned int calcElementoDiagonaisInf(unsigned int tamSl, unsigned int nDiagonais); 
double multiplicaVetores(double *vetA, double *vetB, unsigned int n);

// void copiaVetor (double *a, double *b, unsigned int N);

// void calcularAtxB(SistLinear_t *SL, double *b); 
// void calcularMatrizAtxA(SistLinear_t *SL, double **A);

#endif // __SISLIN_H__

