#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.
#define ALIGNMENT 16


// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  double *A; // coeficientes nao nulos da matriz, seguindo a ordem [diagonal_superior, diagonais_inferiores, diagonais_superiores]
  double *b; // termos independentes
  unsigned int n; // tamanho do SL
  unsigned int tam; 
  unsigned int nDiagonais; // numero de diagonais do SL
} SistLinear_t;

unsigned int indexMap(unsigned int i, unsigned int j, unsigned int k);

void calcMatrizTransp(SistLinear_t *SL, SistLinear_t *SLtransp); 
SistLinear_t* alocaSisLin (unsigned int n, unsigned int k);
void liberaSisLin (SistLinear_t *SL);
void iniSisLin(SistLinear_t *SL, unsigned int nDiagonais, double *matT);

void prnVetorArq(double *v, unsigned int n, FILE *arqSaida);
void prnVetor (double *vet, unsigned int n);
double multiplicaVetores( double  * restrict vetA,  double  * restrict vetB, unsigned int n);

void copiaVetor ( double  * restrict a,  double  * restrict b, unsigned int N);
void calcularAtxB(unsigned int k, unsigned int n, double *matrizTransp, double *vb, double *atvb);
void calcularMatrizAtxA( SistLinear_t   *restrict SL ,double *  matT ,SistLinear_t   *restrict novoSisLin );


#endif // __SISLIN_H__

