#include "sislin.h"
#include "utils.h"
#include "resolvedorGradConjug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define UNROLL 8
#define STRIDE 8
#define ALIGNMENT 16
#define ZERO 1e-30


/**
 * @brief Função para calcular o residuo original conforme A e b originais
 *
 * @param A - Coeficientes A originais
 * @param b - Coeficientes B originais
 * @param x - Solução final
 * @param n - Tamanho do Sistema Linear
 * @param residuo - Vetor a receber o residuo
 */
static inline void calculaResiduoOriginal(SistLinear_t *SL, double *x, double *residuo)
{
    memcpy(residuo, SL->b, (SL->n) * sizeof(double)); // residuo = b

    // Unroll and Jam residuo[]
    for (unsigned int i = 0; i < ((SL->n) - ((SL->n) % UNROLL)); i += UNROLL) // r = b - Ax
        for (unsigned int j = 0; j < (SL->n); ++j)
            for (unsigned int k = 0; k < UNROLL; ++k)
                residuo[i + k] -= SL->A[indexMap(i + k + 1, j + 1, SL->nDiagonais)] * x[j+1];

    // remaider loop
    for (unsigned int i = ((SL->n) - ((SL->n) % UNROLL)); i < (SL->n); ++i)
        for (unsigned int j = 0; j < (SL->n); ++j)
            residuo[i] -= SL->A[indexMap(i + 1, j + 1, SL->nDiagonais)] * x[j+1];
}

/******************FUNCOES GRADIENTE CONJUGADO com pre condicionador **********************************/
/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado com uso do pre condicionador de jacobi
 *
 * @param SL - Sistema Linear sendo ja aplicados os metodos para conversão
 * @param maxIt - numero maximo de iterações
 * @param tol - Tolerancia maxima
 * @param matSaida -
 * @param arqSaida - Arquivo para printar as informações
 */
void gradienteConjugadoPreCondic(SistLinear_t *SL, int maxIt, double tol, FILE *arqSaida)
{
	double alpha;			// s
	double beta;			// m
	double normx;			// ||x||
	double relerr;			// erro relativo atual
	double aux0, aux1;		// aux, aux1
	double *xAnt;			// vetor x anterior
	double *v;				// v
	double *z;			    // z
	double *y;			    // y
	unsigned int numiter;	// numero de iteracoes
	unsigned int indxmax;	// indice no qual max(|xatual - xold|)
	double soma;			// variavel auxiliar nos lacos

    double *simmat = SL->A;
    unsigned int n =SL->n;
    double *atvb = SL->b;
    double *res = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));;
    double *vetx =aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	double k = SL->nDiagonais;

    indxmax = 0;

	z = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	xAnt = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	memcpy(xAnt, vetx, (n + 1)*sizeof(double));	// x0 = 0
	memcpy(res, atvb, (n + 1)*sizeof(double));		// r = b
	y = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	v = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	
	
	//! Variáveis:
	double alpha;			//! s
	double beta;			//! m
	double normx;			//! ||x||
	double relerr;			//! erro relativo atual
	double aux0, aux1;		//! aux, aux1
	double *vetxold;		//! vetor x anterior
	double *vetv;			//! v
	double *vetz;			//! z
	double *vety;			//! y
	unsigned int numiter;	//! numero de iteracoes
	unsigned int indxmax;	//! indice no qual max(|xatual - xold|)
	double soma;			//! variavel auxiliar nos lacos
	//!

	indxmax = 0;

	//! Algoritmo:
	vetv = malloc((n + 1) * sizeof(double));
	vetz = malloc((n + 1) * sizeof(double));
	vety = malloc((n + 1) * sizeof(double));
	vetxold = malloc((n + 1) * sizeof(double));
	memcpy(vetxold, vetx, (n + 1)*sizeof(double));	//! x0 = 0
	memcpy(res, atvb, (n + 1)*sizeof(double));		//! r = b
	
	for(unsigned int g = 1; g < (n + 1); g++)		//! v = M-1b, y = M-1r 
	{
		vetv[g] = (atvb[g] / obtemValor(g, g, k, n, simmat));
		vety[g] = (res[g] / obtemValor(g, g, k, n, simmat)); 
	}
	aux0 = 0.0;
	for(unsigned int g = 1; g < (n + 1); g++)		//! aux = ytr
	{
		aux0 += vety[g] * res[g]; 
	}
	
	LIKWID_MARKER_START("op1");
 
	numiter = 0;	
	do 									//! para k = 0 : max, faca
	{
		numiter++;
		
		for(unsigned int a = 1; a < (n + 1); a++)	//! z = Av
		{
			soma = 0.0;
			for(unsigned int b = 1; b < (n + 1); b++)
			{
				soma += (obtemValor(a, b, k, n, simmat)) * vetv[b];
			}
			vetz[a] = soma;
		}

		soma = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++)	//! calcula vtz
		{
			soma += vetv[g] * vetz[g]; 
		}
		//! s = aux/vtz
		alpha = (fabs(soma) < ZERO) ? 0.0 : (aux0 / soma);

		for(unsigned int g = 1; g < (n + 1); g++)	//! xk+1 = xk + sv
		{
			vetx[g] = vetxold[g] + (alpha * vetv[g]); 
		}

		for(unsigned int g = 1; g < (n + 1); g++)	//! r = r - sz
		{
			res[g] = res[g] - (alpha * vetz[g]); 
		}

		for(unsigned int g = 1; g < (n + 1); g++)	//! y = M-1r
		{
			vety[g] = (res[g] / obtemValor(g, g, k, n, simmat)); 
		}

		normx = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++)	//! calcula ||x||
		{
			if (normx < fabs(vetx[g] - vetxold[g]))
			{
				normx = fabs(vetx[g] - vetxold[g]);
				indxmax = g; 
			}
		}		

		aux1 = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++)	//! aux1 = ytr
		{
			aux1 += vety[g] * res[g]; 
		}

		//! m = aux1 / aux
		beta = (fabs(aux0) < ZERO) ? 0.0 : (aux1 / aux0);
		aux0 = aux1;

		for(unsigned int g = 1; g < (n + 1); g++)	//! v = y + mv
		{
			vetv[g] = vety[g] + (beta * vetv[g]); 
		}

		#ifdef DEBUG
			printf("\n");
			printf("||X%02u|| = %-11.7g \n", numiter, normx);
			printf("\n");
		#endif		

		fprintf(arq, "# iter %u: %.15g\n", numiter, normx);

		//! xold = x
		memcpy(vetxold, vetx, (n + 1)*sizeof(double));
		//! relerr = max(|Xi - Xi-1|) / Xi
		relerr = (fabs(vetx[indxmax]) < ZERO) ? 0.0 : (normx / fabs(vetx[indxmax]));
		//!
	} while ((numiter < iter) && (relerr > err));
	
	LIKWID_MARKER_STOP("op1");
 
	free(vetv);
	free(vetz);
	free(vety);
	free(vetxold);
}
