#include "sislin.h"
#include "utils.h"
#include "resolvedorGradConjug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <likwid.h>
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
	double *vetv;				// v
	double *vetz;			    // z
	double *vety;			    // y
	unsigned int numIter;	// numero de iteracoes
	unsigned int indxmax;	// indice no qual max(|xatual - xold|)
	double soma;			// variavel auxiliar nos lacos
    double tMedioIter, tempoResid, tempoPreCond, tempoInicio;
    unsigned int n =SL->n;

    double *res = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));;
    double *vetx =aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));

    indxmax = 0;

	vetz = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	xAnt = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	memcpy(xAnt, vetx, (n + 1)*sizeof(double));	// x0 = 0
	memcpy(res, SL->b, (n + 1)*sizeof(double));		// r = b
	vety = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	vetv = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	

	//!

	indxmax = 0;

	//! Algoritmo:
	vetv = malloc((n + 1) * sizeof(double));
	vetz = malloc((n + 1) * sizeof(double));
	vety = malloc((n + 1) * sizeof(double));
	xAnt = malloc((n + 1) * sizeof(double));
	memcpy(xAnt, vetx, (n + 1)*sizeof(double));	//! x0 = 0
	memcpy(res, SL->b, (n + 1)*sizeof(double));		//! r = b
	
	for(unsigned int i = 1; i < (n + 1); i++)		//! v = M-1b, y = M-1r 
	{
		vetv[i] = (SL->b[i] / SL->A[indexMap(i,i,SL->nDiagonais)]);
		vety[i] = (res[i] / SL->A[indexMap(i,i,SL->nDiagonais)]); 
	}
	aux0 = 0.0;
	for(unsigned int i = 1; i < (n + 1); i++)		//! aux = ytr
	{
		aux0 += vety[i] * res[i]; 
	}
	
	LIKWID_MARKER_START("op1");
 
	numIter = 0;	
    tempoPreCond = timestamp();
	while (numIter < maxIt)
	{
		numIter++;
		
		for(unsigned int i = 1; i < (n + 1); i++)	//! z = Av
		{
			soma = 0.0;
			for(unsigned int j = 1; j < (n + 1); j++)
			{
				soma += (SL->A[indexMap(i,j,SL->nDiagonais)] ) * vetv[j];
			}
			vetz[i] = soma;
		}

		soma = 0.0;
		for(unsigned int i = 1; i < (n + 1); i++)	//! calcula vtz
		{
			soma += vetv[i] * vetz[i]; 
		}
		//! s = aux/vtz
		alpha = (aux0 / soma);

		for(unsigned int i = 1; i < (n + 1); i++)	//! xk+1 = xk + sv
		{
			vetx[i] = xAnt[i] + (alpha * vetv[i]); 
		}

		for(unsigned int i = 1; i < (n + 1); i++)	//! r = r - sz
		{
			res[i] = res[i] - (alpha * vetz[i]); 
		}

		for(unsigned int i = 1; i < (n + 1); i++)	//! y = M-1r
		{
			vety[i] = (res[i] / SL->A[indexMap(i,i,SL->nDiagonais)]); 
		}

		normx = 0.0;
		for(unsigned int i = 1; i < (n + 1); i++)	//! calcula ||x||
		{
			if (normx < fabs(vetx[i] - xAnt[i]))
			{
				normx = fabs(vetx[i] - xAnt[i]);
				indxmax = i; 
			}
		}		

		aux1 = 0.0;
		for(unsigned int i = 1; i < (n + 1); i++)	//! aux1 = ytr
		{
			aux1 += vety[i] * res[i]; 
		}

		//! m = aux1 / aux
		beta = (aux1 / aux0);
		aux0 = aux1;

		for(unsigned int i = 1; i < (n + 1); i++)	//! v = y + mv
		{
			vetv[i] = vety[i] + (beta * vetv[i]); 
		}

		fprintf(arqSaida, "# iter %u: ||%.15g||\n", numIter, normx);

		//! xold = x
		memcpy(xAnt, vetx, (n + 1)*sizeof(double));
		//! relerr = max(|Xi - Xi-1|) / Xi
		relerr = (normx / fabs(vetx[indxmax])); // nao utilizado nessa implementação
		//!
	};
	
	LIKWID_MARKER_STOP("op1");

	LIKWID_MARKER_START("op2");
    tempoResid = timestamp();

	calculaResiduoOriginal(SL, vetx, res); 
    tempoResid = timestamp() - tempoResid; 
    fprintf(arqSaida, "# residuo: || %.15g || \n", normaL2Residuo(res, SL->n));
  	
	// tempo final
  	tempoPreCond = timestamp() - tempoPreCond;

    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tempoPreCond / numIter;
    fprintf(arqSaida, "# Tempo iter: %.15g \n", tMedioIter);

    fprintf(arqSaida, "# Tempo residuo: %.15g \n", tempoResid);
    fprintf(arqSaida, "# \n");
    fprintf(arqSaida, "%d", SL->n);
    prnVetorArq(vetx, SL->n, arqSaida);

	LIKWID_MARKER_STOP("op2");
 
	free(vetv);
	free(vetz);
	free(vety);
	free(xAnt);
}
