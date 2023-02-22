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
    double tMedioIter, tempoResid, tempoPreCond, tempoInicio;

    double *simmat = SL->A;
    unsigned int n =SL->n;
    double *atvb = SL->b;
    double *res = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
    double *vetx =aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	double k = SL->nDiagonais;
    tempoPreCond = timestamp();
    indxmax = 0;

	z = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	xAnt = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	memcpy(xAnt, vetx, (n + 1)*sizeof(double));	// x0 = 0
	memcpy(res, atvb, (n + 1)*sizeof(double));		// r = b
	y = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	v = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	
	aux0 = 0.0;
	for(unsigned int i = 1; i < (n + 1); i++)		// v = M-1b, y = M-1r 
	{
		v[i] = (atvb[i] / simmat[indexMap(i,i,k)]);
		y[i] = (res[i] / simmat[indexMap(i,i,k)]);
		aux0 += y[i] * res[i];	// aux = ytr
	}
	LIKWID_MARKER_START("op1");
	
	numiter = 0;	
	while(numiter < maxIt) 									// para k = 0 : max, faca
	{	
		
		soma =  0.0;
		//Unroll and Jam z = Av e soma =  vtz
		memset(z, 0, (n + 1));
		for(unsigned int i =  1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	// z = Av
		{
			for(unsigned int j = 1; j < (n + 1); j++)
			{
				for(unsigned int c = 0; c < STRIDE; ++c)
				{
					z[i + c] += simmat[indexMap((i + c), j, k)] * v[j];
				}
			}
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				soma += z[i + c] * v[i + c];	// calcula vtz
			}
		}

		//remaider loop
		for(unsigned int i =  (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)
		{
			for(unsigned int j = 1; j < (n + 1); ++j)
			{
				z[i] += simmat[indexMap(i,j,k)] * v[j];
			}
			soma += z[i] * v[i];
		}

		// s = aux/vtz
		alpha =  (fabs(soma) < ZERO) ? 0.0 : (aux0 / soma);

		normx = 0.0;
		// CALCULO DA NORMA
		//Unroll and Jam do calculo da norma e x = xold
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				vetx[i + j] += (alpha * v[i + j]);			// xk+1 = xk + sv
				if (normx < fabs(vetx[i + j] - xAnt[i + j]))	// calcula ||x||
				{
					normx = fabs(vetx[i + j] - xAnt[i + j]);
					indxmax = i+j; 
				}
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//remainder loop
		{
			vetx[i] += (alpha * v[i]);			// xk+1 = xk + sv
			if (normx < fabs(vetx[i] - xAnt[i]))	// calcula ||x||
			{
				normx = fabs(vetx[i] - xAnt[i]);
				indxmax = i; 
			}  
		}
		
		//Unroll and Jam de r = r - sz e y = M-1r
		// CALCULO DO RESID E Y
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				res[i + j] -= (alpha * z[i + j]);
				y[i + j] = (res[i + j] / simmat[indexMap(i+j,i+j,k)]);
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//remainder loop
		{
			res[i] -= (alpha * z[i]);
			y[i] = (res[i] / simmat[indexMap(i,i,k)]);  
		}
		aux1 = 0.0;

		//Unroll and Jam de aux1 = ytr
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				aux1 += y[i + j] * res[i + j];
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//remainder loop
		{
			aux1 += y[i] * res[i];  
		}
		
		// m = aux1 / aux
		beta =  (aux1 / aux0);
		aux0 = aux1;

		//Unroll and Jam de v = y + mv
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				v[i + j] = y[i + j] + (beta * v[i + j]);
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//remainder loop
		{
			v[i] = y[i] + (beta * v[i]);  
		}	
		fprintf(arqSaida, "# iter %u: %.15g\n", numiter, normx);

		// xold = x
		memcpy(xAnt, vetx, (n + 1)*sizeof(double));
		// relerr = max(|Xi - Xi-1|) / Xi
		relerr =(normx / fabs(vetx[indxmax]));
		//
		++numiter;
	};
	
	LIKWID_MARKER_STOP("op1");
 	
	tempoResid = timestamp();
    fprintf(arqSaida, "# residuo: || %.15g || \n", normaL2Residuo(res, SL->n));
    tempoResid = timestamp() - tempoResid; 

    // tempo final
    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tempoPreCond /  numiter;
    fprintf(arqSaida, "# Tempo iter: %.15g \n", tMedioIter);
	LIKWID_MARKER_START("op2");
    calculaResiduoOriginal(SL,vetx, res);
	LIKWID_MARKER_STOP("op2");
    fprintf(arqSaida, "# Tempo residuo: %.15g \n", tempoResid);
    fprintf(arqSaida, "# \n");
    fprintf(arqSaida, "%d", SL->n);
    prnVetorArq(vetx, SL->n, arqSaida);
	
	
	free(v);
	free(res);
	free(z);
	free(y);
	free(vetx);
	free(xAnt);
}
