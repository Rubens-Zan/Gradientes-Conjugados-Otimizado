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

    //! Unroll and Jam residuo[]
    for (unsigned int i = 0; i < ((SL->n) - ((SL->n) % UNROLL)); i += UNROLL) // r = b - Ax
        for (unsigned int j = 0; j < (SL->n); ++j)
            for (unsigned int k = 0; k < UNROLL; ++k)
                residuo[i + k] -= SL->A[indexMap(i + k + 1, j + 1, SL->nDiagonais)] * x[j+1];

    //! remaider loop
    for (unsigned int i = ((SL->n) - ((SL->n) % UNROLL)); i < (SL->n); ++i)
        for (unsigned int j = 0; j < (SL->n); ++j)
            residuo[i] -= SL->A[indexMap(i + 1, j + 1, SL->nDiagonais)] * x[j+1];
}

/***********************METODOS iguais nos dois tipos de resolução*****************************************/
/**
 * @brief Calcula o proximo x<k>
 *
 *  x[i][j]<k> = x[i][j]<k-1> + alpha <k> * p[i][j]<k-1>
 *
 * @param proxX
 * @param xAnt - Valor do x<i-1> anterior
 * @param alpha - valor do alpha<i-1> anterior
 * @param p - valor da direção de busca p<i-1)> anterior
 * @param n - dimensao do sistema linear
 * @return double - Proximo x calculado
 */
static inline void calcX(double *restrict proxX, double *restrict xAnt, double alpha, double *restrict p, int n)
{
    for (int i = 0; i < n; ++i)
    {
        proxX[i] = xAnt[i] + (alpha * p[i]);
    }
}

/**
 * @brief - Calcula o residuo
 *
 * r<k> = r<k-1> - alpha<k-1> * A * p<k-1>
 *
 * @param residuoAnterior - residuo anterior
 * @param alpha - Alpha calculado
 * @param A - matriz A
 * @param p - vetor de direção de busca
 * @param residuo - vetor de residuo
 * @param n - dimensao do sistema linear
 */
static inline void calcResiduo(double *restrict residuoAnterior, double alpha, SistLinear_t *SL, double *restrict p, double *restrict residuo)
{
    // alpha * A * p<k>
    for (int i = 0; i < SL->n; ++i)
    {
        double somaI = 0;
        for (int j = 0; j < SL->n; j++)
        {
            somaI = somaI + SL->A[indexMap(i + 1, j + 1, SL->nDiagonais)] * p[j]; // somaI = somaI + A[i][j] * p[j];
        }

        residuo[i] = residuo[i] - alpha * somaI;
        if (isnan(residuo[i]) || isinf(residuo[i]))
        {
            printf("Erro (calcResiduo): %g é NaN ou +/-Infinito, Linha: %i\n", residuo[i], i);
            //     exit(1);
        }
    }
}

/***********************METODOS pre condicionadores*****************************************/

/**
 * @brief Inicializa a matriz M com a diagonal principal invertida de A
 * @param SL - Sistema Linear
 * @param M - Matriz pre condicionadora a ser armazenada a 1/ diagonal
 */
static inline void inicializaPreCondJacobi(SistLinear_t *SL, double *M)
{
    for (int i = 0; i < SL->n; ++i)
    {
        M[i] = (double)1 / SL->A[(indexMap(i + 1, i + 1, SL->nDiagonais))]; // gera vetor com diagonais principal do SL
    }
}

/**
 * @brief - Função para aplicar o pre condicionamento no sistema linear
 * @param SL - O sistema linear que vai ser aplicado o pre condicionamento
 * @param M - A matriz M^-1 a ser aplicada
 */
static inline void aplicaPreCondicSL(SistLinear_t *SL, double *M)
{
    SL->A[0] = 0.0;
    for (unsigned int i = 1; i <= SL->n; ++i)
    {
        SL->A[indexMap(i, i, SL->nDiagonais)] = SL->A[indexMap(i, i, SL->nDiagonais)] * M[i - 1];
        SL->b[i - 1] = SL->b[i - 1] * M[i - 1];
    }
}

/**
 * @brief Calcula a proxima direcao de busca p<k>
 *
 * p<k> = z<k> + beta<k-1> * p<k-1>
 *
 * @param proxDir - vetor da proxima direcao de busca
 * @param z - z calculado
 * @param beta - Beta calculado beta<i-1>
 * @param direcAnterior - Direcao anterior calculada
 * @return double
 */
static inline void calcProxDirecBusca(double *restrict proxDir, double *restrict z, double beta, double *restrict direcAnterior, int n)
{
    for (int i = 0; i < n; ++i)
    {
        proxDir[i] = beta * proxDir[i] + z[i];
    }
}

/***********************METODOS AUXILIARES CG PRE CONDICIONADOR*****************************************/

/**
 * @brief - Calcula proximo alpha
 *
 * alpha = r<k>^t * z<k> / (p<k>^t * A * p<k>)
 *
 * @param resid - residuo calculado
 * @param A - Matriz original
 * @param p - Matriz de direção de busca calculada
 * @param z - Matriz z
 * @param n - dimensão da matriz
 * @return double
 */
static inline double calcAlphaPreCond(double *restrict resid, SistLinear_t *SL, double *restrict p, double *restrict z)
{
    double alpha = 0;
    double *pTxA = (double *)malloc(sizeof(double) * SL->n); // vetor aux para calculo de p^T * A

    double residTxZ = multiplicaVetores(z, resid, SL->n);

    // Percorre a matriz SL->A e o vetor D
    for (int i = 0; i < SL->n; ++i)
    {
        double soma = 0.0;
        for (int j = 0; j < SL->n; ++j)
        {
            // calcula o iésimo elemento do vetor_resultante (esse vetor é chamado de 'z' no livro M.Cristina C. Cunha)
            soma = soma + SL->A[(indexMap(i + 1, j + 1, SL->nDiagonais))] * p[j]; //  soma + A[i][j] * p[j];
        }
        // O iésimo elemento do vetor resultante recebe soma.
        pTxA[i] = soma;

        if (isnan(pTxA[i]) || isinf(pTxA[i]))
        {
            printf("Erro (calcAlphaPreCond pTxA): %g é NaN ou +/-Infinito, Linha: %i\n", pTxA[i], i);
            //     exit(1);
        }
    }

    // Tendo o vetor resultante, agora é só multiplicar D^t * vetor_resultante.
    // O cálculo entre um vetor transposto e um vetor normal gera uma matriz de 1x1, ou seja, um número real.
    double pTxAxP = 0.0;
    for (int i = 0; i < SL->n; ++i)
    {
        pTxAxP = pTxAxP + p[i] * pTxA[i];
        if (isnan(pTxAxP) || isinf(pTxAxP))
        {
            printf("Erro (calcAlphaPreCond pTxAxP): %g é NaN ou +/-Infinito, Linha: %i\n", pTxAxP, i);
            //     exit(1);
        }
    }

    alpha = residTxZ / pTxAxP;
    free(pTxA);

    return alpha;
}

/**
 * @brief Calcula beta com base no residuo atual e no anterior
 *
 * beta<k> = (res<k+1>^T * z<k+1>)/(res<k>^T * z<k>)
 *
 * @param resid - vetor de residuo
 * @param residAnt - vetor de residuo anterior
 * @param z - vetor z
 * @param zAnt - vetor z anterior
 * @param n - tamanho do SL
 * @return double
 */
static inline double calcBetaPreCond(double *restrict resid, double *restrict residAnt, double *restrict z, double *restrict zAnt, int n)
{
    double beta = 0;
    double residTxZ = multiplicaVetores(z, resid, n);
    double zAntxResidAnt = multiplicaVetores(zAnt, residAnt, n);
    beta = residTxZ / zAntxResidAnt;

    if (isnan(beta) || isinf(beta))
    {
        printf("Erro (calcBetaPreCond): %g é NaN ou +/-Infinito\n", beta);
        //     exit(1);
    }
    return beta;
}

/**
 * @brief - Calcula o próximo Z
 *
 * z = M^-1 * residuo
 *
 * @param z - z a ser calculado
 * @param matPreConj
 * @param residuo
 * @param n - Tamanho do vetor
 */
static inline void calcZ(double *z, SistLinear_t *SL, double *residuo, unsigned int n)
{
    for(unsigned int a = 0; a < (1 + (n ) - ((n + 1) % UNROLL)); a+=UNROLL)	//! z = Av
    {
        for(unsigned int b = 0; b < (n); b++)
        {
            for(unsigned int c = 0; c < UNROLL; ++c)
            {
                z[a + c] += SL->A[indexMap(a+c+1,b+1,SL->nDiagonais)] * residuo[b];
            }
        }
    }

    //!remaider loop
    for(unsigned int a = (1 + (n) - ((n + 1) % UNROLL)); a < n; ++a)
    {
        for(unsigned int b = 0; b < (n); ++b)
        {
            z[a] += SL->A[indexMap(a+1,b+1,SL->nDiagonais)] * residuo[b];
        }
    }
}

/******************FUNCOES GRADIENTE CONJUGADO com pre condicionador **********************************/
//// TODO : OTIMIZAR
/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado com uso do pre condicionador de jacobi
 *
 * @param SL - Sistema Linear
 * @param maxIt - numero maximo de iterações
 * @param tol - Tolerancia maxima
 * @param matSaida -
 * @param arqSaida - Arquivo para printar as informações
 */
void gradienteConjugadoPreCondic(SistLinear_t *SL, int maxIt, double tol, FILE *arqSaida)
{
    //! Variáveis:
	double alpha;			//! s
	double beta;			//! m
	double normx;			//! ||x||
	double relerr;			//! erro relativo atual
	double aux0, aux1;		//! aux, aux1
	double *xAnt;		//! vetor x anterior
	double *v;			//! v
	double *z;			    //! z
	double *y;			    //! y
	unsigned int numiter;	//! numero de iteracoes
	unsigned int indxmax;	//! indice no qual max(|xatual - xold|)
	double soma;			//! variavel auxiliar nos lacos
	//!

    double *simmat = SL->A;
    unsigned int n =SL->n;
    double *atvb = SL->b;
    double *res = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));;
    double *vetx =aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	double k = SL->nDiagonais;

    indxmax = 0;

	//! Algoritmo:
	z = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	xAnt = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	memcpy(xAnt, vetx, (n + 1)*sizeof(double));	//! x0 = 0
	memcpy(res, atvb, (n + 1)*sizeof(double));		//! r = b
	y = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	v = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	
	aux0 = 0.0;
	for(unsigned int i = 1; i < (n + 1); i++)		//! v = M-1b, y = M-1r 
	{
		v[i] = (atvb[i] / simmat[indexMap(i,i,k)]);
		y[i] = (res[i] / simmat[indexMap(i,i,k)]);
		aux0 += y[i] * res[i];	//! aux = ytr
	}
	// LIKWID_MARKER_START("op1");
	
	numiter = 0;	
	while(numiter < maxIt) 									//! para k = 0 : max, faca
	{
		
		soma =  0.0;
		//!Unroll and Jam z = Av e soma =  vtz
		memset(z, 0, (n + 1));
		for(unsigned int i =  1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	//! z = Av
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
				soma += z[i + c] * v[i + c];	//! calcula vtz
			}
		}

		//!remaider loop
		for(unsigned int i =  (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)
		{
			for(unsigned int j = 1; j < (n + 1); ++j)
			{
				z[i] += simmat[indexMap(i,j,k)] * v[j];
			}
			soma += z[i] * v[i];
		}

		//! s = aux/vtz
		alpha =  (fabs(soma) < ZERO) ? 0.0 : (aux0 / soma);

		normx = 0.0;
		// CALCULO DA NORMA
		//!Unroll and Jam do calculo da norma e x = xold
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				vetx[i + j] += (alpha * v[i + j]);			//! xk+1 = xk + sv
				if (normx < fabs(vetx[i + j] - xAnt[i + j]))	//! calcula ||x||
				{
					normx = fabs(vetx[i + j] - xAnt[i + j]);
					indxmax = i+j; 
				}
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//!remainder loop
		{
			vetx[i] += (alpha * v[i]);			//! xk+1 = xk + sv
			if (normx < fabs(vetx[i] - xAnt[i]))	//! calcula ||x||
			{
				normx = fabs(vetx[i] - xAnt[i]);
				indxmax = i; 
			}  
		}
		
		//!Unroll and Jam de r = r - sz e y = M-1r
		// CALCULO DO RESID E Y
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				res[i + j] -= (alpha * z[i + j]);
				y[i + j] = (res[i + j] / simmat[indexMap(i+j,i+j,k)]);
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//!remainder loop
		{
			res[i] -= (alpha * z[i]);
			y[i] = (res[i] / simmat[indexMap(i,i,k)]);  
		}
		aux1 = 0.0;

		//!Unroll and Jam de aux1 = ytr
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				aux1 += y[i + j] * res[i + j];
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//!remainder loop
		{
			aux1 += y[i] * res[i];  
		}
		
		//! m = aux1 / aux
		beta =  (aux1 / aux0);
		aux0 = aux1;

		//!Unroll and Jam de v = y + mv
		for(unsigned int i = 1; i < (1 + (n + 1) - ((n + 1) % STRIDE)); i+=STRIDE)	
		{
			for(unsigned int j = 0; j < STRIDE; ++j)
			{
				v[i + j] = y[i + j] + (beta * v[i + j]);
			} 
		}
		for(unsigned int i = (1 + (n + 1) - ((n + 1) % STRIDE)); i < (n + 1); ++i)	//!remainder loop
		{
			v[i] = y[i] + (beta * v[i]);  
		}	
		fprintf(arqSaida, "# iter %u: %.15g\n", numiter, normx);

		//! xold = x
		memcpy(xAnt, vetx, (n + 1)*sizeof(double));
		//! relerr = max(|Xi - Xi-1|) / Xi
		relerr =(normx / fabs(vetx[indxmax]));
		//!
		++numiter;
	};

	// LIKWID_MARKER_STOP("op1");
	
	free(v);
	free(z);
	free(y);
	free(xAnt);
}
