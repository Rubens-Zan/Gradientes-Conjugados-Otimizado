#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include "utils.h"
#include "sislin.h"
#include "math.h"

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

// double timestamp(void)
// {
//   struct timespec tp;
//   clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
//   return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
// }

/**
 * @brief Função para fazer o tratamento da entrada
 * @param argc 
 * @param argv 
 * @param comando 
 */
void tratamentoEntrada(int argc, char **argv, tComando *comando){
    comando->erroMax = 0; // se o erro for 0, vai ser desconsiderado
    comando->nIter = 0;
    comando->dimensao = 0;
    
    for (int i=0;i<argc;++i){
        if (strcmp ( argv[i], "-n") == 0 && (i+1 < argc)){
            comando->dimensao = atoi(argv[i+1]); // dimensao do sl
            i++; 
        }else if (strcmp ( argv[i], "-k") == 0 && (i+1 < argc)){
           comando->nDiagonais = atoi(argv[i+1]); // numero de diagonais
           i++;
        }else if (strcmp ( argv[i], "-p") == 0 && (i+1 < argc)){
            comando->usarPreCondicionador = atoi(argv[i+1]) > 0 ? true : false; // sem precondicionador ou de jacobi 
            i++;
        }else if(strcmp ( argv[i], "-i") == 0 && (i+1 < argc)){
            comando->nIter = atoi(argv[i+1]);  // numero de iteraçoes
            i++;
        }else if(strcmp ( argv[i], "-e") == 0 && (i+1 < argc)){
            comando->erroMax = atof(argv[i+1]);  // erro maximo
            i++;           
        }else if (strcmp ( argv[i], "-o") == 0 && (i+1 < argc)){
           strcpy(comando->saida, argv[i+1]); //arquivo de saida
           i++;
        }  
    }

    if (!comando->dimensao || !comando->nDiagonais || !comando->nIter){
        fprintf (stderr, "Erro: argumento invalido\n");
		exit(1);
    }

}

/******************NORMAS**********************************/

/**
 * @brief Calcula a norma do residuo
 * normaL2Residuo = sqrt(residuo^2)
 * @param residuo - vetor de residuo
 * @param n - tamanho do vetor
 * @return double 
 */
double normaL2Residuo( double *residuo, unsigned int n)
{

    double soma = 0.0;
    double raiz;

    // Pecorre o vetor de soluções
    for (int i = 1; i <= n; ++i) {

        // Soma o quadrado dos elementos das soluções
        soma = soma + residuo[i]*residuo[i];

    }
    raiz = sqrt(soma);

    // Retorna a raíz quadrada da soma.
    return raiz;
}

/**
 * @brief - Calcula a norma max da solução relativa
 * max (|(x - xAnt)| / |x|)
 * @param x - Vetor com a solução atual
 * @param xAnt - Vetor com a solução antiga
 * @param n - tamanho do vetor
 * @return double 
 */
double normaMaxRelat( double  * restrict x,  double  * restrict xAnt, unsigned int n,  double  * restrict maiorErroAbs)
{

    double maiorErro = ABS(x[0] - xAnt[0]) / ABS(x[0]);
    *maiorErroAbs = ABS(x[0] - xAnt[0]);
    for (int i = 1; i < n; ++i)
    {
        if (ABS(x[i] - xAnt[i]) / ABS(x[i]) > maiorErro)
        {
            maiorErro = ABS(x[i] - xAnt[i]) / ABS(x[i]);      

            if (ABS(x[i] - xAnt[i] ) > *maiorErroAbs)
                *maiorErroAbs = ABS(x[i] - xAnt[i]);
        }
    }

    // Retorna ao final o maiorErro erro absoluto.
    return maiorErro;
}

