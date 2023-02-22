#include "sislin.h"
#include "utils.h"
#include "resolvedorGradConjug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <likwid.h>



/**
 * @brief Função para calcular o residuo original conforme A e b originais
 * 
 * @param A - Coeficientes A originais
 * @param b - Coeficientes B originais
 * @param x - Solução final
 * @param n - Tamanho do Sistema Linear
 * @param residuo - Vetor a receber o residuo
 */
void calculaResiduoOriginal(double**A,  double*b, double *x, int n,double *residuo)
{

    // Percorre a matriz,  
    for (int i = 0; i < n; ++i)
    {
        double soma = 0.0;
        for (int j = 0; j < n; ++j)
        {
            // realiza a soma das linhas, 
            soma = soma + A[i][j] * x[j];

            // (Teste para ver se não foi gerado um NaN ou um número infinito)
            // if (isnan(soma) || isinf(soma))
            // {
            //     fprintf(stderr, "Erro soma(calculaResiduoOriginal): %g é NaN ou +/-Infinito\n", soma);
            //     exit(1);
            // } 
        }
        //e tira do resíduo.
        residuo[i] = b[i] - soma;
    }
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
void calcX(double *proxX, double *xAnt, double alpha, double *p, int n)
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
void calcResiduo(double *residuoAnterior, double alpha, double **A, double *p, double *residuo, int n)
{
    // alpha * A * p<k>
    for (int i = 0; i < n; ++i)
    {
        double somaI = 0;
        for (int j = 0; j < n; j++)
        {
            somaI = somaI + A[i][j] * p[j];
        }

        residuo[i] = residuo[i] - alpha * somaI;
    }
}

/***********************METODOS pre condicionadores*****************************************/

/**
 * @brief Inicializa a matriz M com a diagonal principal invertida de A
 * @param SL - Sistema Linear
 * @param M - Matriz pre condicionadora a ser armazenada a 1/ diagonal 
 */
void inicializaPreCondJacobi(SistLinear_t *SL, double *M)
{
    for (int i = 0; i < SL->n; ++i)
    {
        M[i] = (double)1 / SL->A[i][i]; // gera vetor com diagonais do SL
        // if (isnan(M[i]) || isinf(M[i]))
        // {
        //     fprintf(stderr, "Erro M[i](inicializaPreCondJacobi): %g é NaN ou +/-Infinito, Linha: %i\n", M[i], i);
        //     exit(1);
        // }
    }
}

/**
 * @brief - Função para aplicar o pre condicionamento no sistema linear
 * @param SL - O sistema linear que vai ser aplicado o pre condicionamento
 * @param M - A matriz M^-1 a ser aplicada
 */
void aplicaPreCondicSL(SistLinear_t *SL, double *M)
{
    for (int i = 0; i < SL->n; i++)
    {
        SL->b[i] = SL->b[i] * M[i];
        for (int j = 0; j < SL->n; j++)
        {
            SL->A[i][j] = SL->A[i][j] * M[i]; // aplica pre condicionador se for na diagonal
        }
    }
}

/**
 * @brief Inicializa o chute inicial com 0
 * @param x - vetor solução
 * @param n - Tamanho do SL
 */
void inicializaSol(double *x, unsigned int n)
{
    for (int i = 0; i < n; ++i)
        x[i] = 0;
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
double calcAlphaPreCond(double *resid, double **A, double *p, double *z, int n)
{
    double alpha = 0;
    double *pTxA = (double *)malloc(sizeof(double) * n); // vetor aux para calculo de p^T * A

    double residTxZ = multiplicaVetores(z, resid, n);


    //Percorre a matriz SL->A e o vetor D
    for (int i = 0; i < n; ++i)  
    {
        double soma = 0.0;
        for (int j = 0; j < n; ++j)
        {
            //calcula o iésimo elemento do vetor_resultante (esse vetor é chamado de 'z' no livro M.Cristina C. Cunha)
            soma = soma + A[i][j] * p[j];

            // Teste para ver se não foi gerado um NaN ou um número infinito
            // if (isnan(soma) || isinf(soma))
            // {
            //     fprintf(stderr, "Erro soma(calcularDenominadorEscalarA): %g é NaN ou +/-Infinito\n", soma);
            //     exit(1);
            // }
        }
        // O iésimo elemento do vetor resultante recebe soma.
        pTxA[i] = soma;
    }

    // Tendo o vetor resultante, agora é só multiplicar D^t * vetor_resultante.
    // O cálculo entre um vetor transposto e um vetor normal gera uma matriz de 1x1, ou seja, um número real.  
    double pTxAxP = 0.0;
    for (int i = 0; i < n; ++i)
    {
        pTxAxP = pTxAxP + p[i] * pTxA[i];
        // Teste para ver se não foi gerado um NaN ou um número infinito
        // if (isnan(pTxAxP) || isinf(pTxAxP)){
        //     fprintf(stderr, "Erro pTxAxP(calcAlphaPreCond): %g é NaN ou +/-Infinito\n", pTxAxP);
        //     exit(1);
        // }
    }

    alpha = residTxZ / pTxAxP;
    // Verificação se resultou em NaN ou +/- infinito
    // if (isnan(alpha) || isinf(alpha))
    // {
    //     fprintf(stderr, "Erro variavel invalida: alpha(calcAlphaPreCond): %g é NaN ou +/-Infinito\n", alpha);
    //     exit(1);
    // }
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
double calcBetaPreCond(double *resid, double *residAnt, double *z, double *zAnt, int n)
{
    double beta = 0;
    double residTxZ = multiplicaVetores(z, resid, n);
    double zAntxResidAnt = multiplicaVetores(zAnt, residAnt, n);

    beta = residTxZ / zAntxResidAnt;
    // Verificação se resultou em NaN ou +/- infinito
    // if (isnan(beta) || isinf(beta))
    // {
    //     fprintf(stderr, "Erro variavel invalida: beta(calcbetaPreCond): %g é NaN ou +/-Infinito\n", beta);
    //     exit(1);
    // }
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
void calcZ(double *z, double *matPreConj, double *residuo, unsigned int n)
{
    for (int i = 0; i < n; ++i)
    {
        z[i] = matPreConj[i] * residuo[i];
    }
}

/******************FUNCOES GRADIENTE CONJUGADO com pre condicionador **********************************/

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
    double **Aoriginal = alocarMatriz(SL->n, SL->n); // Matriz original
    double *boriginal = (double *)malloc(sizeof(double) * SL->n);  // coeficientes b originais

    double **novoA = alocarMatriz(SL->n, SL->n); // Matriz original
    double *novoB = (double *)malloc(sizeof(double) * SL->n);  // coeficientes b originais

    double alpha, beta;
    double *resid = (double *)malloc(sizeof(double) * SL->n);    // matriz de residuo
    double *residAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de residuo anterior
    double *direc = (double *)malloc(sizeof(double) * SL->n);    // matriz de direcao de busca
    double *direcAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de direcao anterior
    double *xAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de chute anterior
    double *z = (double *)malloc(sizeof(double) * SL->n);        // matriz de z
    double *zAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de z antigo
    double *x = (double *)malloc(sizeof(double) * SL->n);        // matriz de soluções
    double *matJacobiInvert = (double *)malloc(sizeof(double) * SL->n);

    double tMedioIter, tempoResid, tempoPreCond, tempoInicio;
    tMedioIter = 0;
    tempoInicio = timestamp();

    int it;

    /**/
    inicializaPreCondJacobi(SL, matJacobiInvert);


    copiaMatriz(SL->A,Aoriginal, SL->n);
    copiaVetor(SL->b, boriginal, SL->n); 

    calcularAtxB(SL,novoB);
    calcularMatrizAtxA(SL, novoA);

    copiaMatriz(novoA,SL->A, SL->n);
    copiaVetor(novoB, SL->b, SL->n); 

    
    // (M^-1) * A
    // (M^-1) * b
    aplicaPreCondicSL(SL, matJacobiInvert);

    tempoPreCond = timestamp() - tempoInicio;
    inicializaSol(x, SL->n);
    // residuo = b
    copiaVetor(SL->b, resid, SL->n); // como x inicial igual a 0, desconsidero o r = b - (A * x)
    // calcula z
    calcZ(z, matJacobiInvert, resid, SL->n);
    // direc = z
    copiaVetor(z, direc, SL->n);

    LIKWID_MARKER_START ("op1");
    for (it = 0; it < maxIt; ++it)
    {
        double tIterInicio = timestamp();

        // ALPHA = z * r / pT * A * p
        alpha = calcAlphaPreCond(resid, SL->A, direc, z, SL->n);
        // calcula novo x
        copiaVetor(x, xAnt, SL->n); // xant = x
        calcX(x, xAnt, alpha, direc, SL->n);
        
        double maiorErroAbs; 
        double normaMaxRel = normaMaxRelat(x, xAnt, SL->n, &maiorErroAbs);
        fprintf(arqSaida, "# iter %d: ||%.15g||\n", it, maiorErroAbs);
        // calcula novo residuo
        copiaVetor(resid, residAnt, SL->n);
        calcResiduo(residAnt, alpha, SL->A, direc, resid, SL->n);

        // if erro
        // o erro aproximado em x após a k-ésima iteração (max|xi - xi-1|);
        if (tol != 0 && tol > normaMaxRel)
        {
            tMedioIter += timestamp() - tIterInicio;
            break;
        }
        // z0 = z
        copiaVetor(resid, zAnt, SL->n);
        // calcula z
        calcZ(z, matJacobiInvert, resid, SL->n);
        // beta = rt * z / rtant *zant
        // tb usa
        beta = calcBetaPreCond(resid, residAnt, z, zAnt, SL->n);
        // calcula prox direcao
        copiaVetor(direc, direcAnt, SL->n); // dAnt = d
        calcProxDirecBusca(direc, z, beta, direcAnt, SL->n);

        tMedioIter += timestamp() - tIterInicio;
    }
    LIKWID_MARKER_STOP ("op1");

    // restaura A e b originais
    copiaMatriz(SL->A,Aoriginal, SL->n);
    copiaVetor(SL->b, boriginal, SL->n); 

    tempoResid = timestamp();
    LIKWID_MARKER_START ("op2");
    calculaResiduoOriginal(SL->A, SL->b, x, SL->n, resid);
    LIKWID_MARKER_STOP ("op2");
    
    fprintf(arqSaida, "# residuo: || %.15g || \n", normaL2Residuo(resid, SL->n));
    tempoResid = timestamp() - tempoResid; 

    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tMedioIter / it;
    fprintf(arqSaida, "# Tempo iter: %.15g \n", tMedioIter);

    fprintf(arqSaida, "# Tempo residuo: %.15g \n", tempoResid);
    fprintf(arqSaida, "# \n");
    fprintf(arqSaida, "%d", SL->n);
    prnVetorArq(x, SL->n, arqSaida);

    liberarMatriz(Aoriginal);
    free(boriginal); 
    liberarMatriz(novoA);
    free(novoB); 
    free(resid);
    free(residAnt);
    free(direc);
    free(direcAnt);
    free(xAnt);
    free(x);
    free(z);
    free(zAnt);
    free(matJacobiInvert);
}

/******************FUNCOES AUXILIAR SEM PRE CONDICIONADOR**********************************/

/**
 * @brief  Calcula o alpha
 *
 * alpha = rT*r / direcT * A * direc
 *
 * @param resid - vetor com residuo
 * @param direc - Direcao de busca atual
 * @param SL - Sistema linear
 * @return double
 */
double calcAlpha(double *resid, double *direc, SistLinear_t *SL)
{
    double residTxresid = multiplicaVetores(resid, resid, SL->n);
    double alpha = 0;
    double *vetDirecxA = (double *)malloc(sizeof(double) * SL->n);
    double direcTxAxDirec = 0.0;

    for (int i = 0; i < SL->n; ++i)
    {
        double totalIDirecXa = 0.0;
        for (int j = 0; j < SL->n; ++j)
        {
            totalIDirecXa = totalIDirecXa + SL->A[i][j] * direc[j];

            // Testa valores inválidos
            // if (isnan(totalIDirecXa) || isinf(totalIDirecXa))
            // {
            //     fprintf(stderr, "Erro variavel invalida: totalIDirecXa(calcAlpha): %g é NaN ou +/-Infinito\n", totalIDirecXa);
            //     exit(1);
            // }
        }
        vetDirecxA[i] = totalIDirecXa;
    }

    for (int i = 0; i < SL->n; ++i)
    {
        direcTxAxDirec += direc[i] * vetDirecxA[i];
    }

    alpha = residTxresid / direcTxAxDirec;
    // Testa valores inválidos

    // if (isnan(direcTxAxDirec) || isinf(direcTxAxDirec))
    // {
    //     fprintf(stderr, "Erro variavel invalida: direcTxAxDirec(calcAlpha): %g é NaN ou +/-Infinito\n", direcTxAxDirec);
    //     exit(1);
    // }
    free(vetDirecxA);
    return alpha;
}

/**
 * @brief - Calcula o valor de beta
 * 
 * beta = r^T * r / rAnt^T * rAnt
 * 
 * @param resid - vetor de residuo
 * @param residAnt - vetor de residuo anterior
 * @param n - Tamanho do SL
 * @return double
 */
double calcBeta(double *resid, double *residAnt, unsigned int n)
{
    double beta = multiplicaVetores(resid, resid, n) / multiplicaVetores(residAnt, residAnt, n);
    return beta;
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
void calcProxDirecBusca(double *proxDir, double *z, double beta, double *direcAnterior, int n)
{
    for (int i = 0; i < n; ++i)
    {
        proxDir[i] = beta * proxDir[i] + z[i] ;

    }
}

/******************FUNCAO GRADIENTE CONJUGADO SEM PRE CONDICIONADOR **********************************/

/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado
 *
 * @param SL - Sistema Linear
 * @param maxIt - numero maximo de iterações
 * @param tol - Tolerancia maxima
 * @param matSaida -
 * @param arqSaida - Arquivo para printar as informações
 */
void gradienteConjugado(SistLinear_t *SL, int maxIt, double tol, FILE *arqSaida)
{
    double **Aoriginal = alocarMatriz(SL->n, SL->n); // Matriz original
    double *boriginal = (double *)malloc(sizeof(double) * SL->n);  // coeficientes b originais

    double **novoA = alocarMatriz(SL->n, SL->n); // Matriz original
    double *novoB = (double *)malloc(sizeof(double) * SL->n);  // coeficientes b originais
    
    double tMedioIter, tempoResid, tempoPreCond, tempoInicio;
    double alpha, beta;
    double *resid = (double *)malloc(sizeof(double) * SL->n);    // matriz de residuo
    double *residAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de residuo anterior
    double *direc = (double *)malloc(sizeof(double) * SL->n);    // matriz de direcao de busca
    double *direcAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de direcao anterior
    double *xAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de chute anterior
    double *x = (double *)malloc(sizeof(double) * SL->n);
    int it = 0;
    tMedioIter = 0;
    
    tempoInicio = timestamp();

    copiaMatriz(SL->A,Aoriginal, SL->n);
    copiaVetor(SL->b, boriginal, SL->n); 

    calcularAtxB(SL,novoB);
    calcularMatrizAtxA(SL, novoA);

    copiaMatriz(novoA,SL->A, SL->n);
    copiaVetor(novoB, SL->b, SL->n); 

    // nao usa pre condicionadores
    tempoPreCond = timestamp() - tempoInicio;

    inicializaSol(x, SL->n);
    // residuo = b
    copiaVetor(SL->b, resid, SL->n);
    // direcao = residuo
    copiaVetor(resid, direc, SL->n);

    for (it = 0; it < maxIt; ++it)
    {
        double tIterInicio = timestamp();

        // calcula alpha
        alpha = calcAlpha(resid, direc, SL);
        // calcula novo x
        // x1 = x0 +alpha * p
        copiaVetor(x, xAnt, SL->n);
        calcX(x, xAnt, alpha, direc, SL->n);

        double maiorErroAbs;
        double normaMaxRel = normaMaxRelat(x, xAnt, SL->n, &maiorErroAbs);
        // TODO: AJUSTAR PARA PRINTAR A MAIOR DIFERENÇA ENTRE OS X
        fprintf(arqSaida, "# iter %d: ||%.15g|| \n", it, maiorErroAbs);

        // calcular novo residuo
        // r1= r0 - alpha0 * A *p0
        copiaVetor(resid, residAnt, SL->n);
        calcResiduo(residAnt, alpha, SL->A, direc, resid, SL->n);

        // normaMaxRelat
        // if
        if (tol != 0 && tol > normaMaxRel)
        {
            tMedioIter += timestamp() - tIterInicio;
            break;
        }

        // escalarBeta = vetor residuo^t * vetor residuo / vetor residuo anteriorT * vetor residuo anterior
        beta = calcBeta(resid, residAnt, SL->n);

        // calcula prox direc
        // P = R + b0 * pAnt
        copiaVetor(direc, direcAnt, SL->n);
        calcProxDirecBusca(direc, resid, beta, direcAnt, SL->n);

        tMedioIter += timestamp() - tIterInicio;
    }

    // restaura A e b originais
    copiaMatriz(SL->A,Aoriginal, SL->n);
    copiaVetor(SL->b, boriginal, SL->n); 

    tempoResid = timestamp();
    calculaResiduoOriginal(SL->A, SL->b, x, SL->n, resid);
    fprintf(arqSaida, "# residuo: || %.15g || \n", normaL2Residuo(resid, SL->n));
    tempoResid = timestamp() - tempoResid; 

    // tempo final
    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tMedioIter / it;
    fprintf(arqSaida, "# Tempo iter: %.15g \n", tMedioIter);

    fprintf(arqSaida, "# Tempo residuo: %.15g \n", tempoResid);
    fprintf(arqSaida, "# \n");
    fprintf(arqSaida, "%d", SL->n);
    prnVetorArq(x, SL->n, arqSaida);
    
    liberarMatriz(Aoriginal);
    free(boriginal); 
    liberarMatriz(novoA);
    free(novoB); 
    free(resid);
    free(residAnt);
    free(direc);
    free(direcAnt);
    free(xAnt);
    free(x);
}
