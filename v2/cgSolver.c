#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "resolvedorGradConjug.h"


int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao, comando->nDiagonais);
	unsigned int simk = (2*SL->nDiagonais) - 1;
    SistLinear_t *novoSL = alocaSisLin(comando->dimensao, simk);
   
    FILE *arqSaida;
    
    srand(20222);

	arqSaida = fopen(comando->saida,"w+");
    
    
	fprintf(arqSaida,"#rzl20 Rubens Zandomenighi Laszlo \n");
	fprintf(arqSaida,"# \n");

    double *matT = aligned_alloc(ALIGNMENT, ((SL->n * SL->nDiagonais * sizeof(double)) + (ALIGNMENT - ((SL->n * SL->nDiagonais * sizeof(double)) % ALIGNMENT))));

    iniSisLin(SL, comando->nDiagonais,matT);

    calcularMatrizAtxA(SL,matT,novoSL);
    calcularAtxB(SL->nDiagonais, SL->n, matT, SL->b, novoSL->b); 
    
    if (comando->usarPreCondicionador){
        gradienteConjugadoPreCondic(novoSL, comando->nIter,comando->erroMax,arqSaida);
    } else {
        printf("ERRO: DEVE SER USADO APENAS COM PRE CONDICIONADOR"); 
    }
	fclose(arqSaida);
    free(comando);
    // liberaSisLin(SL);
    return 0;
}