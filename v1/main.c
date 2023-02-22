#include <stdio.h>
#include <stdlib.h>
#include <likwid.h>
#include "utils.h"
#include "sislin.h"
#include "resolvedorGradConjug.h"


int main(int argc, char **argv)
{
    LIKWID_MARKER_INIT;

    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao);
    // SistLinear_t *SL = alocaSisLin(2);
    FILE *arqSaida;
    
    srand(20222);

	arqSaida = fopen(comando->saida,"w+");
    
    
	fprintf(arqSaida,"#rzl20 Rubens Zandomenighi Laszlo \n");
	fprintf(arqSaida,"# \n");


    iniSisLin(SL, comando->nDiagonais);

    if (comando->usarPreCondicionador){
        gradienteConjugadoPreCondic(SL, comando->nIter,comando->erroMax,arqSaida);
    } else {
        gradienteConjugado(SL,comando->nIter,comando->erroMax, arqSaida); 
    }

	fclose(arqSaida);
    free(comando);
    liberaSisLin(SL);
    LIKWID_MARKER_CLOSE;
    return 0;
}