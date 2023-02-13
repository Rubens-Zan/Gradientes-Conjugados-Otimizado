#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"


int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao, comando->nDiagonais);
    // FILE *arqSaida;
    
    srand(20222);

	// arqSaida = fopen(comando->saida,"w+");
    
    
	// fprintf(arqSaida,"#rzl20 Rubens Zandomenighi Laszlo \n");
	// fprintf(arqSaida,"# \n");


    iniSisLin(SL, comando->nDiagonais);
    // if (comando->usarPreCondicionador){
    //     gradienteConjugadoPreCondic(SL, comando->nIter,comando->erroMax,arqSaida);
    // } else {
    //     gradienteConjugado(SL,comando->nIter,comando->erroMax, arqSaida); 
    // }
	// fclose(arqSaida);
    free(comando);
    liberaSisLin(SL);
    return 0;
}