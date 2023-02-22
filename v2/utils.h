#include <stdbool.h>

typedef struct tComando {
    int dimensao,nDiagonais, nIter;
    double erroMax;
    bool usarPreCondicionador;   
    char saida[100];
} tComando;
double normaMaxRelat( double * restrict x,  double * restrict xAnt, unsigned int n,  double * restrict maiorErroAbs); 

double normaL2Residuo( double *residuo, unsigned int n);

double timestamp(void);
void tratamentoEntrada(int argc, char **argv, tComando *comando);

