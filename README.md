# Trabalho 02: Otimização do Desempenho para Sistemas Lineares Esparsos
Rubens Zandomenighi Laszlo GRR20206147

## Na operação de iteração do método de Gradiente Conjugado com Pré-condicionador de Jacobi (op1);

Na inicialização da solução inicial X, alteração da função que efetua a inicialização do vetor com 0 para utilização da função memset. 


## Na operação de cálculo do resíduo (op2): R = b - Ax;

Para uma melhor utilização da cache foi utilizado a técnica de unroll e Jam, sendo utilizado o fator de unroll = 4.


## Geral
Alteração da estrutura do Sistema Linear para o armazenamento apenas das diagonais não nulas da matriz A. 
Sendo feita a substituição da matriz double **A para o vetor 'double *A' para armazenamento dessas, sendo as primeiras n posições alocadas para a diagonal principal (sendo que essa sempre é não nula nesse projeto), e as matrizes subsequentes, posteriormento as diagonais inferiores e depois as superiores.

Alocação apenas das matrizes não nulas 
SL->A = (double *)malloc((n +1)+ (2 * calcElementoDiagonaisInf(n,k))  * sizeof(double));
Tal como: 
A função: calcElementoDiagonaisInf, utiliza o fato que ,caso a diagonal esteja dentro do limite das k-diagonais da matriz de banda, o número de elementos dessa será *n-distancia_diagonal_principal*, como por exemplo:
Caso o número de diagonais seja 5 e a matriz seja de 10 elementos. 
Se estou na posição caso eu esteja na 1a diagonal inferior, ou seja a distância para a diagonal principal será de 1 e a quantidade de elementos nessa diagonal será de 10-1=9, e assim consecutivamente, enquanto *k < (nDiagonais /2)*, pois são divididas entre superiores e inferiores. 
Fusão dos laços na função iniSisLin para a geração das diagonais inferiores e superiores. 

Para a facilitação da utilização de instruções SIMD, foi utilizado a palavra chave _restrict_, assim informando o compilador que não existe dependência de dados, em funções em que parâmetros utilizam ponteiros para o mesmo tipo de dados. 
