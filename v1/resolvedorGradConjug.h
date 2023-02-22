void calcX(double *proxX, double *xAnt, double alpha, double *p, int n);
void calcResiduo(double *residuoAnterior, double alpha, double **A, double *p, double *residuo, int n);
void calcProxDirecBusca(double *proxDir, double *z, double beta, double *direcAnterior, int n);
void calculaResiduoOriginal(double**A,  double*b, double *x, int n,double *residuo); 
void inicializaPreCondJacobi(SistLinear_t *SL, double *M);
void aplicaPreCondicSL(SistLinear_t *SL, double *M);
double calcAlphaPreCond(double *resid, double **A, double *p, double *z, int n);
void calcZ(double *z, double *matPreConj, double *residuo, unsigned int n);
double calcBetaPreCond(double *resid, double *residAnt, double *z, double *zAnt, int n);

double calcAlpha(double *resid, double *direc, SistLinear_t *SL);
double calcBeta(double *resid, double *residAnt, unsigned int n);

void gradienteConjugadoPreCondic(SistLinear_t *SL, int maxIt, double tol,  FILE *arqSaida); 
void gradienteConjugado(SistLinear_t *SL, int maxIt, double tol, FILE *arqSaida);