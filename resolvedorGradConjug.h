// void calcX( double  * restrict proxX,  double  * restrict xAnt, double alpha,  double  * restrict p, int n);
// void calcProxDirecBusca(double *proxDir, double *z, double beta, double *direcAnterior, int n);
// void inicializaPreCondJacobi(SistLinear_t *SL, double *M);
// void aplicaPreCondicSL(SistLinear_t *SL, double *M);
// double calcAlphaPreCond( double  * restrict resid, SistLinear_t *SL,  double  * restrict p,  double  * restrict z);
// void calculaResiduoOriginal(SistLinear_t *SL,  double * restrictx, double * restrict residuo);

// void calcZ( double  * restrict z,  double  * restrict matPreConj,  double  * restrict residuo, unsigned int n);
// double calcBetaPreCond( double  * restrict resid,  double  * restrict residAnt,  double  * restrict z,  double  * restrict zAnt, int n);
// void calcResiduo( double  * restrict residuoAnterior, double alpha,SistLinear_t *SL,  double  * restrict p,  double  * restrict residuo);

static inline void calcZ(double *z, SistLinear_t *SL, double *residuo, unsigned int n);

void gradienteConjugadoPreCondic(SistLinear_t *SL, int maxIt, double tol,  FILE *arqSaida); 