#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density;

  DATA_VECTOR        (N);        /* Observations (counts) */
  DATA_FACTOR        (position); /* Observations */
  DATA_FACTOR        (time);     /* Observations */
  DATA_FACTOR        (haulid);   /* Observations */

  DATA_SPARSE_MATRIX (Q0);       /* For random field */
  DATA_SPARSE_MATRIX (I);        /* For random field */
  DATA_SPARSE_MATRIX (A);        /* Design matrix (standardization) */

  // Fixed effects
  PARAMETER          (logdelta); /* For random field (corr) */
  PARAMETER          (logkappa); /* For random field (scale) */
  PARAMETER          (tphi_time);/* One-step time correlation */
  PARAMETER          (logsigma); /* Nugget error */
  PARAMETER_VECTOR   (beta);     /* For design matrix */

  // Random effects
  PARAMETER_ARRAY    (eta);    // 2D: space x time
  PARAMETER_VECTOR   (etanug); // 1D: haul

  // Parameter transforms
  Type delta = exp(logdelta);
  Type kappa = exp(logkappa);
  Type sigma = exp(logsigma);
  Type phi_time = tphi_time / sqrt( 1.0 + tphi_time * tphi_time );

  // Random fields
  SparseMatrix<Type> Q = kappa * (Q0 + delta * I);
  GMRF_t<Type>      nldens_spatial = GMRF(Q);
  AR1_t<N01<Type> > nldens_time    = AR1(phi_time);

  Type nll = 0; // Negative log likelhood

  // Process likelihood
  nll += SEPARABLE(nldens_time, nldens_spatial)(eta);

  // Nugget
  nll -= dnorm(etanug, Type(0), sigma, true).sum();

  // Measurement likelihood
  vector<Type> predictor = A * beta;
  for(int i=0; i<N.size(); i++){
    nll -= dpois(N(i), 
		 exp(predictor(i) +
		     eta(position(i), time(i)) +
		     etanug(haulid(i))
		     ),
		 true);
  }

  ADREPORT(phi_time);

  return nll;

}
