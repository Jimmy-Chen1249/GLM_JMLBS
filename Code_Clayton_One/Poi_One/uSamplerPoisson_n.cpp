/** 
 * \file uSamplerPoisson_n.cpp
 * \author Felipe Acosta
 * \date 2014-12-30
 * \brief This function performs an MCMC run on the random effects. The arguments are the same arguments used in
 * loglikelihood with an extra argument 'B' which indicates the MCMC sample size and the argument 'u' now 
 * indicates the intial value for the chain.
 */


#include "RcppArmadillo.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double min0(double a, double b) {
  if (a < b)
    return a;
  return b;
}

// [[Rcpp::export]]
/** Evaluate the log density of a multivariate normal distribution with mean vector 0 */
double ldmn(const arma::vec& x, const arma::mat& sigma) {
  int kDIM = sigma.n_cols;
  for (int i = 0; i < kDIM; i++) {
    for (int j = 0; j < kDIM; j++) {
      if(sigma(i, j) < 0) {
        return -INFINITY;
      }
    }
  }
  
  double VALUE = -0.5 * kDIM * log(2 * M_PI) - 0.5 * log(arma::det(sigma));
  
  arma::mat sigmainv;
  sigmainv = inv(sigma);
  //std::cout<<sigmainv(1,1);
  
  double tmp0 = 0;
  NumericVector tmpVector(kDIM); /** stores the product of x^t and sigma^-1 */
for (int i = 0; i < kDIM; i++) {
  for (int j = 0; j < kDIM; j++) {
    tmpVector(i) += x(j) * sigmainv(j, i);
  }
}
for (int i = 0; i < kDIM; i++) {
  tmp0 += tmpVector(i) * x(i);
}
VALUE += - 0.5 * tmp0;
return VALUE;

}

// [[Rcpp::export]]
/** Evaluate the log density of a multivariate t distribution with mean vector 0*/
double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType) {
  if (df <= 0) {
    return -INFINITY;
  }
  int n = sigma.n_cols;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (sigma(i, j) < 0) {
        return -INFINITY;
      }
    }
  }
  double value = 0; /** Density value */
int k = x.size();
//arma::mat sigma0 = as<arma::mat>(sigma);
value += lgamma(0.5 * (df + k)) - lgamma(0.5 * df) - 0.5 * k * log(df) - 0.5 * k * log(M_PI) - 0.5 * log(arma::det(sigma));

/** If sigma is diagonal it's easier to do the multiplication */
double tmp0 = 0;
if (sigmaType == 0) {
  for (int i = 0; i < k; i++) {
    tmp0 += x(i) * x(i) / sigma(i,i);
  }
  value += - 0.5 * (df + k) * log(1 + tmp0/df);
  return value;
}

arma::mat sigmainv;
sigmainv = inv_sympd(sigma);
//std::cout<<sigmainv(1,1);
NumericVector tmpVector(k); /** stores the product of x^t and sigma^-1 */
for (int i = 0; i < k; i++) {
  for (int j = 0; j < k; j++) {
    tmpVector(i) += x(j) * sigmainv(j, i);
  }
}
for (int i = 0; i < k; i++) {
  tmp0 += tmpVector(i) * x(i);
}
value += - 0.5 * (df + k) * log(1 + tmp0/df);
return value;
}

// [[Rcpp::export]]
double loglikelihoodPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  double value = 0; /** The value to be returned */

int nObs = kY.n_elem;
int kP = kX.n_cols;  /** Dimension of Beta */
int kK = kZ.n_cols;  /** Dimension of U */

/** sum of yij * (wij - log(1 + ...))
 *  This corresponds to the 
 */
for (int i = 0; i < nObs; i++) {
  double wij = 0;
  for (int j = 0; j < kP; j++) {
    wij += kX(i, j) * beta(j);
  }
  
  for (int j = 0; j < kK; j++) {
    wij += kZ(i, j) * u(j);
  }
  
  value += -exp(wij) + kY(i) * wij;
}

value += ldmn(u, sigma);
return value;
}


// [[Rcpp::export]]
double logAcceptPoisson_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& ucurrent, 
const arma::vec& uproposed, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  return min0(0.0, loglikelihoodPoissonCpp_n(beta, sigma, uproposed, kY, kX, kZ) - loglikelihoodPoissonCpp_n(beta, sigma, ucurrent, kY, kX, kZ));
}

// [[Rcpp::export]]
arma::mat uSamplerPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0) {
  RNGScope scope;
  int kK = u.n_rows;
  
  arma::mat usample(B, kK);
  arma::vec ucurrent(kK);
  arma::vec uproposed(kK);
  ucurrent = u;
  usample.row(0) = ucurrent.t();
  
  for (int i = 1; i < B; i++){
    // uproposed = rnorm(kK, 0, sd0);
    for (int j = 0; j < kK; j++) {
      uproposed(j) = rnorm(1, 0 , sd0 * sqrt(sigma(j, j)))(0);
    }
    uproposed += ucurrent;
    if (log(R::runif(0, 1)) < logAcceptPoisson_n(beta, sigma, ucurrent, uproposed, kY, kX, kZ)) {
      ucurrent = uproposed;
    }
    usample.row(i) = ucurrent.t();
  }
  
  return usample;
}
