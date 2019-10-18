#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

float Hbeta(arma::mat& D, float beta, arma::vec& P, int idx) {
  P = arma::exp(-D.col(idx) * beta);
  float sumP = sum(P);
  float H;
  if (sumP == 0){
    H = 0;
    P = D.col(idx) * 0;
  } else {
    H = log(sumP) + beta * sum(D.col(idx) % P) / sumP;
    P /= sumP;
  }
  return(H);
}

// [[Rcpp::export]]
arma::vec compute_simpson_index(arma::mat& D, arma::umat& knn_idx, arma::vec& batch_labels, int n_batches,
                                float perplexity = 15, float tol = 1e-5) {
  int n = D.n_cols;
  arma::vec P = arma::zeros<arma::vec>(D.n_rows);
  arma::vec simpson = arma::zeros<arma::vec>(n);
  float logU = log(perplexity);
  
  float hbeta, beta, betamin, betamax, H, Hdiff;
  int tries;
  for (int i = 0; i < n ; i++) {
    beta = 1;
    betamin = -arma::datum::inf;
    betamax = arma::datum::inf;
    H = Hbeta(D, beta, P, i);
    Hdiff = H - logU;
    tries = 0;
    // first get neighbor probabilities
    while(std::abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0){
        betamin = beta;
        if (!arma::is_finite(betamax)) beta *= 2;
        else beta = (beta + betamax) / 2;
      } else{
        betamax = beta;
        if (!arma::is_finite(betamin)) beta /= 2;
        else beta = (beta + betamin) / 2;
      }
      
      H = Hbeta(D, beta, P, i);
      Hdiff = H - logU;
      tries++;
    }
    
    if (H == 0) {
      simpson.row(i) = -1;
      continue;
    }
    
    // then compute Simpson's Index
    for (int b = 0; b < n_batches; b++) {
      arma::uvec q = find(batch_labels.elem(knn_idx.col(i)) == b); // indices of cells belonging to batch (b)
      if (q.n_elem > 0) {
        float sumP = sum(P.elem(q));
        simpson.row(i) += sumP * sumP;         
      }
    }
  }
  return(simpson);
}


