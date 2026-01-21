// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::export]]
double dmvnormArma(arma::vec x, arma::vec mean, arma::mat sigma) {
  int n = x.n_elem;
  double logdet = sum(log(eig_sym(sigma)));
  double constant = -0.5 * n * log(2 * M_PI) - 0.5 * logdet;
  arma::vec x_minus_mean = x - mean;
  double exponent = -0.5 * as_scalar(x_minus_mean.t() * inv_sympd(sigma) * x_minus_mean);
  return constant + exponent;
}



// [[Rcpp::export]]
double rtruncnormRcpp(int n, double a, double b, double mean, double sd){
  // Function that uses an R package
  Rcpp::Environment truncnormEnv = Rcpp::Environment::namespace_env("truncnorm");
  Rcpp::Function rtruncnorm = truncnormEnv["rtruncnorm"];   
  
  // Check if b is Inf, and if so, use a very large value
  if (std::isinf(b)) {
    b = std::numeric_limits<double>::max();
  }
  
  
  // Call the function and retrieve the result
  double result = as<double>(rtruncnorm(Rcpp::Named("n") = n,
                                        Rcpp::Named("a") = a,
                                        Rcpp::Named("b") = b,
                                        Rcpp::Named("mean") = mean,
                                        Rcpp::Named("sd") = sd));
  
  return result;
}



// [[Rcpp::export]]
arma::mat varArma(double tau, double rho, arma::mat Dv, arma::mat nbd_mat) {
  // Calculate the variance matrix
  arma::mat var = tau * arma::inv(Dv - rho * nbd_mat);
  
  return var;
}



// [[Rcpp::export]]
arma::mat gibbs_y(
    int M,
    arma::mat y,                // n x m
    List nbd_index,             // length n, 1-based indices
    arma::vec num_nei,           // length n
    arma::vec true_alpha,        // length m
    double true_rho,
    arma::mat G                  // m x m covariance matrix
) {
  
  int n = y.n_rows;
  int m = y.n_cols;
  
  if ((int)true_alpha.n_elem != m)
    stop("true_alpha must have length m (number of columns of y)");
  
  if ((int)num_nei.n_elem != n)
    stop("num_nei must have length n (number of rows of y)");
  
  // ----------------------------------------------------
  // Symmetrize + safe Cholesky of G
  // ----------------------------------------------------
  arma::mat Gsym = 0.5 * (G + G.t());
  
  arma::mat cholG;
  bool ok = arma::chol(cholG, Gsym, "lower");
  
  if (!ok) {
    double eps = 1e-8;
    arma::mat G_jitter = Gsym + eps * arma::eye(m, m);
    ok = arma::chol(cholG, G_jitter, "lower");
    if (!ok)
      stop("G is not positive definite even after jitter");
  }
  
  arma::rowvec alpha_row = true_alpha.t();
  
  // ----------------------------------------------------
  // Gibbs sampler
  // ----------------------------------------------------
  for (int aa = 0; aa < M; aa++) {
    for (int k = 0; k < n; k++) {
      
      IntegerVector nbd = nbd_index[k];
      int nk = nbd.size();
      
      arma::rowvec mu = alpha_row;
      
      if (nk > 0) {
        arma::rowvec neigh_sum(m, arma::fill::zeros);
        for (int j = 0; j < nk; j++) {
          int id = nbd[j] - 1;   // convert to 0-based
          neigh_sum += y.row(id);
        }
        
        // EXACT equivalent of:
        // colSums(y[nbd_id, ] - true_alpha)
        mu += (true_rho / num_nei[k]) *
          (neigh_sum - nk * alpha_row);
      }
      
      // MVN draw: N(mu, G / num_nei[k])
      arma::vec z = arma::randn<arma::vec>(m);
      y.row(k) = mu + (cholG * z / std::sqrt(num_nei[k])).t();
    }
  }
  
  return y;
}




// [[Rcpp::export]]
double log_mvn_prec(const arma::vec& x,
                    const arma::mat& Q,
                    arma::mat& cholQ) {
  
  // cholQ = chol(Q)
  if (!arma::chol(cholQ, Q)) return -arma::datum::inf;
  
  double quad = arma::dot(x, Q * x);
  double logdet = 2.0 * sum(log(cholQ.diag()));
  
  return 0.5 * logdet - 0.5 * quad;
}

// [[Rcpp::export]]
List estX_fun_k(int M,
                     arma::vec X,
                     arma::mat Dv,
                     arma::mat W,
                     arma::vec initial,
                     double S_tau,
                     double S_rho,
                     bool autotune) {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  std::unordered_set<double> uniq_tau;
  std::unordered_set<double> uniq_rho;
  
 
  arma::vec tau(M), rho(M);
  tau(0) = initial(0);
  rho(0) = initial(1);
  
  uniq_tau.insert(tau(0));
  uniq_rho.insert(rho(0));
    
  
  int acc_tau = 0, acc_rho = 0;
  double rate_tau = 0.0;
  double rate_rho = 0.0;
  
  arma::mat Q, cholQ;
  
  for (int m = 1; m < M; m++) {
    
    // ---- TAU update ----
    double tau_prop = rtruncnormRcpp(1, 0, INFINITY, tau(m-1), S_tau);
    
    Q = (Dv - rho(m-1) * W) / tau_prop;
    double lp_prop = log_mvn_prec(X, Q, cholQ);
    
    Q = (Dv - rho(m-1) * W) / tau(m-1);
    double lp_curr = log_mvn_prec(X, Q, cholQ);
    
    if (std::log(R::runif(0,1)) < lp_prop - lp_curr) {
      tau(m) = tau_prop;
      acc_tau++;
    } else {
      tau(m) = tau(m-1);
    }
    
    // ---- RHO update ----
    double rho_prop = rtruncnormRcpp(1, 0, 1, rho(m-1), S_rho);
    
    Q = (Dv - rho_prop * W) / tau(m);
    lp_prop = log_mvn_prec(X, Q, cholQ);
    
    Q = (Dv - rho(m-1) * W) / tau(m);
    lp_curr = log_mvn_prec(X, Q, cholQ);
    
    if (std::log(R::runif(0,1)) < lp_prop - lp_curr) {
      rho(m) = rho_prop;
      acc_rho++;
    } else {
      rho(m) = rho(m-1);
    }
    
    uniq_tau.insert(tau(m));
    uniq_rho.insert(rho(m));
    
    rate_tau = double(uniq_tau.size()) / (m + 1);
    rate_rho = double(uniq_rho.size()) / (m + 1);
    
    // optional tuning every 50 iters
    if (autotune) {
      double rt = rate_tau;
      double rr = rate_rho;
      if (rt > 0.4) S_tau *= 1.1;
      if (rt < 0.2) S_tau /= 1.1;
      if (rr > 0.4) S_rho *= 1.1;
      if (rr < 0.2) S_rho /= 1.1;
    }
  }
  
  
  // Record the end time
  auto end_time = std::chrono::high_resolution_clock::now();
  
  // Calculate the elapsed time in seconds
  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  
  
  return List::create(
    _["tau_samples"] = tau,
    _["rho_samples"] = rho,
    _["rate_tau"] = rate_tau,
    _["rate_rho"] = rate_rho,
    _["S_tau"] = S_tau,
    _["S_rho"] = S_rho,
    _["elapsed_time_sec"] = elapsed_seconds.count()
  );
}



