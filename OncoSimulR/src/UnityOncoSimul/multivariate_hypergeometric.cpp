// Code taken from several files from R package extraDistr
// (shared.cpp, shared_inline.h, multivariate_hypergeometric.cpp)
// https://github.com/twolodzko/extraDistr
// Author: Tymoteusz Wolodzko
// License: GPL-2
// Downloaded on 2020-03-18

#include "multivariate_hypergeometric.h"

double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}



// inline functions


inline bool tol_equal(double x, double y) {
  return std::abs(x - y) < MIN_DIFF_EPS;
}

inline double phi(double x) {
  return R::dnorm(x, 0.0, 1.0, false);
}

inline double lphi(double x) {
  return R::dnorm(x, 0.0, 1.0, true);
}

inline double Phi(double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}

inline double InvPhi(double x) {
  return R::qnorm(x, 0.0, 1.0, true, false);
}

inline double factorial(double x) {
  return R::gammafn(x + 1.0);
}

inline double lfactorial(double x) {
  return R::lgammafn(x + 1.0);
}

inline double rng_sign() {
  double u = rng_unif();
  return (u > 0.5) ? 1.0 : -1.0;
}

inline bool is_large_int(double x) {
  if (x > std::numeric_limits<int>::max())
    return true;
  return false;
}

inline double to_dbl(int x) {
  return static_cast<double>(x);
}

inline int to_pos_int(double x) {
  if (x < 0.0 || ISNAN(x))
    Rcpp::stop("value cannot be coerced to integer");
  if (is_large_int(x))
    Rcpp::stop("value out of integer range");
  return static_cast<int>(x);
}

inline double trunc_p(double x) {
  return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); 
}



// functions
bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

double finite_max_int(const Rcpp::NumericVector& x) {
  double max_x = 0.0;
  int n = x.length();
  int i = 0;
  do {
    if (x[i] > 0.0 && !is_large_int(x[i])) {
      max_x = x[i];
      break;
    }
    i++;
  } while (i < n);
  while (i < n) {
    if (x[i] > max_x && !is_large_int(x[i])) {
      max_x = x[i];
    }
    i++;
  }
  return max_x;
}

// Renaming of  cpp_rmvhyper in
// https://github.com/twolodzko/extraDistr/blob/master/src/multivariate-hypergeometric-distribution.cpp
Rcpp::NumericMatrix my_rmvhyper(
    const int nn,
    const Rcpp::NumericMatrix n,
    const NumericVector k
  ) {
  
  if (std::min({static_cast<int>(n.nrow()),
                static_cast<int>(n.ncol()),
                static_cast<int>(k.length())}) < 1) {
    Rcpp::warning("NAs produced");
    Rcpp::NumericMatrix out(nn, n.ncol());
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
  int m = n.ncol();
  Rcpp::NumericMatrix x(nn, m);
  std::vector<double> n_otr(m);
  
  bool wrong_values;
  double k_left;
  
  bool throw_warning = false;

  for (int i = 0; i < nn; i++) {
    
    wrong_values = false;
    n_otr[0] = 0.0;
    
    for (int j = 1; j < m; j++) {
      if (!isInteger(GETM(n, i, j), false) ||
          GETM(n, i, j) < 0.0 || ISNAN(GETM(n, i, j))) {
        wrong_values = true;
        break;
      }
      n_otr[0] += GETM(n, i, j);
    }
    
    if (wrong_values || ISNAN(GETV(k, i)) || ISNAN(GETM(n, i, 0)) ||
        !isInteger(GETM(n, i, 0), false) || GETM(n, i, 0) < 0 ||
        (n_otr[0] + GETM(n, i, 0)) < GETV(k, i) ||
        !isInteger(GETV(k, i), false) || GETV(k, i) < 0.0) {
      throw_warning = true;
      for (int j = 0; j < m; j++)
        x(i, j) = NA_REAL;
      continue;
    }
    
    for (int j = 1; j < m; j++)
      n_otr[j] = n_otr[j-1] - GETM(n, i, j);
    
    k_left = GETV(k, i);
    x(i, 0) = R::rhyper(GETM(n, i, 0), n_otr[0], k_left);
    k_left -= x(i, 0);
    
    if (m > 2) {
      for (int j = 1; j < m-1; j++) {
        x(i, j) = R::rhyper(GETM(n, i, j), n_otr[j], k_left);
        k_left -= x(i, j);
      }
    }
    
    x(i, m-1) = k_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

// [[Rcpp::export]]
/*Rcpp::NumericMatrix wrap_rmvhyperconst (int nn, const Rcpp::NumericMatrix n, const NumericVector k){
  return my_rmvhyper(nn, n, k);
}*/
