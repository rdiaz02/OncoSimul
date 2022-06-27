// Code from file shared.h in the R package extraDistr
// https://github.com/twolodzko/extraDistr
// https://github.com/twolodzko/extraDistr/blob/master/src/shared.h
// Author: Tymoteusz Wolodzko
// License: GPL-2
// Downloaded on 2020-03-18


#include <Rcpp.h>
using namespace Rcpp ;
using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

static const double MIN_DIFF_EPS = 1e-8;

// MACROS

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))

Rcpp::NumericMatrix my_rmvhyper(const int nn, const Rcpp::NumericMatrix n, const NumericVector k);
