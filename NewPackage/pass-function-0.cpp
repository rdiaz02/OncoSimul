#include <Rcpp.h>

using namespace Rcpp;

// Simplified from Dirk's book, p. 56

// [[Rcpp::export]]
NumericVector f1(Function g, NumericVector y) {
  return g(y);
}

// This is more verbose, and we do not needed the Function
// NumericVector f1(Function g, NumericVector y) {
//   Function h(g);
//   return h(y);
// }


// R Code
// s1 <- function(x) sort(x, decreasing = TRUE)
// f1(s1, 1:7)
