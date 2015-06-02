#include <Rcpp.h>

using namespace Rcpp;


struct st {
  double x1;
};

// // [[Rcpp::export]]
// NumericVector f2(NumericVector y) {
//   st st1;
//   st1.x1 = 3*as<double>(y);
//   Rcpp::Rcout << " Value of x1 is " << st1.x1 << std::endl;
//   return wrap(st1.x1);
// }




// [[Rcpp::export]]
NumericVector f3(Function g3, NumericVector y) {
  st st1;
  st1.x1 = as<double>(g3(as<double>(y)));
  Rcpp::Rcout << " Value of x1 is " << st1.x1 << std::endl;
  return wrap(st1.x1);
}

// g31 <- function(x) 4*x
// f3(g31, 5)



// [[Rcpp::export]]
NumericVector f4(Function g4, NumericVector y) {
  st st1;
  g4(as<double>(y));
  Rcpp::Rcout << " Value of x1 is " << st1.x1 << std::endl;
  return wrap(st1.x1);
}

