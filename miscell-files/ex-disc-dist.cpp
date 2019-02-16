#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int f1(Rcpp::NumericVector weights_) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::vector<double> weights = Rcpp::as<std::vector<double> >(weights_);
  std::discrete_distribution<int> dt(weights.begin(), weights.end());

  // If you want to see what it thinks.
  // Rcpp::Rcout << "\n Probabilities ";
  // for( double x:dt.probabilities() )
  //   Rcpp::Rcout << x << " ";

  return dt(gen);
}



// x <- c(1, 12.7, 3.5); n <- 1000000; uu <- replicate(n, f1(x))
// x/sum(x)
