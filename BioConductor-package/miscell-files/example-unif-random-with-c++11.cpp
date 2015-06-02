// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <random>

using namespace Rcpp;

// [[Rcpp::export]]
void f1(IntegerVector num) {
  Rcpp::Rcout << "hola" << "\n";
  int n = as<int>(num);
  std::random_device random_engine;
  std::mt19937 gen(random_engine());
  std::uniform_int_distribution<int> unif(0, n);

  for(int i =0; i < 10; ++i) 
    Rcpp::Rcout << unif(gen) << ' ';
  Rcpp::Rcout << "\n";
}


// [[Rcpp::export]]
void f2(IntegerVector num, IntegerVector s) {
  Rcpp::Rcout << "hola" << "\n";
  int n = as<int>(num);
  int seed = as<int>(s);

  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> unif(0, n);

  for(int i =0; i < 10; ++i) 
    Rcpp::Rcout << unif(gen) << ' ';
  Rcpp::Rcout << "\n";
}
