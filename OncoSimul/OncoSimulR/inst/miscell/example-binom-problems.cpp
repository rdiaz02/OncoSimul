#include <Rcpp.h>

using namespace Rcpp ;

// [[Rcpp::export]]
void precissionLoss(NumericVector nin_){
  // We are storing population sizes as doubles.
  // Should not loose any precission up to 2^53 - 1
  // (e.g., http://stackoverflow.com/a/1848762)
  // but double check if optims break it.
  // Note that the original code by Mather stores it as int.

  // Problems are likely to arise soon, with 4.5e15, because
  // of rbinom. See notes in example-binom-problems.cpp
  double a, b, c, d;
  int e, f;
  a = pow(2, 52) + 1.0;	
  b = pow(2, 52); // 2^53 a little over 9*1e15
  c = (9.0 * 1e15) + 1.0;
  d = (9.0 * 1e15);

  e = static_cast<int>(a - b);
  f = static_cast<int>(c - d);

  if( a == b) std::cout << "WARNING!!!! \n Precission loss: a == b\n";
  if( !(a > b)) std::cout << "WARNING!!!! \n Precission loss: !(a > b)\n";
  if(c == d) std::cout << "WARNING!!!! \n Precission loss: c == d\n";
  if( !(c > d)) std::cout << "WARNING!!!! \n Precission loss: !(c > d)\n";
  if( e != 1 ) std::cout << "WARNING!!!! \n Precission loss: e != 1\n";
  if( f != 1 ) std::cout << "WARNING!!!! \n Precission loss: f != 1\n";

  // Now the binomial issue
  double nin, r;
  nin = as<double>(nin_);
  r = floor(nin + 0.5);
  if( r != nin ) {
    std::cout << "\n r != nin \n";
    std::cout << "\n r - nin is " << (r - nin) << "\n";
  }
}
// precissionLoss(1 + 5e15)

// [[Rcpp::export]]
double Bin(NumericVector num, NumericVector p) {
  double n = as<double>(num);
  std::cout << "\n Num = " << n << "\n";
  return ::Rf_rbinom(n, as<double>(p));
}
// [[Rcpp::export]]
double Bin1(NumericVector num, NumericVector p) {
  double n = as<double>(num);
  n = n + 1.0;
  std::cout << "\n Num = " << n << "\n";
  return ::Rf_rbinom(n, as<double>(p));
}
// very strange: Bin can do 9e15, but this fails with 1e15 and above.
// but not if we add large numbers, say 10000. See below
// with Bin2

// [[Rcpp::export]]
double Bin11(NumericVector num, NumericVector p) {
  double n = as<double>(num);
  double oldn = n;
  ++n;
  std::cout << "\n Num = " << n << "\n";
  std::cout << "\n Num - oldn " << (n - oldn) << "\n";

  return ::Rf_rbinom(n, as<double>(p));
}
// [[Rcpp::export]]
double Bin2(NumericVector num, NumericVector a, NumericVector p) {
  double n = num(0); //as<double>(num);
  double aa = a(0);
  double oldn = n;
  n += aa;
  std::cout << "\n n = " << n << "\n";
  std::cout << "\n n - oldn = " << (n - oldn) << "\n";
  double b1 = ::Rf_rbinom(n, as<double>(p));
  double b2 = ::Rf_rbinom(oldn, as<double>(p));

  std::cout << "\n with small number " << b2;
  std::cout << "\n with large number " << b1 <<"\n";

  // The binom code in binom.c
  double nin, r;
  nin = n;
  r = floor(nin + 0.5);
  if( r != nin ) {
    std::cout << "\n r != nin \n";
    std::cout << "\n r - nin is " << (r - nin) << "\n";
  }
  return b1;
}




// Problem is in largest integer storable,
// but then we add a 0.5, ...
// So we can only go up to 2^52.
// run precissionLoss with 4e15, 5e15, 2^52, 2^53


// Notice this: Bin2(4e15, 1 + 1e15, .4)
// But this works Bin2(5e15, 2, .4)
// Though this breaks Bin2(5e15, 1, .4)


// See rbinom.c for the cause of the problem.
// It is     r = floor(nin + 0.5);
//    if (r != nin) ML_ERR_return_NAN;

// This is with odd numbers, not even. E.g.: if x > 5e15 and odd, it fails.
// 6e15 -> x; y <- (x + 1); floor(y + 0.5) == y 

// [[Rcpp::export]]
double NBin(NumericVector num, NumericVector p) {
  return ::Rf_rnbinom(as<double>(num), as<double>(p));
}


// check we do not get overflows as follows
// a <- Bin1(4e15, 1); b <- Bin(4e15, 1); a - b

