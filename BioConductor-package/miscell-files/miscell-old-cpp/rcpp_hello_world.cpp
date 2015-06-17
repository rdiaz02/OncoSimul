// I will start playing with Algo2, and get it to work.


#include "rcpp_hello_world.h"

SEXP rcpp_hello_world(){
    using namespace Rcpp ;
    
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    
    return z ;
}






// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
#include <limits>
                     #include <iostream>

// declarations
extern "C" {
SEXP file6889607c2f1( SEXP R, SEXP W, SEXP death, SEXP growth, SEXP n, SEXP ti) ;
}

// definition

SEXP file6889607c2f1( SEXP R, SEXP W, SEXP death, SEXP growth, SEXP n, SEXP ti ){
BEGIN_RCPP


double xdeath = as<double>(death);
double xgrowth = as<double>(growth);
double xn = as<double>(n);
double xR = as<double>(R);
double xW = as<double>(W);
double xti = as<double>(ti);
double eq11;
double r;
double rr;

// W < 0 is a signal that mutation is zero, and thus ti is Inf
if(xW <= -99.0) {
 xti = std::numeric_limits<double>::infinity();
 } else {

// using version with log
RNGScope scope;
r = ::Rf_runif(0.0, 1.0);
rr = exp((1/xn) * log(r));

std::cout << "r = " << r << " rr = " << rr << std::endl;

// eq.12
if( ((xR - xW + 2 * xdeath)/(xR + xW - 2 * xgrowth)) < rr ) {


std::cout << "numerator = " << (xR - xW + 2 * xgrowth) - xW - xR + 2 * xdeath << std::endl;
std::cout << "denominator = " << (-xR -xW + 2 * xgrowth) - xW + xR + 2 * xdeath << std::endl;


   //eq. 11
   xti = (1/xR) * log( (rr * (xR - xW + 2 * xgrowth) - xW - xR + 2 * xdeath) /
                     (rr * (-xR -xW + 2 * xgrowth) - xW + xR + 2 * xdeath));
   if(xti < 0.0) {
          throw std::range_error("ti: eq.11 < 0");
          }
   // return(wrap(xti));
   } else {
    xti = std::numeric_limits<double>::infinity();
 //return(wrap(xti));
  }
}
return(wrap(xti));

END_RCPP
}
