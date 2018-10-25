#ifndef _DEBUG_COMMON_H__
#define _DEBUG_COMMON_H__

#include<Rcpp.h>

//#define DEBUGZ
// #define DEBUGV
#define DEBUGW

#define DP1(x) {Rcpp::Rcout << "\n DEBUG2: I am at " << x << std::endl;}
#define DP2(x) {Rcpp::Rcout << "\n DEBUG2: Value of " << #x << " = " << x << std::endl;}
#define DP3(x, t){                                                    \
    Rcpp::Rcout <<"\n DEBUG2:" ;                                      \
    for(int xut = 0; xut < t; ++xut) Rcpp::Rcout << "\t ";            \
    Rcpp::Rcout << "  I am at " << x << std::endl;}
#define DP4(x, t){                                                    \
    Rcpp::Rcout <<"\n DEBUG2:" ;                                      \
    for(int xut = 0; xut < t; ++xut) Rcpp::Rcout << "\t ";            \
    Rcpp::Rcout << "  Value of " << #x << " = " << x << std::endl; }

/* void here(std::string x) { */
/*   Rcpp::Rcout << "\n DEBUG: HERE at " << x << std::endl; */
/* } */


// Windows compiler in BioC is pre 4.8.0, so no to_string
// From http://stackoverflow.com/a/5590404
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
       ( std::ostringstream() << std::dec << x ) ).str()




#ifdef DEBUGW
#define ASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << "\n\n";		\
    }									\
  }
#else
#define ASSERT(x);
#endif


#ifdef DEBUGW
#define STOPASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << std::endl;		\
	throw std::out_of_range("STOPASSERT");				\
    }									\
  }
#else
#define STOPASSERT(x);
#endif


#endif


