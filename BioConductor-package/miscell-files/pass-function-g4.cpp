#include <Rcpp.h>

using namespace Rcpp;

struct st {
  double x1;
};

void g4_cpp(st& st1, const double z){
  st1.x1 = 4 * z;
}


typedef void (*funcPtr)(st& st1, const double z);


// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
    if (fstr == "g4")
        return(XPtr<funcPtr>(new funcPtr(&g4_cpp)));
    else
        return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
void callg4(const double zz, std::string funname) {
  st st1;
  XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
  funcPtr fun = *xpfun;
  fun(st1, zz);
  Rcpp::Rcout << " Value of x1 is " << st1.x1 << std::endl;
}


// http://gallery.rcpp.org/articles/passing-cpp-function-pointers/
