#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void list_to_vector_of_int_vectors(Rcpp::List vlist) {
  std::vector < std::vector<int> > vv(vlist.size());

  for(int i = 0; i != vlist.size(); ++i) {
    vv[i] = Rcpp::as<std::vector<int> >(vlist[i]);
  }

  // Check 
  for(int ii = 0; ii != vv.size(); ++ii) {
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << " list position " << ii + 1 << ": ";
    for(int jj = 0; jj != vv[ii].size(); ++jj ) {
      Rcpp::Rcout << vv[ii][jj] << " ";
    }
  }
  Rcpp::Rcout << "\n";
  // we will later return the vv
}





// x <- c(1, 12.7, 3.5); n <- 1000000; uu <- replicate(n, f1(x))
// x/sum(x)
