#include <Rcpp.h>
using namespace Rcpp;


std::vector < std::vector<int> > list_to_vector_of_int_vectors(Rcpp::List vlist) {
  std::vector < std::vector<int> > vv(vlist.size());
  for(int i = 0; i != vlist.size(); ++i) {
    vv[i] = Rcpp::as<std::vector<int> >(vlist[i]);
  }
  return vv;
}

// [[Rcpp::export]]
void wrap_list_to_vector_of_int_vectors(Rcpp::List vlist) {
  std::vector < std::vector<int> > vo(vlist.size());
  vo = list_to_vector_of_int_vectors(vlist);
  for(int ii = 0; ii != vo.size(); ++ii) {
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << " list position " << ii + 1 << ": ";
    for(int jj = 0; jj != vo[ii].size(); ++jj ) {
      Rcpp::Rcout << vo[ii][jj] << " ";
    }
  }
  Rcpp::Rcout << "\n";
}





// x <- c(1, 12.7, 3.5); n <- 1000000; uu <- replicate(n, f1(x))
// x/sum(x)
