#include <Rcpp.h>
#include <iterator>
#include <algorithm>

bool check_order(std::vector<int> O, std::vector<int> G) {

  if(G.size() < O.size()) return false;
  
  std::vector<int>::iterator p;
  std::vector<size_t> vdist;
  
  auto itb = G.begin();
  
  for(auto o : O) {
    p = find(G.begin(), G.end(), o);
    if( p == G.end() ) {
      return false;
    } else {
      vdist.push_back(std::distance( itb, p ));
    }
  }
  // Rcpp::Rcout << "     ";
  // for(auto vv : vdist ) {
  //   Rcpp::Rcout << vv << " ";
  // }
  // Rcpp::Rcout << std:: endl;
  if( is_sorted(vdist.begin(), vdist.end()) ) {
    return true;
  } else {
    return false;
  }
}


// [[Rcpp::export]]
void f1(){
  std::vector<int> O1 = {2, 1, 3};
  std::vector<int> G1 = {1, 2, 3, 4, 5};
  bool res;

  res = check_order(O1, G1);
  Rcpp::Rcout << " A. It is false: " << res << std::endl;

  O1 = {1, 2};
  G1 = {};
  res = check_order(O1, G1);
  Rcpp::Rcout << " B. It is false: " << res << std::endl;

  O1 = {1, 2, 3};
  G1 = {1, 2};
  res = check_order(O1, G1);
  Rcpp::Rcout << " C. It is false: " << res << std::endl;

  O1 = {1, 2, 4};
  G1 = {1, 2, 6, 8, 9, 7, 3, 11, 21, 4};
  res = check_order(O1, G1);
  Rcpp::Rcout << " D. It is true: " << res << std::endl;


  O1 = {1, 4};
  G1 = {1, 2, 6, 8, 9, 7, 3, 11, 21, 4};
  res = check_order(O1, G1);
  Rcpp::Rcout << " E. It is true: " << res << std::endl;

  O1 = {1, 3, 2, 4};
  G1 = {1, 2, 6, 8, 9, 7, 3, 11, 21, 4};
  res = check_order(O1, G1);
  Rcpp::Rcout << " F. It is false: " << res << std::endl;

  O1 = {3, 2};
  G1 = {1, 2};
  res = check_order(O1, G1);
  Rcpp::Rcout << " G. It is false: " << res << std::endl;

  O1 = {3, 2};
  G1 = {2, 3};
  res = check_order(O1, G1);
  Rcpp::Rcout << " H. It is false: " << res << std::endl;

  O1 = {3, 2};
  G1 = {3, 2};
  res = check_order(O1, G1);
  Rcpp::Rcout << " I. It is true: " << res << std::endl;

}
