#include <Rcpp.h>

using namespace Rcpp ;


struct Deps {
  std::vector<int> module;
  double s;
  double sh;
};

std::multimap<int, Deps> restrictTable;


// [[Rcpp::export]]
void f1() {
  //int n = 5;
  std::vector<int> v1;
  for(int i = 0; i != 5; ++i) {
    for(int j = 0; j != (i + 1); ++j)
      v1.push_back(j);
    Deps d1;
    d1.module = v1;
    d1.s = 0.1 + i;
    d1.sh = 0.6 + i;

    restrictTable.insert(std::make_pair(i, d1));
    v1.clear();
    
  }
  std::cout << "\n The multimap contains \n";
  std::multimap<int, Deps>::iterator it;
  for(it = restrictTable.begin(); it != restrictTable.end(); ++it ) {
    std::cout << it->first << " => " << " Module = ";
    for(auto vv : it->second.module) std::cout << vv << ' ';
    std::cout << " s = " << it->second.s << "; sh = " << it->second.sh << std::endl;
  }
  restrictTable.clear();
}

