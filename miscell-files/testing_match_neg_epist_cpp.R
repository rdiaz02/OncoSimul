#include <Rcpp.h>

// struct A {
//   std::vector<int> v;
//   int a;
// };


bool match_negative_epist(std::vector<int> E, std::vector<int> G) {
  // When we have things like -1, 2 in epistasis. We need to check 2 is
  // present and 1 is not present. E is the vector of epistatic coeffs,
  // and G the sorted genotype.

  if(G.size() < 1) return false;
  
  for(auto e : E) {
    if(e < 0) {
      if(binary_search(G.begin(), G.end(), -e))
	return false;
    } else {
      if(!binary_search(G.begin(), G.end(), e))
	return false;
    }
  }
  return true;
}


// [[Rcpp::export]]
void f1() {
  std::vector<int> e1 = {-2, 1, 3};
  std::vector<int> E1 = {1, 2, 4};
  
  Rcpp::Rcout << "A: false. Gives " << match_negative_epist(e1, E1) << "\n";

  e1 = {-2, 1, 3};
  E1 = {5, 6, 7};

  Rcpp::Rcout << "B: false.Gives  " << match_negative_epist(e1, E1) << "\n";

  e1 = {1, 2, 3};
  E1 = {1, 2, 3, 4};

  Rcpp::Rcout << "C: true. Gives " << match_negative_epist(e1, E1) << "\n";

  e1 = {-2, 1, 3};
  E1 = {1, 4};

  Rcpp::Rcout << "D: false. Gives " << match_negative_epist(e1, E1) << "\n"; 

  e1 = {-3, 1};
  E1 = {1};

  Rcpp::Rcout << "E: true. Gives  " << match_negative_epist(e1, E1) << "\n"; 

  e1 = {-3, 1};
  E1 = {2};

  Rcpp::Rcout << "F: false. Gives " << match_negative_epist(e1, E1) << "\n"; 


  e1 = {-3, 1};
  E1 = {};

  Rcpp::Rcout << "G: false. Gives " << match_negative_epist(e1, E1) << "\n"; 

  e1 = {-2, 1};
  E1 = {1, 4};

  Rcpp::Rcout << "H: true. Gives " << match_negative_epist(e1, E1) << "\n"; 

  
}


// // [[Rcpp::export]]
// void f1(){
//   std::vector<int> A = {1, 2, 3, 4};
//   std::vector<int> B;
//   std::vector<int> D = {3};

//   int C = 3;


//   if( binary_search(D.begin(), D.end(), C))
//     Rcpp::Rcout << "In D\n";

//   if( binary_search(B.begin(),  B.end(), C))
//     Rcpp::Rcout << "In B\n";

  

  
//   if( binary_search(A.begin(), A.end(), C))
//     Rcpp::Rcout << "In A\n";
  
  
//   // A A1;

//   // A1.v = {1, 2, 3, 4};
//   // A1.a = 78;


//   // A A2 = A1;
//   // A2.v.push_back(96);


//   // A A3 = A2;
//   // A3.a = 135;

//   // A2.v.push_back(969);
//   // A2.v.push_back(96999);
  


//   // Rcpp::Rcout << "\n A1 ";
//   // for(auto it = A1.v.begin(); it != A1.v.end(); ++it) {
//   //   Rcpp::Rcout << (*it) << ", ";
//   // }
//   // Rcpp::Rcout << "\n A.a " << A1.a << std::endl;

//   // Rcpp::Rcout << "\n A2 ";
//   // for(auto it = A2.v.begin(); it != A2.v.end(); ++it) {
//   //   Rcpp::Rcout << (*it) << ", ";
//   // }
//   // Rcpp::Rcout << "\n A2.a " << A2.a << std::endl;

//   // Rcpp::Rcout << "\n A3 ";
//   // for(auto it = A3.v.begin(); it != A3.v.end(); ++it) {
//   //   Rcpp::Rcout << (*it) << ", ";
//   // }
//   // Rcpp::Rcout << "\n A3.a " << A3.a << std::endl;
// }
  
