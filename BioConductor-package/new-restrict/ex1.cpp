#include <Rcpp.h>

using std::vector;
using std::back_inserter;

// bool match_negative_epist(std::vector<int> E, std::vector<int> G) {
//   // When we have things like -1, 2 in epistasis. We need to check 2 is
//   // present and 1 is not present. E is the vector of epistatic coeffs,
//   // and G the sorted genotype.

//   if(G.size() < 1) return false;
  
//   for(auto e : E) {
//     if(e < 0) {
//       if(binary_search(G.begin(), G.end(), -e))
// 	return false;
//     } else {
//       if(!binary_search(G.begin(), G.end(), e))
// 	return false;
//     }
//   }
//   return true;
// }


// // [[Rcpp::export]]
// void f1() {
//   std::vector<int> e1 = {-2, 1, 3};
//   std::vector<int> E1 = {1, 2, 4};
  
//   Rcpp::Rcout << "A: false. Gives " << match_negative_epist(e1, E1) << "\n";

//   e1 = {-2, 1, 3};
//   E1 = {5, 6, 7};

//   Rcpp::Rcout << "B: false.Gives  " << match_negative_epist(e1, E1) << "\n";

//   e1 = {1, 2, 3};
//   E1 = {1, 2, 3, 4};

//   Rcpp::Rcout << "C: true. Gives " << match_negative_epist(e1, E1) << "\n";

//   e1 = {-2, 1, 3};
//   E1 = {1, 4};

//   Rcpp::Rcout << "D: false. Gives " << match_negative_epist(e1, E1) << "\n"; 

//   e1 = {-3, 1};
//   E1 = {1};

//   Rcpp::Rcout << "E: true. Gives  " << match_negative_epist(e1, E1) << "\n"; 

//   e1 = {-3, 1};
//   E1 = {2};

//   Rcpp::Rcout << "F: false. Gives " << match_negative_epist(e1, E1) << "\n"; 


//   e1 = {-3, 1};
//   E1 = {};

//   Rcpp::Rcout << "G: false. Gives " << match_negative_epist(e1, E1) << "\n"; 

//   e1 = {-2, 1};
//   E1 = {1, 4};

//   Rcpp::Rcout << "H: true. Gives " << match_negative_epist(e1, E1) << "\n"; 

  
// }

// // [[Rcpp::export]]
// double f4(double x) {
//   return 3 * x;
// }

// // [[Rcpp::export]]
// int f3(Rcpp::IntegerVector v1) {
//   std::vector<int> v2 = Rcpp::as<std::vector<int> > (v1);
//   return v2[0];
// }


inline double prodDeathFitness(vector<double> s) {
  double f = 1.0;
  for(auto si : s) {
    if( si <= -90.0 ) {
      return 99.0;
    } else {
      f *= (1 - si);
    }
  }
  return f;
}



// [[Rcpp::export]]
double f1(Rcpp::NumericVector v1) {
  std::vector<double> v2 = Rcpp::as<std::vector<double> > (v1);
  return(prodDeathFitness(v2));

}

// [[Rcpp::export]]
void fi(Rcpp::IntegerVector v1) {
  std::vector<int> v2 = Rcpp::as<std::vector<int> > (v1);
  for(auto i : v2)
    Rcpp::Rcout << " i = " << i;
  Rcpp::Rcout << std::endl;
}


// [[Rcpp::export]]
void fi2(Rcpp::IntegerVector v1) {
  std::vector<int> v2 = Rcpp::as<std::vector<int> > (v1);
  for(int i : v2)
    Rcpp::Rcout << " i = " << i;
  Rcpp::Rcout << std::endl;
}



// [[Rcpp::export]]
void fi3(Rcpp::IntegerVector v1) {
  std::vector<int> v2 = Rcpp::as<std::vector<int> > (v1);
  for(int i : v2)
    Rcpp::Rcout << " i = " << i;
  Rcpp::Rcout << std::endl;
}




// void fs(Rcpp::CharacterVector st) {
//   const std::string typeFitness = Rcpp::as<std::string>(st);
//   Rcpp::Rcout << "You entered " << typeFitness << std::endl;
//   // Rcpp::Rcout << "You entered " << st << std::endl;
// }


