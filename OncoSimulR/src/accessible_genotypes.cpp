//     Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
// using namespace Rcpp;

inline int HammingDistance(const Rcpp::IntegerVector& x, const Rcpp::IntegerVector& y) {
  Rcpp::NumericVector diff = Rcpp::abs( x - y );
  return std::accumulate(diff.begin(), diff.end(), 0);
}


// eventually remove this. Left now for testing
// [[Rcpp::export]]
Rcpp::IntegerVector accessibleGenotypes_former(Rcpp::IntegerMatrix y,
					       Rcpp::NumericVector f,
					       Rcpp::IntegerVector numMut, //
					       double th) {
  // Return just the indices. Could preserve fitness, but would need
  // another matrix.
  int ng = y.nrow(); //it counts the wt
  Rcpp::IntegerMatrix adm(ng, ng);
  int numMutdiff = 0;
  // I would have thought this would be faster. It ain't.
  // The last genotype never accesses anything.
  // for(int i = 0; i < (ng - 1); ++i) {
  //   // Candidate genotypes to be accessed from i are always of larger
  //   // mutation by 1. And candidates can thus not have smaller index
  //   for(int j = (i + 1); j < ng; ++j) {
  //     if( (numMut(j) == (numMut(i) + 1)) &&
  // 	  ( (f(j) - f(i)) >= th) &&
  // 	  (HammingDistance(y(i, _), y(j, _)) == 1) ) {
  // 	adm(i, j) = 1;
  //     } else if( (numMut(j) > (numMut(i) + 1)) ) {
  // 	break;
  //     }
  //   }
  // }

  // The last genotype never accesses anything.
  for(int i = 0; i < (ng - 1); ++i) {
    // Candidate genotypes to be accessed from i are always of larger
    // mutation by 1. And candidates can thus not have smaller index
    for(int j = (i + 1); j < ng; ++j) {
      numMutdiff = numMut(j) - numMut(i);
      if( numMutdiff > 1) { // no more to search
  	break; 
      } else if(numMutdiff == 1) {
  	// f(j) - f(i) is faster than HammingDistance
  	// but might lead to more evals?
  	// or fewer, depending on landscape
  	if( ( (f(j) - f(i)) >= th) &&
	    (HammingDistance(y(i, Rcpp::_), y(j, Rcpp::_)) == 1)
  	    ) {
  	  adm(i, j) = 1;
	  // Rcpp::Rcout << "i = " << i << " j = " << j << " adm " << adm(i,j) << "\n"; 
  	}
      }
    }
  }



  // Slightly different logic from R: Do not resize object; set the row to
  // 0.
  int colsum = 0;
  // int indicator = 0; // indicator != 0 means we set one row to 0
  // so we need to iterate at least once more.
  
  // accessible is the genotype number, not the column!  WT is 1,
  // etc. This makes it easy to keep track of which are accessible.
  Rcpp::IntegerVector accessible = Rcpp::seq_len(ng);

  // This is doable in one pass
  // while (true) {
  //   indicator = 0;
    for(int k = 1; k < ng; ++k) {
      if(accessible(k) > 0) {
	colsum = std::accumulate(adm(Rcpp::_, k).begin(),
				 adm(Rcpp::_, k).end(), 0);
	if(colsum == 0) { // This genotype ain't reachable
	  // Nothing can be reached from this genotype; fill with 0.
	  adm(k, Rcpp::_) = Rcpp::IntegerVector(ng);
	  accessible(k) = -9;
	  // indicator = 1;
	}
      }
    }
  //   if(indicator == 0) break;
  // }
  return accessible;
}






// [[Rcpp::export]]
Rcpp::NumericMatrix genot2AdjMat(Rcpp::IntegerMatrix y,
				 Rcpp::NumericVector f,
				 Rcpp::IntegerVector numMut) {
  // Return just the indices. Could preserve fitness, but would need
  // another matrix.
  int ng = y.nrow(); //it counts the wt
  Rcpp::NumericMatrix adm(ng, ng);
  
  // fill with NAs: https://stackoverflow.com/a/23753626
  // Filling with NAs and in general having NAs might lead to performance
  // penalties. But I use the NAs in a lot of the code for accessible
  // genotypes, etc.
  std::fill( adm.begin(), adm.end(), Rcpp::NumericVector::get_na() ) ;
  int numMutdiff = 0;
  // I would have thought this would be faster. It ain't.
  // The last genotype never accesses anything.
  // for(int i = 0; i < (ng - 1); ++i) {
  //   // Candidate genotypes to be accessed from i are always of larger
  //   // mutation by 1. And candidates can thus not have smaller index
  //   for(int j = (i + 1); j < ng; ++j) {
  //     if( (numMut(j) == (numMut(i) + 1)) &&
  // 	  ( (f(j) - f(i)) >= th) &&
  // 	  (HammingDistance(y(i, _), y(j, _)) == 1) ) {
  // 	adm(i, j) = 1;
  //     } else if( (numMut(j) > (numMut(i) + 1)) ) {
  // 	break;
  //     }
  //   }
  // }

  // The last genotype never accesses anything.
  for(int i = 0; i < (ng - 1); ++i) {
    // Candidate genotypes to be accessed from i are always of larger
    // mutation by 1. And candidates can thus not have smaller index
    for(int j = (i + 1); j < ng; ++j) {
      numMutdiff = numMut(j) - numMut(i);
      if( numMutdiff > 1) { // no more to search
  	break; 
      } else if(numMutdiff == 1) {
  	if( HammingDistance(y(i, Rcpp::_), y(j, Rcpp::_)) == 1) {
  	  adm(i, j) =  (f(j) - f(i));
  	}
      }
    }
  }
  return adm;
}


Rcpp::IntegerMatrix integerAdjMat(Rcpp::IntegerMatrix y,
				  Rcpp::NumericVector f,
				  Rcpp::IntegerVector numMut, //
				  double th) {
  // Return a genotype adjacency matrix with a 1 if genotype j is
  // accessible (fitness >, within th) from i.
  int ng = y.nrow(); //it counts the wt
  Rcpp::IntegerMatrix adm(ng, ng);
  int numMutdiff = 0;
  // I would have thought this would be faster. It ain't.
  // The last genotype never accesses anything.
  // for(int i = 0; i < (ng - 1); ++i) {
  //   // Candidate genotypes to be accessed from i are always of larger
  //   // mutation by 1. And candidates can thus not have smaller index
  //   for(int j = (i + 1); j < ng; ++j) {
  //     if( (numMut(j) == (numMut(i) + 1)) &&
  // 	  ( (f(j) - f(i)) >= th) &&
  // 	  (HammingDistance(y(i, _), y(j, _)) == 1) ) {
  // 	adm(i, j) = 1;
  //     } else if( (numMut(j) > (numMut(i) + 1)) ) {
  // 	break;
  //     }
  //   }
  // }

  // The last genotype never accesses anything.
  for(int i = 0; i < (ng - 1); ++i) {
    // Candidate genotypes to be accessed from i are always of larger
    // mutation by 1. And candidates can thus not have smaller index
    for(int j = (i + 1); j < ng; ++j) {
      numMutdiff = numMut(j) - numMut(i);
      if( numMutdiff > 1) { // no more to search
  	break; 
      } else if(numMutdiff == 1) {
  	// f(j) - f(i) is faster than HammingDistance
  	// but might lead to more evals?
  	// or fewer, depending on landscape
  	if( ( (f(j) - f(i)) >= th) &&
	    (HammingDistance(y(i, Rcpp::_), y(j, Rcpp::_)) == 1)
  	    ) {
  	  adm(i, j) = 1;
	  // Rcpp::Rcout << "i = " << i << " j = " << j << " adm " << adm(i,j) << "\n"; 
  	}
      }
    }
  }
  return adm;
}


// used in both peaks and accessible genotypes
Rcpp::IntegerVector accessibleGenotypesPeaksLandscape(Rcpp::IntegerMatrix y,
						      Rcpp::NumericVector f,
						      Rcpp::IntegerVector numMut, //
						      double th,
						      bool returnpeaks) {
  // Return the indices. This is like accessibleGenotypes, but we do an
  // extra loop
  int ng = y.nrow(); //it counts the wt
  Rcpp::IntegerMatrix adm(ng, ng);

  adm = integerAdjMat(y, f, numMut, th);
  
  // int numMutdiff = 0;

  // Slightly different logic from R: Do not resize object; set the row to
  // 0.
  int colsum = 0;
  // int indicator = 0; // indicator != 0 means we set one row to 0
  // so we need to iterate at least once more.
  
  // accessible is the genotype number, not the column!  WT is 1,
  // etc. This makes it easy to keep track of which are accessible.
  Rcpp::IntegerVector accessible = Rcpp::seq_len(ng);
  // This is doable in one pass
  // while (true) {
  //   indicator = 0;
    for(int k = 1; k < ng; ++k) {
      if(accessible(k) > 0) {
	colsum = std::accumulate(adm(Rcpp::_, k).begin(),
				 adm(Rcpp::_, k).end(), 0);
	if(colsum == 0) { // This genotype ain't reachable
	  // Nothing can be reached from this genotype; fill with 0.
	  adm(k, Rcpp::_) = Rcpp::IntegerVector(ng);
	  accessible(k) = -9;
	  // indicator = 1;
	}
      }
    }
  //   if(indicator == 0) break;
  // }
    if(!returnpeaks) {
      return accessible;
    } else  {
      // BEWARE: this will not work if several connected genotypes
      // have the same fitness and are maxima
      int rowsum = 0;
      Rcpp::IntegerVector peaks;
      for(int k = 0; k < ng; ++k) {
	if(accessible(k) > 0) {
	  rowsum = std::accumulate(adm(k, Rcpp::_).begin(),
				   adm(k, Rcpp::_).end(), 0);
	  if(rowsum == 0) { // This genotype doesn't have children
	    peaks.push_back(k + 1); // k is index. But in R, WT is in pos 1
	  }
	}
      }
      return peaks;
    }
}



// [[Rcpp::export]]
Rcpp::IntegerVector accessibleGenotypes(Rcpp::IntegerMatrix y,
					Rcpp::NumericVector f,
					Rcpp::IntegerVector numMut, //
					double th) {
  return accessibleGenotypesPeaksLandscape(y, f, numMut, th, false);
}

// [[Rcpp::export]]
Rcpp::IntegerVector peaksLandscape(Rcpp::IntegerMatrix y,
				   Rcpp::NumericVector f,
				   Rcpp::IntegerVector numMut, //
				   double th) {
  return accessibleGenotypesPeaksLandscape(y, f, numMut, th, true);
}



// // This would make it easier returning the actual accessible genotypes easily
// // preserving the fitness if needed
// // Not being used now
// // [[Rcpp::export]]
// IntegerVector acc_ge(Rcpp::IntegerMatrix y, Rcpp::NumericVector f,
// 		     Rcpp::IntegerVector numMut,
// 		     int ng, //it counts the wt
// 		     double th) {
  
//   IntegerMatrix adm(ng, ng);
//   int numMutdiff = 0;
  
//   for(int i = 0; i < (ng - 1); ++i) {
//     // Candidates are always of larger mutation by 1
//     for(int j = (i + 1); j < ng; ++j) {
//       numMutdiff = numMut(j) - numMut(i);
//       if(numMutdiff > 1) { // no more to search
// 	break; 
//       } else if(numMutdiff == 1) {
// 	if( ( (f(j) - f(i)) >= th) &&
// 	    (HammingDistance(y(i, _), y(j, _)) == 1) ) {
// 	  adm(i, j) = 1;
// 	}
//       }
//     }
//   }
//   // Keeps root in Rows
//   IntegerMatrix admtmp = adm(Range(0, ng - 1), Range(1, ng - 1));

//   // Slightly different logic from R: Do not resize object; set the row to
//   // 0.
//   int colsum = 0;
//   int indicator = 0; // indicator != 0 means we set one row to 0
//   // so we need to iterate at least once more.
  
//   // accessible is the genotype number, not the column!  WT is 1,
//   // etc. This makes it easy to keep track of which are accessible.
//   IntegerVector accessible = seq_len(ng - 1) + 1;
  
//   while (true) {
//     indicator = 0;
//     for(int k = 0; k < (ng - 1); ++k) {
//       if(accessible(k) > 0) {
// 	colsum = std::accumulate(admtmp(_, k).begin(),
// 				 admtmp(_, k).end(), 0);
// 	if(colsum == 0) { // This genotype ain't reachable
// 	  // Recall row keeps Root.
// 	  // Nothing can be reached from this genotype; fill with 0.
// 	  admtmp(k + 1, _) = IntegerVector(ng - 1);
// 	  accessible(k) = -9;
// 	  indicator = 1;
// 	}
//       }
//     }
//     if(indicator == 0) break;
//   }
//   return accessible;
// }





// // [[Rcpp::export]]
// Rcpp::IntegerVector accessibleGenotypes(Rcpp::IntegerMatrix y,
// 					Rcpp::NumericVector f,
// 					Rcpp::IntegerVector numMut, //
// 					double th) {
  
//   // Return just the indices. Could preserve fitness, but would need
//   // another matrix.
//   int ng = y.nrow(); //it counts the wt
//   Rcpp::IntegerMatrix adm(ng, ng);

//   adm = integerAdjMat(y, f, numMut, th);
  
//   int numMutdiff = 0;

//   // Slightly different logic from R: Do not resize object; set the row to
//   // 0.
//   int colsum = 0;
//   // int indicator = 0; // indicator != 0 means we set one row to 0
//   // so we need to iterate at least once more.
  
//   // accessible is the genotype number, not the column!  WT is 1,
//   // etc. This makes it easy to keep track of which are accessible.
//   Rcpp::IntegerVector accessible = Rcpp::seq_len(ng);

//   // This is doable in one pass
//   // while (true) {
//   //   indicator = 0;
//     for(int k = 1; k < ng; ++k) {
//       if(accessible(k) > 0) {
// 	colsum = std::accumulate(adm(Rcpp::_, k).begin(),
// 				 adm(Rcpp::_, k).end(), 0);
// 	if(colsum == 0) { // This genotype ain't reachable
// 	  // Nothing can be reached from this genotype; fill with 0.
// 	  adm(k, Rcpp::_) = Rcpp::IntegerVector(ng);
// 	  accessible(k) = -9;
// 	  // indicator = 1;
// 	}
//       }
//     }
//   //   if(indicator == 0) break;
//   // }
//   return accessible;
// }




// // [[Rcpp::export]]
// Rcpp::IntegerVector peaksLandscape(Rcpp::IntegerMatrix y,
// 				   Rcpp::NumericVector f,
// 				   Rcpp::IntegerVector numMut, //
// 				   double th) {
//   // Return the indices. This is like accessibleGenotypes, but we do an
//   // extra loop
//   int ng = y.nrow(); //it counts the wt
//   Rcpp::IntegerMatrix adm(ng, ng);

//   adm = integerAdjMat(y, f, numMut, th);
  
//   int numMutdiff = 0;

//   // Slightly different logic from R: Do not resize object; set the row to
//   // 0.
//   int colsum = 0;
//   // int indicator = 0; // indicator != 0 means we set one row to 0
//   // so we need to iterate at least once more.
  
//   // accessible is the genotype number, not the column!  WT is 1,
//   // etc. This makes it easy to keep track of which are accessible.
//   Rcpp::IntegerVector accessible = Rcpp::seq_len(ng);
//   // This is doable in one pass
//   // while (true) {
//   //   indicator = 0;
//     for(int k = 1; k < ng; ++k) {
//       if(accessible(k) > 0) {
// 	colsum = std::accumulate(adm(Rcpp::_, k).begin(),
// 				 adm(Rcpp::_, k).end(), 0);
// 	if(colsum == 0) { // This genotype ain't reachable
// 	  // Nothing can be reached from this genotype; fill with 0.
// 	  adm(k, Rcpp::_) = Rcpp::IntegerVector(ng);
// 	  accessible(k) = -9;
// 	  // indicator = 1;
// 	}
//       }
//     }
//   //   if(indicator == 0) break;
//   // }

//     int rowsum = 0;
//     Rcpp::IntegerVector peaks;
//     for(int k = 1; k < ng; ++k) {
//       if(accessible(k) > 0) {
// 	rowsum = std::accumulate(adm(k, Rcpp::_).begin(),
// 				 adm(k, Rcpp::_).end(), 0);
// 	if(rowsum == 0) { // This genotype doesn't have children
// 	  peaks.push_back(k);
// 	}
//       }
//     }

   
//   return peaks;
// }
