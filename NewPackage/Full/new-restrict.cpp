#include <Rcpp.h>

using namespace Rcpp ;


// struct Deps {
//   std::vector<int> module;
//   double s;
//   double sh;
// };

// std::multimap<int, Deps> Poset;


// FIXME: check later penalty for using a string for type
struct geneDeps {
  // int typeDep; // smaller, predictable size. A lot less readable, though.
  std::string typeDep; // smaller, predictable size. A lot less readable, though.
  double s;
  double sh;
  // std::vector<int> deps;
  std::vector< std::vector<int> > deps;
};

// // [[Rcpp::export]]
// void f4(){
//   std::vector<geneDeps> g3;
//   g3.resize(4);
//   g3[0].s = 0.3;
//   g3[0].sh = 0.6;
//   g3[0].deps.resize(2);
//   // g3[0].deps[0].push_back(std::vector<int>(3));
//   std::cout << "g3[0].deps.size() " << g3[0].deps.size() << std::endl;
//   std::cout << "g3[0].deps[0].size() " << 
//     g3[0].deps[0].size() << std::endl;
//   g3[0].deps[0].push_back(3);
//   std::cout << "g3[0].deps[0].size() " << 
//     g3[0].deps[0].size() << std::endl;

//   std::vector<int>vv(4, 3);
//   g3[1].deps.push_back(vv);
//   std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
//   g3[1].deps.push_back(std::vector<int>(9, 5));
//   std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
  
//   for (auto c : g3[1].deps[0])
//     std::cout << c << ' ';
  
//   for (auto c : g3[1].deps[1])
//     std::cout << c << ' ';
// }

void printPoset(const std::vector<geneDeps>& Poset) {
  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" << std::endl;
  Rcpp::Rcout << "\t Size = " << Poset.size() << std::endl;
  for(size_t i = 0; i != Poset.size(); ++i) {
    Rcpp::Rcout <<"\n\t\t Dependent gene " << (i + 1) << std::endl;
    Rcpp::Rcout <<"\t\t\t typeDep = " << Poset[i].typeDep << " ";
    Rcpp::Rcout <<"\t s = " << Poset[i].s << " ";
    Rcpp::Rcout <<"\t sh = " << Poset[i].sh << std::endl;
    if(isinf(Poset[i].sh))
      ++counterInfs;
    if(isinf(Poset[i].sh) && (Poset[i].sh < 0))
      ++counterNegInfs;
    // here the code for parent modules
    Rcpp::Rcout << "\t\t\t Number of parent modules or genes = " << 
      Poset[i].deps.size() << std::endl;
    for(size_t j = 0; j != Poset[i].deps.size(); ++j) {
      Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": ";
      for (auto c : Poset[i].deps[j])
       	Rcpp::Rcout << c << ' ';
      Rcpp::Rcout << std::endl;
    }
    Rcpp::Rcout << std::endl;
  }
  if(counterInfs) {
    Rcpp::Rcout << "In sh there were " << counterNegInfs 
	     << " negative infinites and "
	     << (counterInfs - counterNegInfs) 
	     << " positive infinites" << std::endl;
  }
}

void rTable_to_Poset(Rcpp::List rt,
			  std::vector<geneDeps>& Poset) { 
  int ndeps;
  Poset.resize(rt.size());

  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;

  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    Poset[i].typeDep = as<std::string>(rt_element["typeDep"]);
    Poset[i].s = as<double>(rt_element["s"]);
    Poset[i].sh = as<double>(rt_element["sh"]);
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    if(isinf(Poset[i].s))
      Rcpp::Rcout << "WARNING: at least one s is infinite" 
		  << std::endl;
    if(isinf(Poset[i].sh) && (Poset[i].sh > 0))
      Rcpp::Rcout << "WARNING: at least one sh is positive infinite" 
		  << std::endl;

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      Poset[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
    }
  }
}


// [[Rcpp::export]]
void wrap_test_rt(Rcpp::List rtR) {
  std::vector<geneDeps> Poset;
  rTable_to_Poset(rtR, Poset);
  printPoset(Poset);
}




//  to be used directly from R, since cannot export: second argument of
//  type unknown to Rcpp // [[Rcpp::export]]

// // [[Rcpp::export]]
// void rTable_to_Poset0(Rcpp::List rt){
//   std::vector<geneDeps> Poset;
//   int ndeps;
//   Poset.resize(rt.size());
//   Rcpp::List rt_element;
//   Rcpp::List parent_list;
//   Rcpp::IntegerVector module;

//   for(int i = 0; i != rt.size(); ++i) {
//     rt_element = rt[i];
//     Poset[i].typeDep = as<std::string>(rt_element["typeDep"]);
//     Poset[i].s = as<double>(rt_element["s"]);
//     Poset[i].sh = as<double>(rt_element["sh"]);
//     parent_list = rt_element["parent"];
//     ndeps = parent_list.size();

//     for(int j = 0; j != ndeps; ++j) {
//       module = as<IntegerVector>(parent_list[j]);
//       Poset[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
//     }
//   }
//   // printPoset(Poset);

// }

// stretching vector of mutatedDrv unlikely to speed up;
// access (seeing if a module gene in there) is O(1) but I do
// that for every gene and then sum over the vector (linear in
// number of positions?)

// Could be made much faster if we could assume lower numbered
// positions cannot depend on higher numbered ones. Nope, not that much.
// static void checkConstraints(const std::vector<int>& Drv,
// 			     const std::vector<geneDeps>& Poset,
// 			     std::vector<double>& s_vector,
// 			     std::vector<double>& sh_vector) {
//   size_t numDeps;
//   size_t sumDepsMet = 0;
//   int module_mutated = 0;

//   for(std::vector<int>::const_iterator gene = Drv.begin();
//       gene != Drv.end(); ++gene) {
//     if( (Poset[(*gene)].deps.size() == 1) &&
// 	(Poset[(*gene)].deps[0][0] == 0) ) { //Depends only on root
//       s_vector.push_back(Poset[(*gene)].s);
//     } else {
//       numDeps = Poset[(*gene)].deps.size();
//       // for(std::vector<std::vector<int> >::const_iterator Module = Poset[(*gene)].deps.begin();
//       // 	  Module !=  Poset[(*gene)].deps.end(); ++Module) {
//       // and then replace Poset[(*gene)].deps[i] by Module
//       for(size_t i = 0; i != numDeps; ++i) {
// 	// any of the module entries are mutated?
// 	for(std::vector<int>::const_iterator m = Poset[(*gene)].deps[i].begin();
// 	    m != Poset[(*gene)].deps[i].end(); ++m) {
// 	  // if sorted, could use binary search
// 	  module_mutated = (std::find(mutatePositions.begin(), 
// 				      Drv.end(), 
// 				      (*m)) != Drv.end());
// 	  if(module_mutated) {
// 	    ++sumDepsMet;
// 	    break;
// 	  }
// 	}
//       }
//       if( ((Poset[(*gene)].type == "SM") && (sumDepsMet)) ||
// 	  ((Poset[(*gene)].type == "MN") && (sumDepsMet == numDeps)) ) {
// 	s_vector.push_back(Poset[(*gene)].s);
//       } else {
// 	sh_vector.push_back(Poset[(*gene)].sh);
//       }
//     }
//   }
// }



// FIXME: be very careful with numberings. Do mutations start at 0?
// No, cannot.

// when mutation happens, if a passenger, insert in mutatedPass and give p_vector.
  // if(mutatedPos >= numDrivers) { //the new mutation is a passenger
  //   // do something: iterate through the passenger part of genotype
  //   // keeping the rest of the parent stuff.
  //   return;
  // } else {


// if driver, insert in Drv.
// and check restrictions

