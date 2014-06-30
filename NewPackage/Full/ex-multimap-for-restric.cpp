#include <Rcpp.h>

using namespace Rcpp ;


// struct Deps {
//   std::vector<int> module;
//   double s;
//   double sh;
// };

// std::multimap<int, Deps> restrictTable;


// FIXME: check later penalty for using a string for typeDep
struct geneDeps {
  int typeDep; // smaller, predictable size. A lot less readable, though.
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

void printRestrictTable(const std::vector<geneDeps>& restrictTable) {
  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" << std::endl;
  Rcpp::Rcout << "\t Size = " << restrictTable.size() << std::endl;
  for(size_t i = 0; i != restrictTable.size(); ++i) {
    Rcpp::Rcout <<"\n\t\t Dependent node " << (i + 1) << std::endl;
    Rcpp::Rcout <<"\t\t\t typeDep = " << restrictTable[i].typeDep << " ";
    Rcpp::Rcout <<"\t s = " << restrictTable[i].s << " ";
    Rcpp::Rcout <<"\t sh = " << restrictTable[i].sh << std::endl;
    if(isinf(restrictTable[i].sh))
      ++counterInfs;
    if(isinf(restrictTable[i].sh) && (restrictTable[i].sh < 0))
      ++counterNegInfs;
    // here the code for parent modules
    Rcpp::Rcout << "\t\t\t Number of parent modules or genes = " << 
      restrictTable[i].deps.size() << std::endl;
    for(size_t j = 0; j != restrictTable[i].deps.size(); ++j) {
      Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": ";
      for (auto c : restrictTable[i].deps[j])
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

void restrictTable_to_cpp(Rcpp::List rt,
			  std::vector<geneDeps>& restrictTable) { 
  int ndeps;
  restrictTable.resize(rt.size());

  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;

  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    restrictTable[i].typeDep = rt_element["type"];
    restrictTable[i].s = rt_element["s"];
    restrictTable[i].sh = rt_element["sh"];
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    if(isinf(restrictTable[i].s))
      Rcpp::Rcout << "WARNING: at least one s is Infinite" 
		  << std::endl;
    if(isinf(restrictTable[i].sh) && (restrictTable[i].sh > 0))
      Rcpp::Rcout << "WARNING: at least one sh is positive Infinite" 
		  << std::endl;

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      restrictTable[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
    }
  }
}


// [[Rcpp::export]]
void wrap_test_rt(Rcpp::List rtR) {
  std::vector<geneDeps> restrictTable;
  restrictTable_to_cpp(rtR, restrictTable);
  printRestrictTable(restrictTable);
}




//  to be used directly from R, since cannot export: second argument of
//  type unknown to Rcpp // [[Rcpp::export]]

// [[Rcpp::export]]
void restrictTable_to_cpp0(Rcpp::List rt){
  std::vector<geneDeps> restrictTable;
  int ndeps;
  restrictTable.resize(rt.size());
  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;

  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    restrictTable[i].typeDep = rt_element["type"];
    restrictTable[i].s = rt_element["s"];
    restrictTable[i].sh = rt_element["sh"];
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      restrictTable[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
    }
  }

  printRestrictTable(restrictTable);

}


static void checkThisConstraint(const int& thisPos, 
				const std::vector<int>& mutatedPositions,
				const geneDeps& thisDeps) {
  // const std::vector<geneDeps>& restrictTable) {
  size_t numDependencies;
  size_t sumDependenciesMet = 0;
  size_t sizeModule = 0;
  if( (thisDeps.deps.size() == 1) &&
      (thisDeps.deps[0][0] == 0) ) {
    return thisDeps.s;
  } else {
    //FIXME: I am casting here. Is this OK? Doing it twice
    numDependencies = thisDeps.deps.size();
    for(size_t i = 0; i != numDependencies; ++i) {
      sizeModule = thisDeps.deps[i].size();
      
    }
    
  }

}




static void checkConstraints(const int& mutatedPos, 
			     const int& numDrivers,
			     const std::vector<geneDeps>& restrictTable,
			     // const std::string& typeCBN,
			     const Genotype64& newGenotype) {
  //      **** Are driver constraints met? ***
  using namespace Rcpp;

  int numDependencies;
  int sumDriversMet = 0;
  int sumDriversNoMet = 0;
  int sumDependenciesMet = 0;

  // Two cases: same s, sh, sp or different ones. If same, return three
  // integers: sumDriversMet, sumDriversNoMet, sumPassengers.  If
  // different, return three vectors, filled with the non-zero
  // entries. These vectors then are combined as dictated by the fintness
  // functions.

  // If same single s, sh, sp: function takes three integers. O.w. it
  // takes three integer vectors.
  

  if(mutatedPos >= numDrivers) { //the new mutation is a passenger
    // do something: iterate through the passenger part of genotype.
    return;
  } else {
    for(int m = 0; m < numDrivers; ++m) {
      if( newGenotype[m] ) { // this m is mutated
	const Rcpp::IntegerMatrix::Column thisRestrict = 
	  restrictTable(_, m);
	numDependencies = thisRestrict[1];
	if(!numDependencies) { // this driver has no dependencies
	  sumDriversMet++;
#ifdef DEBUGZ
	  Rcpp::Rcout << "\n No dependencies:  ";
	  DP2(sumDriversMet);
#endif
	}
	else {
	  sumDependenciesMet = 0;
	  for(int i = 2; i < (2 + numDependencies); i++) {
	    sumDependenciesMet += newGenotype[ thisRestrict[i] ];
	  }
	  if( ( (typeCBN == "Multiple") && (sumDependenciesMet) ) ||
	      ( (typeCBN == "CBN") && (sumDependenciesMet == numDependencies) )) {
	    sumDriversMet++;   
	  } else {
	    sumDriversNoMet++;
	  }
	}
      }
    }
  }

#ifdef DEBUGZ
  DP2(sumDriversMet);
  DP2(sumDriversNoMet);
  DP2(sh);
  DP2(typeFitness);
#endif
}
