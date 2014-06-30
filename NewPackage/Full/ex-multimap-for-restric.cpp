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

// stretching vector of mutatedDrv unlikely to speed up;
// access (seeing if a module gene in there) is O(1) but I do
// that for every gene and then sum over the vector (linear in
// number of positions?)

// Could be made much faster if we could assume lower numbered
// positions cannot depend on higher numbered ones. Nope, not that much.
static void checkThisConstraint(const std::vector<int>& mutatedDrv,
				const geneDeps& thisDeps,
				std::vector<double>& s_vector,
				std::vector<double>& sh_vector) {
  // const std::vector<geneDeps>& restrictTable) {
  size_t numDepends;
  size_t sumDependsMet = 0;
  size_t sizeModule = 0;
  int module_mutated = 0;
  if( (thisDeps.deps.size() == 1) &&
      (thisDeps.deps[0][0] == 0) ) { //Depends only on root
    s_vector.push_back(thisDeps.s);
  } else {
    numDepends = thisDeps.deps.size();
    for(size_t i = 0; i != numDependencies; ++i) {
      // any of the module entries are mutated?
      for(std::vector<int>::const_iterator g = thisDeps.deps[i].begin();
	  g != thisDeps.deps[i].end(); ++g) {
	// if sorted, could use binary search
	module_mutated = (std::find(mutatePositions.begin(), 
				    mutatedDrv.end(), 
				    (*g)) != mutatedDrv.end());
	if(module_mutated) {
	  ++sumDependensMet;
	  break;
	}
      }
    }
    if( ((thisDeps.type == "SM") && (sumDependensMet)) ||
	((thisDeps.type == "MN") && (sumDependensMet == numDepends)) ) {
      s_vector.push_back(thisDeps.s);
    } else {
      sh_vector.push_back(thisDeps.sh);
    }
  }
}


// Now, for each clone we have to keep vector s, sh, sp
static void checkConstraints(const int& numDrivers,
			     const std::vector<geneDeps>& restrictTable,
			     const std::vector<int>& mutatedDrv,
			     const Genotype64& newGenotype) {
  for(std::vector<int>::const_iterator mg = mutatedDrv.begin();
      mg != mutatedDrv.end(); ++mg) {
    checkThisConstraint(mutatedDrv, restrictTable[(*mg)],
			s_vector, sh_vector);
  }
}

// FIXME: be very careful with numberings. Do mutations start at 0?
// No, cannot.

// when mutation happens, if a passenger, insert in mutatedPass and give p_vector.
  // if(mutatedPos >= numDrivers) { //the new mutation is a passenger
  //   // do something: iterate through the passenger part of genotype
  //   // keeping the rest of the parent stuff.
  //   return;
  // } else {


// if driver, insert in mutatedDrv.
// and check restrictions

