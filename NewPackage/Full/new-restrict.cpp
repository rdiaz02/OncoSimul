#include <Rcpp.h>

using namespace Rcpp ;




struct geneDeps {
  // int typeDep; // smaller, predictable size. A lot less readable, though.
  std::string Name;
  std::string typeDep; 
  double s;
  double sh;
  // std::vector<int> deps;
  std::vector<std::string > deps; // as module names
};

struct geneToModule {
  int Gene;
  int Module;
};

struct geneToModuleLong {
  int Gene;
  int Module;
  std::string GeneName;
  std::string ModuleName;
};

struct Poset_and_Modules {
  std::vector<geneDeps> Poset;
  std::vector<geneToModule> geneModule;
};

// change this: we have to numeric ids and two strings
// two things to keep:
// the numeric IDs pair and the strings IDs pair.
// the second only used for checking, not for real
// a vector of tuple or a vector of struct

static void rGM_GeneModule(Rcpp::DataFrame rGM,
			   std::vector<geneToModule> >& geneModule,
			   std::vector<geneToModuleLong>& geneModuleLong,
){
  Rcpp::IntegerVector GeneID = rGM["GeneNumID"];
  Rcpp::IntegerVector ModuleID = rGM["ModuleNumID"];
  Rcpp::CharacterVector Gene = rGM["Gene"];
  Rcpp::CharacterVector Module = rGM["Module"];
  geneModule.resize(id.size());
  geneModuleLong.resize(id.size());

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != GeneID[i])
      throw std::logic_error(" i != GeneID");
    geneModule[i].Gene = GeneID[i];
    geneModule[i].Module = ModuleID[i];
    geneModuleLong[i].Gene = GeneID[i];
    geneModuleLong[i].Module = ModuleID[i];
    geneModuleLong[i].GeneName = Gene[i];
    geneModuleLong[i].ModuleName = Module[i];
  }
}


static void printGeneToModule(const 
			      std::vector<std::pair<std::string, std::string> >& geneModule) {
  Rcpp::Rcout << "\n\ngeneModule table: row number, Gene, Module\n";
  int i = 0;
  for(auto it = geneModule.begin(); it != geneModule.end(); ++it) {
    Rcpp::Rcout << i << ' ' << it->first << ' ' << it->second << std::endl;
    ++i;
  }
}

// For users: if something depends on 0, that is it. No further deps.
// And do not touch the 0 in GeneToModule.

static void rTable_to_Poset(Rcpp::List rt,
			    std::vector<geneDeps>& Poset) { 
  int ndeps;
  // The restriction table, or Poset, has a first element
  // with nothing, so that all references by mutated gene
  // are simply accessing the Poset[mutated gene] without
  // having to remember to add 1, etc. Or things like that.
  Poset.resize(rt.size() + 1);
  Poset[0].Name = "0";
  Poset[0].typeDep = "NA";
  Poset[0].s = std::numeric_limits<double>::quiet_NaN();
  Poset[0].sh = std::numeric_limits<double>::quiet_NaN();
  Poset[0].deps.resize(0);
    

  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::CharacterVector module;

  for(int i = 1; i != (rt.size() + 1); ++i) {
    rt_element = rt[i - 1];
    // FIXME: check this again!!! we are now using numeric IDs
    // if(as<int>(rt_element["child"]) != i) {
    //   // Rcpp::Rcout << "\n child = " << as<int>(rt_element["child"]);
    //   // Rcpp::Rcout << "\n i = " << i << std::endl;
    //   throw std::logic_error("child != index");
    // }
    Poset[i].Name = as<std::string>(rt_element["child"]);
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
      module = as<CharacterVector>(parent_list[j]);
      Poset[i].deps.push_back(Rcpp::as< std::string >(module));
    }
  }
}



void printPoset(const std::vector<geneDeps>& Poset) {
  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" << std::endl;
  Rcpp::Rcout << "\t Size = " << (Poset.size() - 1) << std::endl;
  for(size_t i = 1; i != Poset.size(); ++i) {
    // We do not show the Poset[0]
    Rcpp::Rcout <<"\n\t\t Dependent Module or gene " << i 
		<< ". Name: " << Poset[i].Name << std::endl;
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
      Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": " 
		  << Poset[i].deps[j] << std::endl;
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

// [[Rcpp::export]]
void wrap_test_rt(Rcpp::List rtR,
		  Rcpp::DataFrame rGM) {
  Poset_and_Modules PosetModule;
  rGM_GeneModule(rGM, PosetModule.geneToModule);
  
  printGeneToModule(PosetModule.geneToModule);
  rTable_to_Poset(rtR, PosetModule.Poset);
  printPoset(PosetModule.Poset);
}




// stretching vector of mutatedDrv unlikely to speed up;
// access (seeing if a module gene in there) is O(1) but I do
// that for every gene and then sum over the vector (linear in
// number of positions?)

// Could be made much faster if we could assume lower numbered
// positions cannot depend on higher numbered ones. Nope, not that much.


static void DrvToModule(const std::vector<int>& Drv,
			const 
			std::vector<std::pair<std::string, std::string> >& geneToModule,
			std::vector<std::string>& mutatedModules) {
  
  for(auto it = Drv.begin(); it != Drv.end(); ++it) {
    mutatedModules.push_back(geneModule[(*it)].second);
  }
  sort( mutatedModules.begin(), mutatedModules.end() );
  mutatedModules.erase( unique( mutatedModules.begin(), 
				mutatedModules.end() ), 
			mutatedModules.end() );
}


static void checkConstraints(const std::vector<int>& Drv,
			     Poset_and_Modules& PosetModule, 
    //			     const std::vector<geneDeps>& Poset,
			     std::vector<double>& s_vector,
			     std::vector<double>& sh_vector) {
  size_t numDeps;
  size_t sumDepsMet = 0;
  int module_mutated = 0;

  std::vector<std::string> mutatedModules;
  // FIXME
  // Think this!! Use the numericIDs always!!!
  // Only at end we go back from Genes to IDs (or maybe in R)

  DrvToModules(Drv, PosetModule.geneToModule, mutatedModules);


  // for(std::vector<int>::const_iterator gene = Drv.begin();
  //     gene != Drv.end(); ++gene) {

  for(auto module = mutatedModules.begin();
      module != mutatedModules.end(); ++module) {
    if( (PosetModule.Poset[(*module)].deps.size() == 1) &&
	(PosetModule.Poset[(*module)].deps[0] == "0") ) { //Depends only on root
      s_vector.push_back(Poset[(*module)].s);
    } else {
      sumDepsMet = 0;
      numDeps = Poset[(*gene)].deps.size();
      // for(std::vector<std::vector<int> >::const_iterator Module = Poset[(*gene)].deps.begin();
      // 	  Module !=  Poset[(*gene)].deps.end(); ++Module) {
      // and then replace Poset[(*gene)].deps[i] by Module
      for(size_t i = 0; i != numDeps; ++i) {
	// any of the module entries are mutated?
	for(std::vector<int>::const_iterator m = Poset[(*gene)].deps[i].begin();
	    m != Poset[(*gene)].deps[i].end(); ++m) {
	  // if sorted, could use binary search
	  module_mutated = (std::find(Drv.begin(), 
				      Drv.end(), 
				      (*m)) != Drv.end());
	  if(module_mutated) {
	    ++sumDepsMet;
	    break;
	  }
	}
      }
      if( ((Poset[(*gene)].typeDep == "SM") && (sumDepsMet)) ||
	  ((Poset[(*gene)].typeDep == "MN") && (sumDepsMet == numDeps)) ) {
	s_vector.push_back(Poset[(*gene)].s);
      } else {
	sh_vector.push_back(Poset[(*gene)].sh);
      }
    }
  }
}



// // [[Rcpp::export]]
// SEXP wrap_test_checkRestriction(Rcpp::List rtR, 
// 				Rcpp::IntegerVector genotype) {
//   std::vector<geneDeps> Poset;
//   std::vector<int> Genotype;
//   std::vector<double> s_vector;
//   std::vector<double> sh_vector;

//   rTable_to_Poset(rtR, Poset);
//   Genotype = Rcpp::as< std::vector<int> >(genotype);

//   checkConstraints(Genotype, Poset, s_vector, sh_vector);

//   // Rcpp::Rcout << "s_vector:\n ";
//   // for (auto c : s_vector)
//   //   Rcpp::Rcout << c << ' ';
//   // Rcpp::Rcout << std::endl;
//   // Rcpp::Rcout << "sh_vector:\n ";
//   // for (auto c : sh_vector)
//   //   Rcpp::Rcout << c << ' ';
//   // Rcpp::Rcout << std::endl;

//   return
//     List::create(Named("s_vector") = s_vector,
// 		 Named("sh_vector") = sh_vector);
// }


// when mutation happens, if a passenger, insert in mutatedPass and give p_vector.
  // if(mutatedPos >= numDrivers) { //the new mutation is a passenger
  //   // do something: iterate through the passenger part of genotype
  //   // keeping the rest of the parent stuff.
  //   return;
  // } else {


// if driver, insert in Drv.
// and check restrictions








// static void rTable_to_Poset(Rcpp::List rt,
// 			  std::vector<geneDeps>& Poset) { 
//   int ndeps;
//   // The restriction table, or Poset, has a first element
//   // with nothing, so that all references by mutated gene
//   // are simply accessing the Poset[mutated gene] without
//   // having to remember to add 1, etc.
//   Poset.resize(rt.size() + 1);
//   Poset[0].typeDep = "NA";
//   Poset[0].s = std::numeric_limits<double>::quiet_NaN();
//   Poset[0].sh = std::numeric_limits<double>::quiet_NaN();
//   Poset[0].deps.resize(0);
    

//   Rcpp::List rt_element;
//   Rcpp::List parent_list;
//   Rcpp::IntegerVector module;

//   for(int i = 1; i != (rt.size() + 1); ++i) {
//     rt_element = rt[i - 1];
//     if(as<int>(rt_element["child"]) != i) {
//       // Rcpp::Rcout << "\n child = " << as<int>(rt_element["child"]);
//       // Rcpp::Rcout << "\n i = " << i << std::endl;
//       throw std::logic_error("child != index");
//     }
//     Poset[i].typeDep = as<std::string>(rt_element["typeDep"]);
//     Poset[i].s = as<double>(rt_element["s"]);
//     Poset[i].sh = as<double>(rt_element["sh"]);
//     parent_list = rt_element["parent"];
//     ndeps = parent_list.size();

//     if(isinf(Poset[i].s))
//       Rcpp::Rcout << "WARNING: at least one s is infinite" 
// 		  << std::endl;
//     if(isinf(Poset[i].sh) && (Poset[i].sh > 0))
//       Rcpp::Rcout << "WARNING: at least one sh is positive infinite" 
// 		  << std::endl;

//     for(int j = 0; j != ndeps; ++j) {
//       module = as<IntegerVector>(parent_list[j]);
//       Poset[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
//     }
//   }
// }




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

// void printPoset(const std::vector<geneDeps>& Poset) {
//   int counterInfs = 0;
//   int counterNegInfs = 0;
//   Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" << std::endl;
//   Rcpp::Rcout << "\t Size = " << (Poset.size() - 1) << std::endl;
//   for(size_t i = 1; i != Poset.size(); ++i) {
//     // We do not show the Poset[0]
//     Rcpp::Rcout <<"\n\t\t Dependent gene " << i << std::endl;
//     Rcpp::Rcout <<"\t\t\t typeDep = " << Poset[i].typeDep << " ";
//     Rcpp::Rcout <<"\t s = " << Poset[i].s << " ";
//     Rcpp::Rcout <<"\t sh = " << Poset[i].sh << std::endl;
//     if(isinf(Poset[i].sh))
//       ++counterInfs;
//     if(isinf(Poset[i].sh) && (Poset[i].sh < 0))
//       ++counterNegInfs;
//     // here the code for parent modules
//     Rcpp::Rcout << "\t\t\t Number of parent modules or genes = " << 
//       Poset[i].deps.size() << std::endl;
//     for(size_t j = 0; j != Poset[i].deps.size(); ++j) {
//       Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": ";
//       for (auto c : Poset[i].deps[j])
//        	Rcpp::Rcout << c << ' ';
//       Rcpp::Rcout << std::endl;
//     }
//     Rcpp::Rcout << std::endl;
//   }
//   if(counterInfs) {
//     Rcpp::Rcout << "In sh there were " << counterNegInfs 
// 	     << " negative infinites and "
// 	     << (counterInfs - counterNegInfs) 
// 	     << " positive infinites" << std::endl;
//   }
// }

























