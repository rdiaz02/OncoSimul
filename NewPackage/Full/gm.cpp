#include <Rcpp.h>

using namespace Rcpp ;


struct geneDeps {
  // int typeDep; // smaller, predictable size. A lot less readable, though.
  std::string typeDep; 
  double s;
  double sh;
  // std::vector<int> deps;
  std::vector<std::string > deps; // as module names
};


struct Poset_and_Modules {
  std::vector<geneDeps> Poset;
  std::vector<std::pair<std::string, std::string> > geneToModule;
};

static void rGM_GeneModule(Rcpp::DataFrame rGM,
			   std::vector<std::pair<std::string, std::string> >& geneModule){
  Rcpp::IntegerVector id = rGM["NumericID"];
  Rcpp::CharacterVector Gene = rGM["Gene"];
  Rcpp::CharacterVector Module = rGM["Module"];
  geneModule.resize(id.size());

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != id[i])
      throw std::logic_error(" i != id");
    geneModule[i].first = Gene[i];
    geneModule[i].second = Module[i];
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
  Poset[0].typeDep = "NA";
  Poset[0].s = std::numeric_limits<double>::quiet_NaN();
  Poset[0].sh = std::numeric_limits<double>::quiet_NaN();
  Poset[0].deps.resize(0);
    

  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::CharacterVector module;

  for(int i = 1; i != (rt.size() + 1); ++i) {
    rt_element = rt[i - 1];
    // if(as<int>(rt_element["child"]) != i) {
    //   // Rcpp::Rcout << "\n child = " << as<int>(rt_element["child"]);
    //   // Rcpp::Rcout << "\n i = " << i << std::endl;
    //   throw std::logic_error("child != index");
    // }

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
    Rcpp::Rcout <<"\n\t\t Dependent Module or gene " << i << std::endl;
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

