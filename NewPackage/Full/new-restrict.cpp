#include <Rcpp.h>

using namespace Rcpp ;

enum class Dependency {monotone, semimonotone, NA}; // monotone, semimonotone

inline static Dependency stringToDep(const std::string& dep) {
  if(dep == "monotone")
    return Dependency::monotone;
  else if(dep == "semimonotone")
    return Dependency::semimonotone;
  else 
    throw std::out_of_range("Not a valid typeDep");
}

inline static std::string depToString(const Dependency dep) {
  switch(dep) {
  case Dependency::monotone:
    return "monotone";
  case Dependency::semimonotone:
    return "semimonotone";
  case Dependency::NA:
    return "NA";
  default:
    throw std::out_of_range("Not a valid dependency");
  }
}

struct geneToModule {
  int GeneID;
  int ModuleID;
};

// this is temporary. Will remove at end??
struct geneToModuleLong {
  int GeneID;
  int ModuleID;
  std::string GeneName;
  std::string ModuleName;
};

struct geneDeps {
  Dependency typeDep;
  int childID; //redundant, but leave
  double s;
  double sh;
  std::vector<int> parentsID;
  // The next two are clearly redundant but a triple check
  std::string child;
  std::vector<std::string> parents;
};


// For users: if something depends on 0, that is it. No further deps.
// And do not touch the 0 in GeneToModule.

static void rTable_to_Poset(Rcpp::List rt,
			    std::vector<geneDeps>& Poset) { 

  // The restriction table, or Poset, has a first element
  // with nothing, so that all references by mutated gene
  // are simply accessing the Poset[mutated gene] without
  // having to remember to add 1, etc. 
  Poset.resize(rt.size() + 1);
  Poset[0].child = "0"; //should this be Root?? I don't think so.
  Poset[0].childID = 0;
  Poset[0].typeDep = Dependency::NA;
  Poset[0].s = std::numeric_limits<double>::quiet_NaN();
  Poset[0].sh = std::numeric_limits<double>::quiet_NaN();
  Poset[0].parents.resize(0);
  Poset[0].parentsID.resize(0);

  Rcpp::List rt_element;
  // std::string tmpname;
  Rcpp::IntegerVector parentsid;
  Rcpp::StringVector parents;

  for(int i = 1; i != (rt.size() + 1); ++i) {
    rt_element = rt[i - 1];
    Poset[i].child = Rcpp::as<std::string>(rt_element["child"]);
    Poset[i].childID = as<int>(rt_element["childID"]);
    Poset[i].typeDep = stringToDep(as<std::string>(rt_element["typeDep"]));
    Poset[i].s = as<double>(rt_element["s"]);
    Poset[i].sh = as<double>(rt_element["sh"]);

    if(i != Poset[i].childID) {
      // Rcpp::Rcout << "\n childID, original = " << as<int>(rt_element["childID"]);
      // Rcpp::Rcout << "\n childID, Poset = " << Poset[i].childID;
      // Rcpp::Rcout << "\n i = " << i << std::endl;
      throw std::logic_error("childID != index");
    }
    // Rcpp::IntegerVector parentsid = as<Rcpp::IntegerVector>(rt_element["parentsID"]);
    // Rcpp::StringVector parents = as<Rcpp::StringVector>(rt_element["parents"]);

    parentsid = as<Rcpp::IntegerVector>(rt_element["parentsID"]);
    parents = as<Rcpp::StringVector>(rt_element["parents"]);

    if( parentsid.size() != parents.size() ) {
      throw std::logic_error("parents size != parentsID size");
    }

    for(int j = 0; j != parentsid.size(); ++j) {
      Poset[i].parentsID.push_back(parentsid[j]);
      Poset[i].parents.push_back( (Rcpp::as< std::string >(parents[j])) );
      // tmpname = Rcpp::as< std::string >(parents[j]);
      // Poset[i].parents.push_back(tmpname);
    }

    if(isinf(Poset[i].s))
      Rcpp::Rcout << "WARNING: at least one s is infinite" 
		  << std::endl;
    if(isinf(Poset[i].sh) && (Poset[i].sh > 0))
      Rcpp::Rcout << "WARNING: at least one sh is positive infinite" 
		  << std::endl;
  }
}



void printPoset(const std::vector<geneDeps>& Poset) {

  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" 
	      << std::endl;
  Rcpp::Rcout << "Size = " << (Poset.size() - 1) << std::endl;
  for(size_t i = 1; i != Poset.size(); ++i) {
    // We do not show the Poset[0]
    Rcpp::Rcout <<"\t Dependent Module or gene (child) " << i 
		<< ". childID: " << Poset[i].childID 
		<< ". child full name: " << Poset[i].child
		<< std::endl;
    Rcpp::Rcout <<"\t\t typeDep = " << depToString(Poset[i].typeDep) << ' ' ;
    Rcpp::Rcout <<"\t s = " << Poset[i].s << " ";
    Rcpp::Rcout <<"\t sh = " << Poset[i].sh << std::endl;
    if(isinf(Poset[i].sh))
      ++counterInfs;
    if(isinf(Poset[i].sh) && (Poset[i].sh < 0))
      ++counterNegInfs;
    Rcpp::Rcout << "\t\t Number of parent modules or genes = " << 
      Poset[i].parents.size() << std::endl;
    Rcpp::Rcout << "\t\t\t Parents IDs: ";
    for(auto c : Poset[i].parentsID)
      Rcpp::Rcout << c << "; ";
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "\t\t\t Parents names: ";
    for(auto c : Poset[i].parents)
      Rcpp::Rcout << '(' << c << ')' << "; ";
    Rcpp::Rcout << std::endl;
    
    // for(size_t j = 0; j != Poset[i].deps.size(); ++j) {
    //   Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": " 
    // 		  << Poset[i].deps[j] << std::endl;
  }
  Rcpp::Rcout << std::endl;

  if(counterInfs) {
    Rcpp::Rcout << "In sh there were " << counterNegInfs 
		<< " negative infinites and "
		<< (counterInfs - counterNegInfs) 
		<< " positive infinites" << std::endl;
  }
}



static void rGM_GeneModule(Rcpp::DataFrame rGM,
			   std::vector<geneToModule>& geneModule,
			   std::vector<geneToModuleLong>& geneModuleLong){

  Rcpp::IntegerVector GeneID = rGM["GeneNumID"];
  Rcpp::IntegerVector ModuleID = rGM["ModuleNumID"];
  Rcpp::CharacterVector GeneName = rGM["Gene"];
  Rcpp::CharacterVector ModuleName = rGM["Module"];
  geneModule.resize(GeneID.size());
  geneModuleLong.resize(GeneID.size()); // remove later

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != GeneID[i])
      throw std::logic_error(" i != GeneID");
    geneModule[i].GeneID = GeneID[i];
    geneModule[i].ModuleID = ModuleID[i];
    // remove these later?
    geneModuleLong[i].GeneID = GeneID[i];
    geneModuleLong[i].ModuleID = ModuleID[i];
    geneModuleLong[i].GeneName = GeneName[i];
    geneModuleLong[i].ModuleName = ModuleName[i];
  }
}


static void printGeneToModule(const 
			      std::vector<geneToModuleLong>& geneModulesLong) {
  Rcpp::Rcout << 
    "\n\n******** geneModule table inside C++ *******:\ngene ID,\t Gene,\t Module ID,\t Module\n";
  for(auto it = geneModulesLong.begin(); it != geneModulesLong.end(); ++it) {
    Rcpp::Rcout << '\t' << it->GeneID << '\t' << it->GeneName << '\t' 
		<< it->ModuleID << '\t' << it->ModuleName << std::endl;
  }
}



// [[Rcpp::export]]
void wrap_test_rt(Rcpp::List rtR, Rcpp::DataFrame rGM) {
  std::vector<geneDeps> Poset;
  std::vector<geneToModule> geneModules;
  std::vector<geneToModuleLong> geneModulesLong;

  rGM_GeneModule(rGM, geneModules, geneModulesLong);
  printGeneToModule(geneModulesLong);

  rTable_to_Poset(rtR, Poset);
  printPoset(Poset);
}




// // stretching vector of mutatedDrv unlikely to speed up;
// // access (seeing if a module gene in there) is O(1) but I do
// // that for every gene and then sum over the vector (linear in
// // number of positions?)

// // Could be made much faster if we could assume lower numbered
// // positions cannot depend on higher numbered ones. Nope, not that much.


static void DrvToModule(const std::vector<int>& Drv,
			const 
			std::vector<geneToModule>& geneModules,
			std::vector<int>& mutatedModules) {
  
  for(auto it = Drv.begin(); it != Drv.end(); ++it) {
    mutatedModules.push_back(geneModules[(*it)].ModuleID);
  }
  sort( mutatedModules.begin(), mutatedModules.end() );
  mutatedModules.erase( unique( mutatedModules.begin(), 
				mutatedModules.end() ), 
			mutatedModules.end() );
  // That is sorted. So use binary search below. But shouldn't I use a set?
}

// FIXME: can we make it faster if we know each module a single gene?
// FIXME: if genotype is always kept sorted, with drivers first, can it be
// faster? As well, note that number of drivers is automatically known
// from this table of constraints.

static void checkConstraints(const std::vector<int>& Drv,
			     const std::vector<geneDeps>& Poset,
			     const std::vector<geneToModule>& geneToModules,
			     std::vector<double>& s_vector,
			     std::vector<double>& sh_vector) {
  size_t numDeps;
  size_t sumDepsMet = 0;
  int parent_module_mutated = 0;
  std::vector<int> mutatedModules;

  // Map mutated genes to modules (mutatedModules) and then examine, for
  // each mutatedModule, if its dependencies (in terms of modules) are met.

  DrvToModule(Drv, geneToModules, mutatedModules);

  for(auto it_mutatedModule = mutatedModules.begin();
      it_mutatedModule != mutatedModules.end(); ++it_mutatedModule) {
    if( (Poset[(*it_mutatedModule)].parentsID.size() == 1) &&
	(Poset[(*it_mutatedModule)].parentsID[0] == 0) ) { //Depends only on root.
      // FIXME: isn't it enough to check the second condition?
      s_vector.push_back(Poset[(*it_mutatedModule)].s);
    } else {
      sumDepsMet = 0;
      numDeps = Poset[(*it_mutatedModule)].parentsID.size();
      for(auto it_Parents = Poset[(*it_mutatedModule)].parentsID.begin();
	  it_Parents != Poset[(*it_mutatedModule)].parentsID.end();
	  ++it_Parents) {
	// if sorted, could use binary search
	// FIXME: try a set or sort mutatedModules?

	//parent_module_mutated =
	//std::binary_search(mutatedModules.begin(), mutatedModules.end(),
	//(*it_Parents))

	parent_module_mutated = 
	  (std::find(mutatedModules.begin(), 
		     mutatedModules.end(), 
		     (*it_Parents)) != mutatedModules.end());
	if(parent_module_mutated)  {
	  ++sumDepsMet;
	  if( Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone)
	    break;
	}
      }
      if( ((Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone) && 
	   (sumDepsMet)) ||
	  ((Poset[(*it_mutatedModule)].typeDep == Dependency::monotone) && 
	   (sumDepsMet == numDeps)) ) {
	s_vector.push_back(Poset[(*it_mutatedModule)].s);
      } else {
	sh_vector.push_back(Poset[(*it_mutatedModule)].sh);
      }
    }
  }
}





// SEXP wrap_test_checkRestriction(Rcpp::List rtR, 
// 				Rcpp::DataFrame rGM,
// 				Rcpp::IntegerVector genotype) {


// Turn the following into a function called from R naturally.  Allows for
// testing AND allows users to understand the consequences of an rt and a
// genotype. The R function is evalGenotype


// [[Rcpp::export]]
SEXP eval_Genotype(Rcpp::List rtR, 
		   Rcpp::DataFrame rGM,
		   Rcpp::IntegerVector genotype) {

  std::vector<geneDeps> Poset;
  std::vector<geneToModule> geneModules;
  std::vector<geneToModuleLong> geneModulesLong;

  std::vector<int> Genotype;
  std::vector<double> s_vector;
  std::vector<double> sh_vector;

  rTable_to_Poset(rtR, Poset);
  rGM_GeneModule(rGM, geneModules, geneModulesLong);


  Genotype = Rcpp::as< std::vector<int> >(genotype);

  checkConstraints(Genotype, Poset, geneModules,  s_vector, sh_vector);

  // Rcpp::Rcout << "s_vector:\n ";
  // for (auto c : s_vector)
  //   Rcpp::Rcout << c << ' ';
  // Rcpp::Rcout << std::endl;
  // Rcpp::Rcout << "sh_vector:\n ";
  // for (auto c : sh_vector)
  //   Rcpp::Rcout << c << ' ';
  // Rcpp::Rcout << std::endl;

  return
    List::create(Named("s_vector") = s_vector,
		 Named("sh_vector") = sh_vector);
}


// when mutation happens, if a passenger, insert in mutatedPass and give p_vector.
  // if(mutatedPos >= numDrivers) { //the new mutation is a passenger
  //   // do something: iterate through the passenger part of genotype
  //   // keeping the rest of the parent stuff.
  //   return;
  // } else {


// if driver, insert in Drv.
// and check restrictions



































