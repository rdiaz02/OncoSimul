// FIXME later, put all fitness effects in a single structure

#include <Rcpp.h>
#include <iomanip> // for setw
#include <algorithm>    
#include <random>

#define DP2(x) {Rcpp::Rcout << "\n DEBUG2: Value of " << #x << " = " << x << std::endl;}

using namespace Rcpp ;

// FIXME: move this later
int seed = 1; 
std::mt19937 ran_gen(seed);

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
  int GeneNumID;
  int ModuleNumID;
};

// this is temporary. Will remove at end??
struct geneToModuleLong {
  std::string GeneName;
  std::string ModuleName;
  int GeneNumID;
  int ModuleNumID;
};

struct geneDeps {
  Dependency typeDep;
  int childNumID; //redundant, but leave
  double s;
  double sh;
  std::vector<int> parentsNumID;
  // The next two are clearly redundant but a triple check
  std::string child;
  std::vector<std::string> parents;
};


// We use same structure for epistasis and order effects. With order
// effects, NumID is NOT sorted, but reflects the order of the
// restriction. And checking is done using that fact.

struct epistasis {
  double s;
  std::vector<int> NumID; //a set instead? nope.using includes with epistasis
  std::vector<std::string> names; // will remove later
};

// struct orderEffects {
//   std::vector<int> NumID;
//   std::vector<std::string> names; // will remove later
//   double s;
// };

struct genesWithoutInt {
  int shift; // access the s as s[index of mutation or index of mutated
	     // gene in genome - shift]. shift is the min. of NumID, given
	     // how that is numbered from R. We assume mutations always
	     // indexed 1 to something. Not 0 to something.
  // If shift is -9, no elements
  // The next first two are not really needed. Will remove later
  std::vector<int> NumID;
  std::vector<std::string> names;
  std::vector<double> s;
};


struct fitnessEffectsAll {
  bool gMOneToOne;
  // Here to place new mutations in their correct place. Only one is needed.
  // Use the one that is presumably always shorter
  std::vector<int> allOrderG;
  // std::vector<int> allEpistRTG;
  std::vector<geneDeps> Poset;
  std::vector<epistasis> Epistasis;
  std::vector<epistasis> orderE;
  std::vector<geneToModule> geneModules;
  std::vector<geneToModuleLong> geneModulesLong;
  genesWithoutInt genesNoInt;
};


// No, there are no shared genes in order and epist.  Yes, any gene in
// orderEff can also be in the posets, but orderEff is for only for those
// that have order effects.
struct Genotype {
  std::vector<int> orderEff;
  std::vector<int> epistRtEff; //always sorted
  std::vector<int> rest; // always sorted
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
  Poset[0].childNumID = 0;
  Poset[0].typeDep = Dependency::NA;
  Poset[0].s = std::numeric_limits<double>::quiet_NaN();
  Poset[0].sh = std::numeric_limits<double>::quiet_NaN();
  Poset[0].parents.resize(0);
  Poset[0].parentsNumID.resize(0);

  Rcpp::List rt_element;
  // std::string tmpname;
  Rcpp::IntegerVector parentsid;
  Rcpp::CharacterVector parents;

  for(int i = 1; i != (rt.size() + 1); ++i) {
    rt_element = rt[i - 1];
    Poset[i].child = Rcpp::as<std::string>(rt_element["child"]);
    Poset[i].childNumID = as<int>(rt_element["childNumID"]);
    Poset[i].typeDep = stringToDep(as<std::string>(rt_element["typeDep"]));
    Poset[i].s = as<double>(rt_element["s"]);
    Poset[i].sh = as<double>(rt_element["sh"]);

    if(i != Poset[i].childNumID) {
      // Rcpp::Rcout << "\n childNumID, original = " << as<int>(rt_element["childNumID"]);
      // Rcpp::Rcout << "\n childNumID, Poset = " << Poset[i].childNumID;
      // Rcpp::Rcout << "\n i = " << i << std::endl;
      throw std::logic_error("childNumID != index");
    }
    // Rcpp::IntegerVector parentsid = as<Rcpp::IntegerVector>(rt_element["parentsNumID"]);
    // Rcpp::CharacterVector parents = as<Rcpp::CharacterVector>(rt_element["parents"]);

    parentsid = as<Rcpp::IntegerVector>(rt_element["parentsNumID"]);
    parents = as<Rcpp::CharacterVector>(rt_element["parents"]);

    if( parentsid.size() != parents.size() ) {
      throw std::logic_error("parents size != parentsNumID size");
    }

    for(int j = 0; j != parentsid.size(); ++j) {
      Poset[i].parentsNumID.push_back(parentsid[j]);
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


static void rGM_GeneModule(Rcpp::List rGM,
			   std::vector<geneToModule>& geneModule,
			   std::vector<geneToModuleLong>& geneModuleLong){

  Rcpp::IntegerVector GeneNumID = rGM["GeneNumID"];
  Rcpp::IntegerVector ModuleNumID = rGM["ModuleNumID"];
  Rcpp::CharacterVector GeneName = rGM["Gene"];
  Rcpp::CharacterVector ModuleName = rGM["Module"];
  geneModule.resize(GeneNumID.size());
  geneModuleLong.resize(GeneNumID.size()); // remove later

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != GeneNumID[i])
      throw std::logic_error(" i != GeneNumID");
    geneModule[i].GeneNumID = GeneNumID[i];
    geneModule[i].ModuleNumID = ModuleNumID[i];
    // remove these later?
    geneModuleLong[i].GeneNumID = GeneNumID[i];
    geneModuleLong[i].ModuleNumID = ModuleNumID[i];
    geneModuleLong[i].GeneName = GeneName[i];
    geneModuleLong[i].ModuleName = ModuleName[i];
  }
}


static void DrvToModule(const std::vector<int>& Drv,
			const 
			std::vector<geneToModule>& geneModules,
			std::vector<int>& mutatedModules) {
  
  for(auto it = Drv.begin(); it != Drv.end(); ++it) {
    mutatedModules.push_back(geneModules[(*it)].ModuleNumID);
  }
  sort( mutatedModules.begin(), mutatedModules.end() );
  mutatedModules.erase( unique( mutatedModules.begin(), 
				mutatedModules.end() ), 
			mutatedModules.end() );
  // That is sorted. So use binary search below. But shouldn't I use a set?
}


static void convertNoInts(Rcpp::List nI,
			  genesWithoutInt& genesNoInt) {
  
  Rcpp::CharacterVector names = nI["Gene"];
  Rcpp::IntegerVector id = nI["GeneNumID"];
  Rcpp::NumericVector s1 = nI["s"];
  
  genesNoInt.names = Rcpp::as<std::vector<std::string> >(names);
  genesNoInt.NumID = Rcpp::as<std::vector<int> >(nI["GeneNumID"]);
  genesNoInt.s = Rcpp::as<std::vector<double> >(nI["s"]);
  genesNoInt.shift = genesNoInt.NumID[0]; // FIXME: we assume mutations always
				      // indexed 1 to something. Not 0 to
				      // something.
}

static void convertEpiOrderEff(Rcpp::List ep,
			       std::vector<epistasis>& Epistasis,
			       std::vector<int>& allG) {
  Rcpp::List element;
  // For epistasis, the numID must be sorted, but never with order effects.
  // Things come sorted (or not) from R.
  Epistasis.resize(ep.size());
  for(int i = 0; i != ep.size(); ++i) {
    element = ep[i];
    Epistasis[i].NumID = Rcpp::as<std::vector<int> >(element["NumID"]);
    Epistasis[i].names = Rcpp::as<std::vector<std::string> >(element["ids"]);
    Epistasis[i].s = as<double>(element["s"]);
    for(auto g : Epistasis[i].NumID) {
      allG.push_back(g);
    }
  }
}

static void convertEpiOrderEff0(Rcpp::List ep,
				std::vector<epistasis>& Epistasis) {
  Rcpp::List element;
  // For epistasis, the numID must be sorted, but never with order effects.
  // Things come sorted (or not) from R.
  Epistasis.resize(ep.size());
  for(int i = 0; i != ep.size(); ++i) {
    element = ep[i];
    Epistasis[i].NumID = Rcpp::as<std::vector<int> >(element["NumID"]);
    Epistasis[i].names = Rcpp::as<std::vector<std::string> >(element["ids"]);
    Epistasis[i].s = as<double>(element["s"]);
  }
}

static void convertFitnessEffects(Rcpp::List rFE,
				  fitnessEffectsAll& fe) {

  // Yes, some of the things below are data.frames in R, but for
  // us that is used just as a list.
  Rcpp::List rrt = rFE["long.rt"];
  Rcpp::List re = rFE["long.epistasis"];
  Rcpp::List ro = rFE["long.orderEffects"];
  Rcpp::List rgi = rFE["long.geneNoInt"];
  Rcpp::List rgm = rFE["geneModule"];
  bool rone = rFE["gMOneToOne"];
  

  if(rrt.size()) {
    rTable_to_Poset(rrt, fe.Poset);
  } 
  if(re.size()) {
    convertEpiOrderEff0(re, fe.Epistasis);
    // convertEpiOrderEff(re, fe.Epistasis, fe.allEpistRTG);
  } 
  if(ro.size()) {
    convertEpiOrderEff(ro, fe.orderE, fe.allOrderG);
  } 
  if(rgi.size()) {
    convertNoInts(rgi, fe.genesNoInt);
  } else {
    fe.genesNoInt.shift = -9L;
  }
  rGM_GeneModule(rgm, fe.geneModules, fe.geneModulesLong);
  sort(fe.allOrderG.begin(), fe.allOrderG.end());
  // sort(fe.allEpistRTG.begin(), fe.allEpistRTG.end());
  fe.gMOneToOne = rone;
}


// It is simple to write specialized functions for when
// there are no restrictions or no order effects , etc.

static Genotype createNewGenotype(const Genotype& parent,
				  const std::vector<int>& mutations,
				  const fitnessEffectsAll& fe,
				  std::mt19937& ran_gen) {
  Genotype newGenot = parent;
  std::vector<int> tempOrder;
  bool sort_rest = false;
  bool sort_epist = false;
  for(auto g : mutations) {
    if( g >= fe.genesNoInt.shift ) {
      newGenot.rest.push_back(g);
      sort_rest = true;
    } else {
      if( binary_search(fe.allOrderG.begin(), fe.allOrderG.end(), g) ) {
	tempOrder.push_back(g);
      } else {
	newGenot.epistRtEff.push_back(g);
	sort_epist = true;
      }
    }
  }

  // With chromothripsis and order, we just randomly insert them
  if(tempOrder.size() > 1)
    shuffle(tempOrder.begin(), tempOrder.end(), ran_gen);
  for(auto g : tempOrder)
    newGenot.orderEff.push_back(g);

  if(sort_rest)
    sort(newGenot.rest.begin(), newGenot.rest.end());
  if(sort_epist)
    sort(newGenot.epistRtEff.begin(), newGenot.epistRtEff.end());
  
  return newGenot;
}


static void convertGenotype(Rcpp::List rGE,
			    Genotype& g) {
  Rcpp::IntegerVector oe = rGE["orderEffGenes"];
  Rcpp::IntegerVector ert = rGE["epistRTGenes"];
  Rcpp::IntegerVector rest = rGE["noInteractionGenes"];

  g.orderEff = Rcpp::as<std::vector<int> > (oe);
  g.epistRtEff = Rcpp::as<std::vector<int> > (ert);
  g.rest = Rcpp::as<std::vector<int> > (rest);
  sort(g.epistRtEff.begin(), g.epistRtEff.end());
  sort(g.rest.begin(), g.rest.end());
}



static void genotypeFitness(const Genotype& Ge,
			    const fitnessEffectsAll& fe,
			    std::vector<double>& s) {

  std::vector<int> r1 (Ge.epistRtEff);
  r1.insert( r1.end(), Ge.orderEff.begin(), Ge.orderEff.end());
  sort(r1.begin(), r1.end()); // yes, for rT checking

  if(fe.genesNoInt.shift > 0) {
    int shift = fe.genesNoInt.shift;
    for(auto  r : Ge.rest ) {
      s.push_back(fe.genesNoInt.s[r - shift]);
    }
  }

  // rT checking: with the r1 vector, call something like checkConstraints

  // finally, order checking
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
    if( (Poset[(*it_mutatedModule)].parentsNumID.size() == 1) &&
	(Poset[(*it_mutatedModule)].parentsNumID[0] == 0) ) { //Depends only on root.
      // FIXME: isn't it enough to check the second condition?
      s_vector.push_back(Poset[(*it_mutatedModule)].s);
    } else {
      sumDepsMet = 0;
      numDeps = Poset[(*it_mutatedModule)].parentsNumID.size();
      for(auto it_Parents = Poset[(*it_mutatedModule)].parentsNumID.begin();
	  it_Parents != Poset[(*it_mutatedModule)].parentsNumID.end();
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




void printPoset(const std::vector<geneDeps>& Poset) {

  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Restriction table (internal) *******" 
	      << std::endl;
  if(!Poset.size()) {
    Rcpp::Rcout << "No posets: restriction table of size 0"<< std::endl;
  } else {
    Rcpp::Rcout << "Size = " << (Poset.size() - 1) << std::endl;
    for(size_t i = 1; i != Poset.size(); ++i) {
      // We do not show the Poset[0]
      Rcpp::Rcout <<"\t Dependent Module or gene (child) " << i 
		  << ". childNumID: " << Poset[i].childNumID 
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
      for(auto c : Poset[i].parentsNumID)
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
}

static void printGeneToModule(const 
			      std::vector<geneToModuleLong>& geneModulesLong,
			      const bool gMOneToOne) {
  // Rcpp::Rcout << 
  //   "\n\n******** geneModule table (internal) *******:\nGene name\t Gene NumID\t Module name\t Module NumID\n";
  // for(auto it = geneModulesLong.begin(); it != geneModulesLong.end(); ++it) {
  //   Rcpp::Rcout << '\t' << it->GeneName << '\t' << it->GeneNumID << '\t' 
  // 		<< it->ModuleName << '\t' << it->ModuleNumID << std::endl;
  // }

  Rcpp::Rcout << 
    "\n\n******** geneModule table (internal) *******:\n" <<
    std::setw(14) << std::left << "Gene name" << std::setw(14) << "Gene NumID" << std::setw(14)
	      << "Module name" << std::setw(14) << "Module NumID" << "\n";
  for(auto it = geneModulesLong.begin(); it != geneModulesLong.end(); ++it) {
    Rcpp::Rcout << std::setw(14) << std::left << it->GeneName << std::setw(14)
		<< it->GeneNumID << std::setw(14) << it->ModuleName
		<< std::setw(14) << it->ModuleNumID << std::endl;
  }


  if(gMOneToOne)
    Rcpp::Rcout << "This is a dummy module table: each module is one gene."
		<< std::endl;
}


static void printOtherEpistasis(const std::vector<epistasis>& Epistasis,
				const std::string effectName,
				const std::string sepstr) {
  Rcpp::Rcout << "\n **********  General epistatic effects (internal) *******"
	      << std::endl;
  if(!Epistasis.size()) {
    Rcpp::Rcout << "No general " << effectName << std::endl;
  } else {
    Rcpp::Rcout << " Number of " << effectName <<"s = " << Epistasis.size();
    for(size_t i = 0; i != Epistasis.size(); ++i) {
      Rcpp::Rcout << "\n\t " << effectName << " " << i + 1 << ": " <<
	". Modules or Genes (names) = " << Epistasis[i].names[0];
      for(size_t j = 1; j != Epistasis[i].NumID.size(); ++j) {
	Rcpp::Rcout << sepstr << Epistasis[i].names[j] ;
      }
      Rcpp::Rcout << ".\t Modules or Genes (NumID) = " << Epistasis[i].NumID[0];
      for(size_t j = 1; j != Epistasis[i].NumID.size(); ++j) {
	Rcpp::Rcout << sepstr << Epistasis[i].NumID[j] ;
      }
      Rcpp::Rcout << ".\t s = " << Epistasis[i].s;
    }
  }
  Rcpp::Rcout << std::endl;
}

static void printNoInteractionGenes(const genesWithoutInt& genesNoInt) {
  Rcpp::Rcout << "\n **********  All remaining genes without interactions (internal) *******"
	      << std::endl;
  
  if(genesNoInt.shift <= 0) {
    Rcpp::Rcout << "No other genes without interactions" << std::endl;
  } else {
    Rcpp::Rcout << std::setw(14) << std::left << "Gene name" << std::setw(14)
		<< "Gene NumID" << std::setw(14) << "s" << std::endl;
    for(size_t i = 0; i != genesNoInt.NumID.size(); ++i) {
      Rcpp::Rcout << std::setw(14) << std::left << genesNoInt.names[i]
		  << std::setw(14) << genesNoInt.NumID[i]	
		  << std::setw(14) << genesNoInt.s[i] << '\n';
    }
  }
}

static void printFitnessEffects(const fitnessEffectsAll& fe) {
  printGeneToModule(fe.geneModulesLong, fe.gMOneToOne);
  printPoset(fe.Poset);
  printOtherEpistasis(fe.orderE, "order effect", " > ");
  printOtherEpistasis(fe.Epistasis, "epistatic interaction", ", ");
  printNoInteractionGenes(fe.genesNoInt);
}

// // [[Rcpp::export]]
// void wrap_test_rt(Rcpp::List rtR, Rcpp::DataFrame rGM) {
//   std::vector<geneDeps> Poset;
//   std::vector<geneToModule> geneModules;
//   std::vector<geneToModuleLong> geneModulesLong;

//   rGM_GeneModule(rGM, geneModules, geneModulesLong);
//   printGeneToModule(geneModulesLong);

//   rTable_to_Poset(rtR, Poset);
//   printPoset(Poset);
// }




// // stretching vector of mutatedDrv unlikely to speed up;
// // access (seeing if a module gene in there) is O(1) but I do
// // that for every gene and then sum over the vector (linear in
// // number of positions?)

// // Could be made much faster if we could assume lower numbered
// // positions cannot depend on higher numbered ones. Nope, not that much.








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


// static void printOrderEffects(const std::vector<epistasis>& orderE) {
//   // do something
// }



// [[Rcpp::export]]
void readFitnessEffects(Rcpp::List rFE,
			bool echo) {
  fitnessEffectsAll fitnessEffects;
  convertFitnessEffects(rFE, fitnessEffects);
  if(echo) {
     printFitnessEffects(fitnessEffects);
  }
}





// // [[Rcpp::export]]
// void wrap_FitnessEffects(Rcpp::List rtR,
// 			Rcpp::List longEpistasis,
// 			Rcpp::List longOrderE,
// 			Rcpp::DataFrame rGM,
// 			Rcpp::DataFrame geneNoInt,
// 			Rcpp::IntegerVector gMOneToOne) {

   
//   std::vector<geneDeps> Poset;
//   std::vector<geneToModule> geneModules;
//   std::vector<geneToModuleLong> geneModulesLong;

//   convertFitnessEffects(rtR, longEpistasis, longOrderE,
// 			rGM, geneNoInt, gMOneToOne,
// 			Poset, geneModules, geneModulesLong);
//   printGeneToModule(geneModulesLong);
//   printPoset(Poset);
// }



// void readFitnessEffects(Rcpp::List rtR,
// 			       Rcpp::List longEpistasis,
// 			       Rcpp::List longOrderE,
// 			       Rcpp::DataFrame rGM,
// 			       Rcpp::DataFrame geneNoInt,
// 			       Rcpp::IntegerVector gMOneToOne,
// 			       std::vector<geneDeps>&  Poset,
// 			       std::vector<geneToModule>& geneModules,
// 			       std::vector<geneToModuleLong>& geneModulesLong
// 			       // return objects
// 			       // something else
// 			       ) {
//    rTable_to_Poset(rtR, Poset);
//    rGM_GeneModule(rGM, geneModules, geneModulesLong);
 
// }


// epistasis and order can also be of modules, as we map modules to genes
// for that.



// How we store a given genotype and how we store the restrictions need to
// have a 1-to-1 mapping.


// checkEpistasis: for every epistatic set, check if there?

// checkOrderEffects: we find each member in the vector, but the indices
// must be correlative.






















