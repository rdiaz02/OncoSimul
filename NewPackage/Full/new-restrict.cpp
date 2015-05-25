#include <Rcpp.h>
#include <iomanip> // for setw
#include <algorithm>    
#include <random>

#define DP2(x) {Rcpp::Rcout << "\n DEBUG2: Value of " << #x << " = " << x << std::endl;}

using namespace Rcpp ;
using std::vector;
using std::back_inserter;

// FIXME: move this later
int seed = 1; 
std::mt19937 ran_gen(seed);

enum class Dependency {monotone, semimonotone, xmpn, NA}; // monotone, semimonotone

inline static Dependency stringToDep(const std::string& dep) {
  if(dep == "monotone")
    return Dependency::monotone;
  else if(dep == "semimonotone")
    return Dependency::semimonotone;
  else if(dep == "xmpn")
    return Dependency::xmpn;
  else 
    throw std::out_of_range("Not a valid typeDep");
}

inline static std::string depToString(const Dependency dep) {
  switch(dep) {
  case Dependency::monotone:
    return "CMPN or monotone";
  case Dependency::semimonotone:
    return "DMPN or semimonotone";
  case Dependency::xmpn:
    return "XMPN (XOR)";
  case Dependency::NA:
    return "NA";
  default:
    throw std::out_of_range("Not a valid dependency");
  }
}

struct Gene_Module_struct {
  std::string GeneName;
  std::string ModuleName;
  int GeneNumID;
  int ModuleNumID;
};

struct Poset_struct {
  Dependency typeDep;
  int childNumID; //Not redundant
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
  
  // We use allOrderG or allEpistRTG to place new mutations in their
  // correct place (orderEff or epistRtEff). Only one is needed.  Use the
  // one that is presumably always shorter which is allOrderG. And this is
  // sorted.
  std::vector<int> allOrderG; // Modules or genes if one-to-one.
  // std::vector<int> allEpistRTG;

  // This makes it faster to run evalPosetConstraints
  std::vector<int> allPosetG; //Modules or genes if one-to-one
  std::vector<Poset_struct> Poset;
  std::vector<epistasis> Epistasis;
  std::vector<epistasis> orderE;
  // std::vector<Gene_Module_struct> Gene_Module_tabl;
  std::vector<Gene_Module_struct> Gene_Module_tabl;
  genesWithoutInt genesNoInt;
};


// There are no shared genes in order and epist.  Any gene in orderEff can
// also be in the posets or general epistasis, but orderEff is only for
// those that have order effects.

// For all genes for which there are no order effects, any permutation of
// the same mutations is the same genotype, and has the same fitness. That
// is why we separate orderEff, which is strictly in the order in which
// mutations accumulate, and thus usorted, from the other effects, that
// are always kept sorted.

// rest are those genes that have no interactions. Evaluating their
// fitness is simple, and there can be no modules here.
struct Genotype {
  std::vector<int> orderEff;
  std::vector<int> epistRtEff; //always sorted
  std::vector<int> rest; // always sorted
};


// For users: if something depends on 0, that is it. No further deps.
// And do not touch the 0 in Gene_Module_table.

std::vector<Poset_struct> rTable_to_Poset(Rcpp::List rt) { 

  // The restriction table, or Poset, has a first element
  // with nothing, so that all references by mutated gene
  // are simply accessing the Poset[mutated gene] without
  // having to remember to add 1, etc.

  std::vector<Poset_struct> Poset;
  
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

    // if(i != Poset[i].childNumID) {
    // Nope: this assumes we only deal with posets.
    //   // Rcpp::Rcout << "\n childNumID, original = " << as<int>(rt_element["childNumID"]);
    //   // Rcpp::Rcout << "\n childNumID, Poset = " << Poset[i].childNumID;
    //   // Rcpp::Rcout << "\n i = " << i << std::endl;
    //   throw std::logic_error("childNumID != index");
    // }
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

    // Should not be needed if R always does what is should. Disable later?
    if(! is_sorted(Poset[i].parentsNumID.begin(), Poset[i].parentsNumID.end()) )
      throw std::logic_error("ParentsNumID not sorted");
    
    if(isinf(Poset[i].s))
      Rcpp::Rcout << "WARNING: at least one s is infinite" 
		  << std::endl;
    if(isinf(Poset[i].sh) && (Poset[i].sh > 0))
      Rcpp::Rcout << "WARNING: at least one sh is positive infinite" 
		  << std::endl;
  }
  return Poset;
}


std::vector<Gene_Module_struct> R_GeneModuleToGeneModule(Rcpp::List rGM) {

  std::vector<Gene_Module_struct> geneModule;

  Rcpp::IntegerVector GeneNumID = rGM["GeneNumID"];
  Rcpp::IntegerVector ModuleNumID = rGM["ModuleNumID"];
  Rcpp::CharacterVector GeneName = rGM["Gene"];
  Rcpp::CharacterVector ModuleName = rGM["Module"];
  // geneModule.resize(GeneNumID.size());
  geneModule.resize(GeneNumID.size()); // remove later

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != GeneNumID[i])
      throw std::logic_error(" i != GeneNumID");
    // geneModule[i].GeneNumID = GeneNumID[i];
    // geneModule[i].ModuleNumID = ModuleNumID[i];
    // remove these later?
    geneModule[i].GeneNumID = GeneNumID[i];
    geneModule[i].ModuleNumID = ModuleNumID[i];
    geneModule[i].GeneName = GeneName[i];
    geneModule[i].ModuleName = ModuleName[i];
  }
 
  return geneModule;
}


std::vector<int> GeneToModule(const std::vector<int>& Drv,
			     const 
			      std::vector<Gene_Module_struct>& Gene_Module_tabl,
			      const bool sortout, const bool uniqueout) {
  
  std::vector<int>  mutatedModules;
  
  for(auto it = Drv.begin(); it != Drv.end(); ++it) {
    mutatedModules.push_back(Gene_Module_tabl[(*it)].ModuleNumID);
  }
  // sortout and uniqueout returns a single element of each. uniqueout only removes
  // successive duplicates. sortout without unique is just useful for knowing
  // what happens for stats, etc. Neither sortout nor uniqueout for keeping
  // track of order of module events.
  if(sortout) {
    sort( mutatedModules.begin(), mutatedModules.end() );
  }
  if(uniqueout) {
    mutatedModules.erase( unique( mutatedModules.begin(), 
				  mutatedModules.end() ), 
			  mutatedModules.end() );
  }
  return mutatedModules;
}




genesWithoutInt convertNoInts(Rcpp::List nI) {

  genesWithoutInt genesNoInt;
  // FIXME: I think I want to force Gene in long.geneNoInt to be a char.vector
  // to avoid this transformation.
  Rcpp::CharacterVector names = nI["Gene"];
  // Rcpp::IntegerVector id = nI["GeneNumID"];
  // Rcpp::NumericVector s1 = nI["s"];
  
  genesNoInt.names = Rcpp::as<std::vector<std::string> >(names);
  genesNoInt.NumID = Rcpp::as<std::vector<int> >(nI["GeneNumID"]);
  genesNoInt.s = Rcpp::as<std::vector<double> >(nI["s"]);
  genesNoInt.shift = genesNoInt.NumID[0]; // FIXME: we assume mutations always
				      // indexed 1 to something. Not 0 to
				      // something.
  return genesNoInt;
}


std::vector<epistasis> convertEpiOrderEff(Rcpp::List ep) {

  std::vector<epistasis> Epistasis;
  
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
  return Epistasis;
}

std::vector<int> sortedAllOrder(std::vector<epistasis>& E) {
  
  std::vector<int> allG;
  for(auto ec : E) {
    for(auto g : ec.NumID) {
      allG.push_back(g);
    }
  }
  sort(allG.begin(), allG.end());
  allG.erase( unique( allG.begin(), allG.end()),
		      allG.end());
  return allG;
}

std::vector<int> sortedAllPoset(std::vector<Poset_struct>& Poset) {
  // Yes, this could be done inside rTable_to_Poset but this is cleaner
  // and will only add very little time. 
  std::vector<int> allG;
  for(auto p : Poset) {
    allG.push_back(p.childNumID);
  }
  sort(allG.begin(), allG.end());
  allG.erase( unique( allG.begin(), allG.end()),
		      allG.end());
  return allG; 
}

fitnessEffectsAll convertFitnessEffects(Rcpp::List rFE) {
  // Yes, some of the things below are data.frames in R, but for
  // us that is used just as a list.

  fitnessEffectsAll fe;
  
  Rcpp::List rrt = rFE["long.rt"];
  Rcpp::List re = rFE["long.epistasis"];
  Rcpp::List ro = rFE["long.orderEffects"];
  Rcpp::List rgi = rFE["long.geneNoInt"];
  Rcpp::List rgm = rFE["geneModule"];
  bool rone = rFE["gMOneToOne"];
  

  if(rrt.size()) {
    fe.Poset = rTable_to_Poset(rrt);
  } 
  if(re.size()) {
    fe.Epistasis = convertEpiOrderEff(re);
  } 
  if(ro.size()) {
    fe.orderE = convertEpiOrderEff(ro);
  } 
  if(rgi.size()) {
    fe.genesNoInt = convertNoInts(rgi);
  } else {
    fe.genesNoInt.shift = -9L;
  }
  fe.Gene_Module_tabl = R_GeneModuleToGeneModule(rgm);
  fe.allOrderG = sortedAllOrder(fe.orderE);
  fe.allPosetG = sortedAllPoset(fe.Poset);
  fe.gMOneToOne = rone;

  return fe;
}


// It is simple to write specialized functions for when
// there are no restrictions or no order effects , etc.
Genotype createNewGenotype(const Genotype& parent,
			   const std::vector<int>& mutations,
			   const fitnessEffectsAll& fe,
			   std::mt19937& ran_gen) {
  Genotype newGenot = parent;
  std::vector<int> tempOrder; // holder for multiple muts if order.
  bool sort_rest = false;
  bool sort_epist = false;

  // Order of ifs: I suspect order effects rare. No idea about
  // non-interaction genes, but if common the action is simple.
  for(auto g : mutations) {
    if( (fe.genesNoInt.shift < 0) || (g < fe.genesNoInt.shift) ) { // Gene with int
      // We can be dealing with modules
      int m; 
      if(fe.gMOneToOne) {
	m = g; 
      } else {
	m = fe.Gene_Module_tabl[g].ModuleNumID;
      }
      if( !binary_search(fe.allOrderG.begin(), fe.allOrderG.end(), m) ) {
	newGenot.epistRtEff.push_back(g);
	sort_epist = true;
      } else {
	tempOrder.push_back(g);
      }
    } else {
      // No interaction genes so no module stuff
      newGenot.rest.push_back(g);
      sort_rest = true;
    }
  }    

  // If there is order but multiple simultaneous mutations
  // (chromothripsis), we randomly insert them
  if(tempOrder.size() > 1)
    shuffle(tempOrder.begin(), tempOrder.end(), ran_gen);
  for(auto g : tempOrder)
    newGenot.orderEff.push_back(g);

  // Sorting done at end, in case multiple mutations
  if(sort_rest)
    sort(newGenot.rest.begin(), newGenot.rest.end());
  if(sort_epist)
    sort(newGenot.epistRtEff.begin(), newGenot.epistRtEff.end());
  
  return newGenot;
}

// FIXME: Prepare specialized functions: 
// Specialized functions:
// Never interactions: push into rest and sort. Identify by shift == 1.
// Never no interactions: remove the if. shift == -9.


Genotype convertGenotypeFromR(Rcpp::IntegerVector rG,
			      fitnessEffectsAll& fe) {
  // A genotype is of one kind or another depending on what genes are of
  // what type.

  std::vector<int> gg = Rcpp::as<std::vector<int> > (rG);
  Genotype newGenot;


  // Very similar to logic in createNewGenotype for placing each gene in
  // its correct place, which needs to look at module mapping.
  for(auto g : gg) {
    if( (fe.genesNoInt.shift < 0) || (g < fe.genesNoInt.shift) ) { // Gene with int
      // We can be dealing with modules
      int m; 
      if(fe.gMOneToOne) {
	m = g; 
      } else {
	m = fe.Gene_Module_tabl[g].ModuleNumID;
      }
      if( !binary_search(fe.allOrderG.begin(), fe.allOrderG.end(), m) ) {
	newGenot.epistRtEff.push_back(g);
      } else {
	newGenot.orderEff.push_back(g);
      }
    } else {
      // No interaction genes so no module stuff
      newGenot.rest.push_back(g);
    }
  }    

  sort(newGenot.rest.begin(), newGenot.rest.end());
  sort(newGenot.epistRtEff.begin(), newGenot.epistRtEff.end());
  
  return newGenot;
}


// Genotype convertGenotypeFromR(Rcpp::List rGE) {

//   Genotype g;
			
//   Rcpp::IntegerVector oe = rGE["orderEffGenes"];
//   Rcpp::IntegerVector ert = rGE["epistRTGenes"];
//   Rcpp::IntegerVector rest = rGE["noInteractionGenes"];

//   g.orderEff = Rcpp::as<std::vector<int> > (oe);
//   g.epistRtEff = Rcpp::as<std::vector<int> > (ert);
//   g.rest = Rcpp::as<std::vector<int> > (rest);
//   sort(g.epistRtEff.begin(), g.epistRtEff.end());
//   sort(g.rest.begin(), g.rest.end());

//   return g;
// }




bool match_order_effects(std::vector<int> O, std::vector<int> G) {
  //As the name sayes: we check if the order effect is matched
  if(G.size() < O.size()) return false;
  
  std::vector<int>::iterator p;
  std::vector<size_t> vdist;
  
  auto itb = G.begin();
  
  for(auto o : O) {
    p = find(G.begin(), G.end(), o);
    if( p == G.end() ) {
      return false;
    } else {
      vdist.push_back(std::distance( itb, p ));
    }
  }
  // Rcpp::Rcout << "     ";
  // for(auto vv : vdist ) {
  //   Rcpp::Rcout << vv << " ";
  // }
  // Rcpp::Rcout << std:: endl;
  if( is_sorted(vdist.begin(), vdist.end()) ) {
    return true;
  } else {
    return false;
  }
}

std::vector<double> evalOrderEffects(const std::vector<int>& mutatedM,
				     const std::vector<epistasis>& OE) {
  std::vector<double> s;
  for(auto o : OE) {
    if(match_order_effects(o.NumID, mutatedM))
      s.push_back(o.s);
  }
  return s;
}


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


std::vector<double> evalEpistasis(const std::vector<int>& mutatedModules,
				  const std::vector<epistasis>& Epistasis) {
  std::vector<double> s;

  // check_disable_later
  if(! is_sorted(mutatedModules.begin(), mutatedModules.end()))
    throw std::logic_error("mutatedModules not sorted in evalEpistasis");
  
  for(auto p : Epistasis ) {
    if(p.NumID[0] > 0 ) {
      if(includes(mutatedModules.begin(), mutatedModules.end(),
		   p.NumID.begin(), p.NumID.end()))
	s.push_back(p.s);
    } else {
      if(match_negative_epist(p.NumID, mutatedModules))
	s.push_back(p.s);
    }
  }
  // An alternative, but more confusing way
  // for(auto p : Epistasis ) {
  //   if(p.NumID[0] > 0 ) {
  //     if(!includes(mutatedModules.begin(), mutatedModules.end(),
  // 		   p.NumID.begin(), p.NumID.end()))
  // 	continue;
  //   } else {
  //     if(!match_negative_epist(p.NumID, mutatedModules))
  // 	continue;
  //   }
  //   s.push_back(p.s);
  // }
  return s;
}

// FIXME: can we make it faster if we know each module a single gene?
// FIXME: if genotype is always kept sorted, with drivers first, can it be
// faster? As well, note that number of drivers is automatically known
// from this table of constraints.

// int getPosetByChild(const int child,
// 		    const std::vector<Poset_struct>& Poset) {
  
// }

std::vector<double> evalPosetConstraints(const std::vector<int>& mutatedModules,
					 const std::vector<Poset_struct>& Poset,
					 const std::vector<int>& allPosetG) {

  // check_disable_later
  if(! is_sorted(mutatedModules.begin(), mutatedModules.end()))
    throw std::logic_error("mutatedModules not sorted in evalPosetConstraints");

  size_t numDeps;
  size_t sumDepsMet = 0;
  Dependency deptype;
  std::vector<int> parent_matches;
  std::vector<double> s;

  //This works reverted w.r.t. to evalOrderEffects and evalEpistasis:
  //there, I examine if the effect is present in the genotype. Here, I
  //examine if the genotype satisfies the constraints.


  // Since the genotype can contain genes not in the poset, first find
  // those that are mutated in the genotype AND are in the poset. Then
  // check if the mutated have restrictions satisfied.
  
  std::vector<int> MPintersect(allPosetG.size());
  std::set_intersection(allPosetG.begin(), allPosetG.end(),
			mutatedModules.begin(), mutatedModules.end(),
			std::back_inserter(MPintersect));
  
  // We know MPintersect is sorted, so we can avoid an O(n*n) loop
  size_t i = 0;
  for(auto m : MPintersect) {
    while ( Poset[i].childNumID != m) ++i;
    // Not to catch the twisted case of an XOR with a 0 and something else
    // as parents
    if( (Poset[i].parentsNumID[0] == 0) &&
	(Poset[i].parentsNumID.size() == 1)) {
      s.push_back(Poset[i].s);
    } else {
      parent_matches.clear();
      std::set_intersection(mutatedModules.begin(), mutatedModules.end(),
			    Poset[i].parentsNumID.begin(),
			    Poset[i].parentsNumID.end(),
			    back_inserter(parent_matches));
      sumDepsMet = parent_matches.size();
      numDeps = Poset[i].parentsNumID.size();
      deptype = Poset[i].typeDep;

      if( ((deptype == Dependency::semimonotone) &&
	   (sumDepsMet)) ||
	  ((deptype == Dependency::monotone) &&
	   (sumDepsMet == numDeps)) ||
	  ((deptype == Dependency::xmpn) &&
	   (sumDepsMet == 1)) ) {
	s.push_back(Poset[i].s);
      } else {
	s.push_back(Poset[i].sh);
      }
    }
  }
      // for(auto parent : Poset[i].parentsNumID ) {
      // 	parent_module_mutated = binary_search(mutatedModules.begin(), 
      // 					      mutatedModules.end(),
      // 					      parent);
      // 	if(parent_module_mutated) {
      // 	  ++sumDepsMet;
      // 	  if( deptype == Dependency::semimonotone ) {
      // 	    known = true;
      // 	    break;
      // 	  }
      // 	  if( (deptype == Dependency::xmpn && (sumDepsMet > 1))) {
      // 	    knwon = true;
      // 	    break;
      // 	  }
      // 	}
      // }
  return s;
}
  


std::vector<double> evalGenotypeFitness(const Genotype& ge,
					const fitnessEffectsAll& F){
  std::vector<double> s;

  // Genes without any restriction or epistasis are just genes. No modules.
  // So simple we do it here.
  if(F.genesNoInt.shift > 0) {
    int shift = F.genesNoInt.shift;
    for(auto  r : ge.rest ) {
      s.push_back(F.genesNoInt.s[r - shift]);
    }
  }

  // For the rest, there might be modules. Three different effects on
  // fitness possible: as encoded in Poset, general epistasis, order effects.
  
  // Epistatis and poset are checked against all mutations. Create single
  // sorted vector with all mutations and map to modules, if needed. Then
  // eval.
  std::vector<int> mutG (ge.epistRtEff);
  mutG.insert( mutG.end(), ge.orderEff.begin(), ge.orderEff.end());
  std::vector<int> mutatedModules;
  if(F.gMOneToOne) {
    sort(mutG.begin(), mutG.end()); 
    mutatedModules = mutG;
  } else {
    mutatedModules = GeneToModule(mutG, F.Gene_Module_tabl, true, true);
  }
  std::vector<double> srt =
    evalPosetConstraints(mutatedModules, F.Poset, F.allPosetG);
  std::vector<double> se =
    evalEpistasis(mutatedModules, F.Epistasis);
  
  // For order effects we need a new vector of mutatedModules:
  if(F.gMOneToOne) {
    mutatedModules = ge.orderEff;
  } else {
    mutatedModules = GeneToModule(ge.orderEff, F.Gene_Module_tabl, false, true);
  }
  std::vector<double> so =
    evalOrderEffects(mutatedModules, F.orderE);

  // I keep s, srt, se, so separate for now for debugging.
  s.insert(s.end(), srt.begin(), srt.end());
  s.insert(s.end(), se.begin(), se.end());
  s.insert(s.end(), so.begin(), so.end());
  
  return s;
}

// No logs because of log(<=0)
inline double prodFitness(vector<double> s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    [](double x, double y) {return (x * (1 + y));});
}

// inline double logSumFitness(vector<double> s) {
//   return accumulate(s.begin(), s.end(), 1.0,
// 		    [](double x, double y) {return (x * (1 + y));})
// }



void printPoset(const std::vector<Poset_struct>& Poset) {

  int counterInfs = 0;
  int counterNegInfs = 0;
  Rcpp::Rcout << "\n **********  Poset or Restriction table (internal) *******" 
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
	Rcpp::Rcout << c << "; ";
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

void printGene_Module_table(const 
		       std::vector<Gene_Module_struct>& Gene_Module_tabl,
		       const bool gMOneToOne) {
  // Rcpp::Rcout << 
  //   "\n\n******** geneModule table (internal) *******:\nGene name\t Gene NumID\t Module name\t Module NumID\n";
  // for(auto it = Gene_Module_tabl.begin(); it != Gene_Module_tabl.end(); ++it) {
  //   Rcpp::Rcout << '\t' << it->GeneName << '\t' << it->GeneNumID << '\t' 
  // 		<< it->ModuleName << '\t' << it->ModuleNumID << std::endl;
  // }

  Rcpp::Rcout << 
    "\n\n******** geneModule table (internal) *******:\n" <<
    std::setw(14) << std::left << "Gene name" << std::setw(14) << "Gene NumID" << std::setw(14)
	      << "Module name" << std::setw(14) << "Module NumID" << "\n";
  for(auto it = Gene_Module_tabl.begin(); it != Gene_Module_tabl.end(); ++it) {
    Rcpp::Rcout << std::setw(14) << std::left << it->GeneName << std::setw(14)
		<< it->GeneNumID << std::setw(14) << it->ModuleName
		<< std::setw(14) << it->ModuleNumID << std::endl;
  }


  if(gMOneToOne)
    Rcpp::Rcout << "This is a dummy module table: each module is one gene."
		<< std::endl;
}


void printOtherEpistasis(const std::vector<epistasis>& Epistasis,
				const std::string effectName,
				const std::string sepstr) {
  Rcpp::Rcout << "\n **********  General " << effectName << "s (internal) *******"
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

void printNoInteractionGenes(const genesWithoutInt& genesNoInt) {
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

void printAllOrderG(const std::vector<int> ge) {
  Rcpp::Rcout << "\n **********  NumID of genes/modules in the order restrict. (internal) *******"
	      << std::endl;
  for(auto g : ge)
    Rcpp::Rcout << g << " ";
  Rcpp::Rcout << std::endl;
}


void printFitnessEffects(const fitnessEffectsAll& fe) {
  printGene_Module_table(fe.Gene_Module_tabl, fe.gMOneToOne);
  printPoset(fe.Poset);
  printOtherEpistasis(fe.orderE, "order effect", " > ");
  printOtherEpistasis(fe.Epistasis, "epistatic interaction", ", ");
  printNoInteractionGenes(fe.genesNoInt);
  printAllOrderG(fe.allOrderG);
}


// [[Rcpp::export]]
void readFitnessEffects(Rcpp::List rFE,
			bool echo) {
  // fitnessEffectsAll fitnessEffects;
  // convertFitnessEffects(rFE, fitnessEffects);
  fitnessEffectsAll fitnessEffects = convertFitnessEffects(rFE);
  if(echo) {
     printFitnessEffects(fitnessEffects);
  }
}



// [[Rcpp::export]]
double evalRGenotype(Rcpp::IntegerVector rG, Rcpp::List rFE, bool verbose) {
  const Rcpp::List rF(rFE);
  fitnessEffectsAll F = convertFitnessEffects(rF);
  Genotype g = convertGenotypeFromR(rG, F);
  vector<double> s = evalGenotypeFitness(g, F);
  if(verbose) {
    Rcpp::Rcout << "\n Individual s terms are :";
    for(auto i : s) Rcpp::Rcout << " " << i;
    Rcpp::Rcout << std::endl;
  }
  return prodFitness(s);
}


// [[Rcpp::export]]
void dummy(Rcpp::List rt) { 

  // The restriction table, or Poset, has a first element
  // with nothing, so that all references by mutated gene
  // are simply accessing the Poset[mutated gene] without
  // having to remember to add 1, etc.

  std::vector<Poset_struct> Poset;
  
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

    // if(i != Poset[i].childNumID) {
    // Nope: this assumes we only deal with posets.
    //   // Rcpp::Rcout << "\n childNumID, original = " << as<int>(rt_element["childNumID"]);
    //   // Rcpp::Rcout << "\n childNumID, Poset = " << Poset[i].childNumID;
    //   // Rcpp::Rcout << "\n i = " << i << std::endl;
    //   throw std::logic_error("childNumID != index");
    // }
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

    // Should not be needed if R always does what is should. Disable later?
    if(! is_sorted(Poset[i].parentsNumID.begin(), Poset[i].parentsNumID.end()) )
      throw std::logic_error("ParentsNumID not sorted");
    
    if(isinf(Poset[i].s))
      Rcpp::Rcout << "WARNING: at least one s is infinite" 
		  << std::endl;
    if(isinf(Poset[i].sh) && (Poset[i].sh > 0))
      Rcpp::Rcout << "WARNING: at least one sh is positive infinite" 
		  << std::endl;
  }
}

// [[Rcpp::export]]
void dummy2(Rcpp::List rFE) {
  // Yes, some of the things below are data.frames in R, but for
  // us that is used just as a list.

  fitnessEffectsAll fe;
  
  Rcpp::List rrt = rFE["long.rt"];
  Rcpp::List re = rFE["long.epistasis"];
  Rcpp::List ro = rFE["long.orderEffects"];
  Rcpp::List rgi = rFE["long.geneNoInt"];
  Rcpp::List rgm = rFE["geneModule"];
  bool rone = rFE["gMOneToOne"];
  

  if(rrt.size()) {
    fe.Poset = rTable_to_Poset(rrt);
  } 
  if(re.size()) {
    fe.Epistasis = convertEpiOrderEff(re);
  } 
  if(ro.size()) {
    fe.orderE = convertEpiOrderEff(ro);
  } 
  // if(rgi.size()) {
  //   fe.genesNoInt = convertNoInts(rgi);
  // } else {
  //   fe.genesNoInt.shift = -9L;
  // }
  
  // fe.Gene_Module_tabl = R_GeneModuleToGeneModule(rgm);
  // fe.allOrderG = sortedAllOrder(fe.orderE);
  // fe.allPosetG = sortedAllPoset(fe.Poset);
  fe.gMOneToOne = rone;
}


// // [[Rcpp::export]]
// void wrap_test_rt(Rcpp::List rtR, Rcpp::DataFrame rGM) {
//   std::vector<Poset_struct> Poset;
//   std::vector<geneToModule> geneModules;
//   std::vector<geneToModuleLong> geneModules;

//   R_GeneModuleToGeneModule(rGM, geneModules, geneModules);
//   printGeneToModule(geneModules);

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


// // [[Rcpp::export]]
// SEXP eval_Genotype(Rcpp::List rtR, 
// 		   Rcpp::DataFrame rGM,
// 		   Rcpp::IntegerVector genotype) {

  
//   // std::vector<geneToModule> geneModules;
//   std::vector<geneToModuleLong> geneModules;

//   std::vector<int> Genotype;
//   std::vector<double> s_vector;
//   std::vector<double> sh_vector;

//   std::vector<Poset_struct> Poset = rTable_to_Poset(rtR);
//   geneModules = R_GeneModuleToGeneModule(rGM, geneModules, geneModules);


//   Genotype = Rcpp::as< std::vector<int> >(genotype);

//   evalPosetConstraints(Genotype, Poset, geneModules,  s_vector, sh_vector);

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


// static void printOrderEffects(const std::vector<epistasis>& orderE) {
//   // do something
// }








// // [[Rcpp::export]]
// void wrap_FitnessEffects(Rcpp::List rtR,
// 			Rcpp::List longEpistasis,
// 			Rcpp::List longOrderE,
// 			Rcpp::DataFrame rGM,
// 			Rcpp::DataFrame geneNoInt,
// 			Rcpp::IntegerVector gMOneToOne) {

   
//   std::vector<Poset_struct> Poset;
//   std::vector<geneToModule> geneModules;
//   std::vector<geneToModuleLong> geneModules;

//   convertFitnessEffects(rtR, longEpistasis, longOrderE,
// 			rGM, geneNoInt, gMOneToOne,
// 			Poset, geneModules, geneModules);
//   printGeneToModule(geneModules);
//   printPoset(Poset);
// }



// void readFitnessEffects(Rcpp::List rtR,
// 			       Rcpp::List longEpistasis,
// 			       Rcpp::List longOrderE,
// 			       Rcpp::DataFrame rGM,
// 			       Rcpp::DataFrame geneNoInt,
// 			       Rcpp::IntegerVector gMOneToOne,
// 			       std::vector<Poset_struct>&  Poset,
// 			       std::vector<geneToModule>& geneModules,
// 			       std::vector<geneToModuleLong>& geneModules
// 			       // return objects
// 			       // something else
// 			       ) {
//    rTable_to_Poset(rtR, Poset);
//    R_GeneModuleToGeneModule(rGM, geneModules, geneModules);
 
// }


// epistasis and order can also be of modules, as we map modules to genes
// for that.



// How we store a given genotype and how we store the restrictions need to
// have a 1-to-1 mapping.


// checkEpistasis: for every epistatic set, check if there?

// checkOrderEffects: we find each member in the vector, but the indices
// must be correlative.



















// static void evalPosetConstraints(const std::vector<int>& Drv,
// 			     const std::vector<Poset_struct>& Poset,
// 			     const std::vector<geneToModule>& geneToModules,
// 			     std::vector<double>& s_vector,
// 			     std::vector<double>& sh_vector) {
//   size_t numDeps;
//   size_t sumDepsMet = 0;
//   int parent_module_mutated = 0;
//   std::vector<int> mutatedModules;
  

//   // Map mutated genes to modules (mutatedModules) and then examine, for
//   // each mutatedModule, if its dependencies (in terms of modules) are met.

//   mutatedModules = DrvToModule(Drv, geneToModules);

//   for(auto it_mutatedModule = mutatedModules.begin();
//       it_mutatedModule != mutatedModules.end(); ++it_mutatedModule) {
//     if( (Poset[(*it_mutatedModule)].parentsNumID.size() == 1) &&
// 	(Poset[(*it_mutatedModule)].parentsNumID[0] == 0) ) { //Depends only on root.
//       // FIXME: isn't it enough to check the second condition?
//       s_vector.push_back(Poset[(*it_mutatedModule)].s);
//     } else {
//       sumDepsMet = 0;
//       numDeps = Poset[(*it_mutatedModule)].parentsNumID.size();
//       for(auto it_Parents = Poset[(*it_mutatedModule)].parentsNumID.begin();
// 	  it_Parents != Poset[(*it_mutatedModule)].parentsNumID.end();
// 	  ++it_Parents) {
// 	// if sorted, could use binary search
// 	// FIXME: try a set or sort mutatedModules?

// 	//parent_module_mutated =
// 	//std::binary_search(mutatedModules.begin(), mutatedModules.end(),
// 	//(*it_Parents))

// 	parent_module_mutated = 
// 	  (std::find(mutatedModules.begin(), 
// 		     mutatedModules.end(), 
// 		     (*it_Parents)) != mutatedModules.end());
// 	if(parent_module_mutated)  {
// 	  ++sumDepsMet;
// 	  if( Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone)
// 	    break;
// 	}
//       }
//       if( ((Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone) && 
// 	   (sumDepsMet)) ||
// 	  ((Poset[(*it_mutatedModule)].typeDep == Dependency::monotone) && 
// 	   (sumDepsMet == numDeps)) ) {
// 	s_vector.push_back(Poset[(*it_mutatedModule)].s);
//       } else {
// 	sh_vector.push_back(Poset[(*it_mutatedModule)].sh);
//       }
//     }
//   }
  
// }



// std::vector<double> evalPosetConstraints0(const std::vector<int>& mutatedModules,
// 					 const std::vector<Poset_struct>& Poset) {

//   size_t numDeps;
//   size_t sumDepsMet = 0;
//   int parent_module_mutated = 0;
//   // std::vector<int> mutatedModules;
  
//   std::vector<double> s;
  
//   // Map mutated genes to modules (mutatedModules) and then examine, for
//   // each mutatedModule, if its dependencies (in terms of modules) are met.

  

//   for(auto it_mutatedModule = mutatedModules.begin();
//       it_mutatedModule != mutatedModules.end(); ++it_mutatedModule) {
//     if( (Poset[(*it_mutatedModule)].parentsNumID.size() == 1) &&
// 	(Poset[(*it_mutatedModule)].parentsNumID[0] == 0) ) { //Depends only on root.
//       // FIXME: isn't it enough to check the second condition?
//       s.push_back(Poset[(*it_mutatedModule)].s);
//     } else {
//       sumDepsMet = 0;
//       numDeps = Poset[(*it_mutatedModule)].parentsNumID.size();
//       for(auto it_Parents = Poset[(*it_mutatedModule)].parentsNumID.begin();
// 	  it_Parents != Poset[(*it_mutatedModule)].parentsNumID.end();
// 	  ++it_Parents) {

// 	parent_module_mutated = binary_search(mutatedModules.begin(), 
// 					      mutatedModules.end(),
// 					      (*it_Parents));
// 	  // (std::find(mutatedModules.begin(), 
// 	  // 	     mutatedModules.end(), 
// 	  // 	     (*it_Parents)) != mutatedModules.end());
// 	if(parent_module_mutated)  {
// 	  ++sumDepsMet;
// 	  if( Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone)
// 	    break;
// 	}
//       }
//       if( ((Poset[(*it_mutatedModule)].typeDep == Dependency::semimonotone) && 
// 	   (sumDepsMet)) ||
// 	  ((Poset[(*it_mutatedModule)].typeDep == Dependency::monotone) && 
// 	   (sumDepsMet == numDeps)) ) {
// 	s.push_back(Poset[(*it_mutatedModule)].s);
//       } else {
// 	s.push_back(Poset[(*it_mutatedModule)].sh);
//       }
//     }
//   }
//   return s;
// }

