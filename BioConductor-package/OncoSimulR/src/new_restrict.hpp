#ifndef _NEW_RESTRICT_H__
#define _NEW_RESTRICT_H__

#include <Rcpp.h>
#include"debug_common.hpp"
#include <limits>

enum class Dependency {monotone, semimonotone, xmpn, single, NA}; 

inline Dependency stringToDep(const std::string& dep) {
  if(dep == "monotone") // AND, CBN, CMPN
    return Dependency::monotone;
  else if(dep == "semimonotone") // OR, SMN, DMPN
    return Dependency::semimonotone;
  else if(dep == "xmpn") // XOR, XMPN
    return Dependency::xmpn;
  else if(dep == "--") // for root, for example
    return Dependency::single;
  else 
    throw std::out_of_range("Not a valid typeDep");
  // We never create the NA from entry data. NA is reserved for Root.
}


inline std::string depToString(const Dependency dep) {
  switch(dep) {
  case Dependency::monotone:
    return "CMPN or monotone";
  case Dependency::semimonotone:
    return "DMPN or semimonotone";
  case Dependency::xmpn:
    return "XMPN (XOR)";
  case Dependency::single:
    return "--";
  default:
    throw std::out_of_range("Not a valid dependency");
  }
}


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
  // If shift is -9, no elements The next first two are not really
  // needed. Will remove later. Nope! we use them to provide nice output.
  std::vector<int> NumID;
  std::vector<std::string> names;
  std::vector<double> s;
};

struct Gene_Module_struct {
  std::string GeneName;
  std::string ModuleName;
  int GeneNumID;
  int ModuleNumID;
};

struct fitnessEffectsAll {
  bool gMOneToOne;
  int genomeSize; 
  // We use allOrderG or allEpistRTG to place new mutations in their
  // correct place (orderEff or epistRtEff). Only one is needed.  Use the
  // one that is presumably always shorter which is allOrderG. And this is
  // sorted.
  std::vector<int> allOrderG; // Modules or genes if one-to-one.
  // std::vector<int> allEpistRTG;

  // This makes it faster to run evalPosetConstraints
  std::vector<int> allPosetG; //Modules or genes if one-to-one. Only
			      //poset. Not epist.
  std::vector<Poset_struct> Poset;
  std::vector<epistasis> Epistasis;
  std::vector<epistasis> orderE;
  // std::vector<Gene_Module_struct> Gene_Module_tabl;
  std::vector<Gene_Module_struct> Gene_Module_tabl;
  std::vector<int> allGenes; //used whenever a mutation created. Genes,
			     //not modules. Sorted.
  std::vector<int> drv; // Sorted.
  genesWithoutInt genesNoInt;
  
};

struct fitness_as_genes {
  // fitnessEffectsAll in terms of genes.  Useful for output
  // conversions. There could be genes that are both in orderG and
  // posetEpistG. In such a case, only in orderG.
  // We only use a small part for now.
  // All are ordered vectors.
  std::vector<int> orderG;
  std::vector<int> posetEpistG;
  std::vector<int> noInt;
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

inline Genotype wtGenotype() {
  // Not needed but to make it explicit
  Genotype g;
  g.orderEff.resize(0);
  g.epistRtEff.resize(0);
  g.rest.resize(0);
  return g;
}
std::vector<int> genotypeSingleVector(const Genotype& ge);

inline bool operator==(const Genotype& lhs, const Genotype& rhs) {
  return (lhs.orderEff == rhs.orderEff) &&
    (lhs.epistRtEff == rhs.epistRtEff) &&
    (lhs.rest == rhs.rest);
}

inline bool operator<(const Genotype& lhs, const Genotype& rhs) {
  std::vector<int> lh = genotypeSingleVector(lhs);
  std::vector<int> rh = genotypeSingleVector(rhs);
  if( lh.size() < rh.size() ) return true;
  else if ( lh.size() > rh.size() ) return false;
  else {
    for(size_t i = 0; i != lh.size(); ++i) {
      if( lh[i] < rh[i] ) return true;
    }
    return false;
  }
}



inline double prodFitness(std::vector<double> s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    [](double x, double y) {return (x * std::max(0.0, (1 + y)));});
}



// I am looping twice
// inline double prodDeathFitness(vector<double> s) {
//   // For Bozic's
//   if( *min_element( s.begin(), s.end()) <= 99.0 ) {
//     return -99.0;
//   } else {
//     return accumulate(s.begin(), s.end(), 1.0,
// 		      [](double x, double y) {return (x * (1 - y));});
//   }
// }

// But using infinity deals with this. See below
// inline double prodDeathFitness(std::vector<double> s) {
//   double f = 1.0;
//   for(auto si : s) {
//     if( si <= -90.0 ) {
//       return 99.0;
//     } else {
//       f *= (1 - si);
//     }
//   }
//   return f;
// }


inline double prodDeathFitness(std::vector<double> s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    [](double x, double y) {return (x * std::max(0.0, (1 - y)));});
}



// inline double logSumFitness(vector<double> s) {
//   return accumulate(s.begin(), s.end(), 1.0,
// 		    [](double x, double y) {return (x * (1 + y));})
// }


void obtainMutations(const Genotype& parent,
		     const fitnessEffectsAll& fe,
		     int& numMutablePosParent,
		     std::vector<int>& newMutations,
		     std::mt19937& ran_gen);

Genotype createNewGenotype(const Genotype& parent,
			   const std::vector<int>& mutations,
			   const fitnessEffectsAll& fe,
			   std::mt19937& ran_gen);

std::vector<double> evalGenotypeFitness(const Genotype& ge,
					const fitnessEffectsAll& F);


fitnessEffectsAll convertFitnessEffects(Rcpp::List rFE);
std::vector<int> getGenotypeDrivers(const Genotype& ge, const std::vector<int>& drv);
void print_Genotype(const Genotype& ge);

fitness_as_genes fitnessAsGenes(const fitnessEffectsAll& fe);

std::map<int, std::string> mapGenesIntToNames(const fitnessEffectsAll& fe);

std::vector<int> getGenotypeDrivers(const Genotype& ge, const std::vector<int>& drv);
#endif

