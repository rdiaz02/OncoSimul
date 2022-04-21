//     Copyright 2013-2021 Ramon Diaz-Uriarte

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


// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include "debug_common.h"
#include "common_classes.h"
#include "new_restrict.h"
#include "exprtk.h"
#include <Rcpp.h>
#include <iomanip>
#include <algorithm>
#include <random>
#include <string>
#include <sstream>
#include <limits>
#include <regex>


using namespace Rcpp;
using std::vector;
using std::back_inserter;

std::string concatIntsString(const std::vector<int>& ints,
			     const std::string sep = ", ") {
  std::string strout;
  std::string comma = "";
  for(auto const &g : ints) {
    strout += (comma + std::to_string(g));
    comma = sep;
  }
  return strout;
}

// Cumulative product of (1 + s).
//      If using fitness landscape, a single s is passed, and that already
//      has a "-1" subtracted. So when a single s is passed, we just
//      return s + 1, which is the fitness (birth rate) as specified in
//      the fitness landscape. This applies to FDF too.
//      See function evalGenotypeFitness.
//      Whenever we call prodFitness (e.g., nr_fitness or
//      initMutantInitialization) it is called on the output of
//      evalGenotypeFitness
double prodFitness(const std::vector<double>& s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    [](double x, double y) {return (x * std::max(0.0, (1 + y)));});
}


// Plays same role as prodFitness, but for Bozic models
double prodDeathFitness(const std::vector<double>& s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    [](double x, double y) {return (x * std::max(0.0, (1 - y)));});
}


// To compute cumulative product of mutator effects
double prodMuts(const std::vector<double>& s) {
  return accumulate(s.begin(), s.end(), 1.0,
		    std::multiplies<double>());
}

// Compute cPDetect for Detection/stopping mechanism
double set_cPDetect(const double n2, const double p2,
		    const double PDBaseline) {
  return ( -log(1.0 - p2) *  (PDBaseline / (n2 - PDBaseline)) );
}

// Probability of detection given size
double probDetectSize(const double n, const double cPDetect,
		      const double PDBaseline) {
  if(n <= PDBaseline) {
    return 0;
  } else {
    return (1 - exp( -cPDetect * ( (n - PDBaseline)/PDBaseline ) ));
  }
}


// return true if this is detected (as given by the probability of detection as a function of size)
bool detectedSizeP(const double n, const double cPDetect,
		   const double PDBaseline, std::mt19937& ran_gen) {
  if(cPDetect < 0) {
    // As we OR, return false if this condition does not apply
    return false;
  } else {
    std::uniform_real_distribution<double> runif(0.0, 1.0);
    double prob = probDetectSize(n, cPDetect, PDBaseline);
    if(prob <= 0.0) return false;
    if(runif(ran_gen) <= prob) {
      return true;
    } else {
      return false;
    }
  }
}



bool operator==(const Genotype& lhs, const Genotype& rhs) {
  return (lhs.orderEff == rhs.orderEff) &&
    (lhs.epistRtEff == rhs.epistRtEff) &&
    (lhs.rest == rhs.rest) &&
    (lhs.flGenes == rhs.flGenes);
}

// Added for completeness, but not used now
// bool operator<(const Genotype& lhs, const Genotype& rhs) {
//   std::vector<int> lh = genotypeSingleVector(lhs);
//   std::vector<int> rh = genotypeSingleVector(rhs);
//   if( lh.size() < rh.size() ) return true;
//   else if ( lh.size() > rh.size() ) return false;
//   else {
//     for(size_t i = 0; i != lh.size(); ++i) {
//       if( lh[i] < rh[i] ) return true;
//     }
//     return false;
//   }
// }


TypeModel stringToModel(const std::string& mod) {
  if(mod == "exp")
    return TypeModel::exp;
  else if(mod == "bozic1")
    return TypeModel::bozic1;
  else if(mod == "mcfarlandlog")
    return TypeModel::mcfarlandlog;
  else if(mod == "mcfarlandlogd")
    return TypeModel::mcfarlandlog_d;
  else if (mod == "arbitrary")
    return TypeModel::arbitrary;
  else if (mod == "constant")
    return TypeModel::constant;
  else
    throw std::out_of_range("Not a valid TypeModel");
}


Dependency stringToDep(const std::string& dep) {
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



// genotype struct -> genotype as a vector of ints
vector<int> genotypeSingleVector(const Genotype& ge) {
  // orderEff in the order they occur. All others are sorted.
  std::vector<int> allgG;
  allgG.insert(allgG.end(), ge.orderEff.begin(), ge.orderEff.end());
  allgG.insert(allgG.end(), ge.epistRtEff.begin(), ge.epistRtEff.end());
  allgG.insert(allgG.end(), ge.rest.begin(), ge.rest.end());
  allgG.insert(allgG.end(), ge.flGenes.begin(), ge.flGenes.end());
  // this should not be unique'd as it aint' sorted
  return allgG;
}


vector<int> allGenesinFitness(const fitnessEffectsAll& F) {
  // Sorted
  std::vector<int> g0;

  if(F.Gene_Module_tabl.size()) {
    if( F.Gene_Module_tabl[0].GeneNumID != 0 )
      throw std::logic_error("\n Gene module table's first element must be 0."
			     " This should have been caught in R.");
    for(decltype(F.Gene_Module_tabl.size()) i = 1;
  	i != F.Gene_Module_tabl.size(); i++) {
      g0.push_back(F.Gene_Module_tabl[i].GeneNumID);
    }
  }
  for(auto const &b: F.genesNoInt.NumID) {
    g0.push_back(b);
  }
  for(auto const &b: F.fitnessLandscape.NumID) {
    g0.push_back(b);
  }
  sort(g0.begin(), g0.end());

  // Can we assume the fitness IDs go from 0 to n? Nope: because of
  // muEF. But we assume in several places that there are no repeated
  // elements in the output from this function.

  // We verify there are no repeated elements. That is to strongly
  // check our assumptions are right. Alternatively, return the "uniqued"
  // vector and do not check anything.
  std::vector<int> g0_cp(g0);
  g0.erase( unique( g0.begin(), g0.end() ), g0.end() );
  if(g0.size() != g0_cp.size())
    throw std::logic_error("\n allGenesinFitness: repeated genes. "
			   " This should have been caught in R.");
  return g0;
}

vector<int> allGenesinGenotype(const Genotype& ge){
  // Like genotypeSingleVector, but sorted
  std::vector<int> allgG;
  for(auto const &g1 : ge.orderEff)
    allgG.push_back(g1);
  for(auto const &g2 : ge.epistRtEff)
    allgG.push_back(g2);
  for(auto const &g3 : ge.rest)
    allgG.push_back(g3);
  for(auto const &g4 : ge.flGenes)
    allgG.push_back(g4);

  sort(allgG.begin(), allgG.end());
  // Remove duplicates see speed comparisons here:
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  // We assume here there are no duplicates. Yes, a gene can have both
  // fitness effects and order effects and be in the DAG. But it will be
  // in only one of the buckets.
  std::vector<int> g0_cp(allgG);
  allgG.erase( unique( allgG.begin(), allgG.end() ), allgG.end() );
  if(allgG.size() != g0_cp.size())
    throw std::logic_error("\n allGenesinGenotype: repeated genes."
			   " This should have been caught in R.");
  return allgG;
}


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

    parentsid = as<Rcpp::IntegerVector>(rt_element["parentsNumID"]);
    parents = as<Rcpp::CharacterVector>(rt_element["parents"]);

    if( parentsid.size() != parents.size() ) {
      throw std::logic_error("parents size != parentsNumID size. Bug in R code.");
    }

    for(int j = 0; j != parentsid.size(); ++j) {
      Poset[i].parentsNumID.push_back(parentsid[j]);
      Poset[i].parents.push_back( (Rcpp::as< std::string >(parents[j])) );
    }

    // Should not be needed if R always does what is should. Disable later?
    if(! is_sorted(Poset[i].parentsNumID.begin(), Poset[i].parentsNumID.end()) )
      throw std::logic_error("ParentsNumID not sorted. Bug in R code.");

    if(std::isinf(Poset[i].s))
      Rcpp::Rcout << "WARNING: at least one s is infinite"
		  << std::endl;
    if(std::isinf(Poset[i].sh) && (Poset[i].sh > 0))
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

  geneModule.resize(GeneNumID.size()); // remove later?

  for(size_t i = 0; i != geneModule.size(); ++i) {
    if( static_cast<int>(i) != GeneNumID[i])
      throw std::logic_error(" i != GeneNumID. Bug in R code.");
    // remove the next two these later?
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


fitnessLandscape_struct convertFitnessLandscape(Rcpp::List flg,
                Rcpp::List fl_df, Rcpp::List full_FDF_spec,
                bool fdb, bool fdd){


  // We assume this:
  // A symbol can only be a symbol if it is
  // in the fitness landscape table. A genotype not in the fitness landscape table
  // always has fitness zero, and noting can depend on it.
  // A genotype of 0 genotype will never appear in an equation.


  fitnessLandscape_struct flS;

  flS.names = Rcpp::as<std::vector<std::string> >(flg["Gene"]);
  flS.NumID = Rcpp::as<std::vector<int> >(flg["GeneNumID"]);

  if (fdb) {

    flS.flbmap.clear(); //Set to 0 flmap
    std::vector<std::string> genotNames =
      Rcpp::as<std::vector<std::string> >(full_FDF_spec["Genotype_as_numbers"]);

    std::vector<std::string> fvarsbvect =
    Rcpp::as<std::vector<std::string> > (full_FDF_spec["Genotype_as_fvarsb"]);

    std::vector<std::string> birth =
      Rcpp::as<std::vector<std::string> > (full_FDF_spec["Birth_as_fvars"]);

    if(fvarsbvect.size() != genotNames.size() )
      throw std::logic_error("fvarsbvect (fitnessLandscapeVariables) and "
          "genotNames (fitnessLandscape_df$Genotypes) "
          "are of different lenght. Should have been caught in R");
    
    // Fill up the map genotypes (as string of ints) to fitness (as string:
    // the expression that exprTk will evaluate).
    // Length given by number of genotypes in fitness landscape
    for(size_t i = 0; i != genotNames.size(); ++i) {
      flS.flFDBmap.insert({genotNames[i], birth[i]});
    }

    // Fill up the map genotypes (as string of ints) to fVars (as string)
    // If a mismatch in $fitnessLandscape_df and
    // $fitnessLandscapeVariables, this does silly things.
    for(size_t i = 0; i != genotNames.size(); ++i) {
      flS.flfVarsBmap.insert({genotNames[i], fvarsbvect[i]});
    }
  } else {

    flS.flFDBmap.clear();//Set to 0 flFDBmap
    flS.flfVarsBmap.clear();//Set to 0 flfVarsBmap

    std::vector<std::string> genotNames =
      Rcpp::as<std::vector<std::string> >(fl_df["Genotype"]);
    Rcpp::NumericVector birth = fl_df["Birth"];

    for(size_t i = 0; i != genotNames.size(); ++i) {
      flS.flbmap.insert({genotNames[i], birth[i]});
    }
  }

  // Check if death is specified
  if (fl_df.containsElementNamed("Death")) {

    if(fdd) {

      flS.fldmap.clear(); //Set to 0 fldmap

      std::vector<std::string> genotNames =
      Rcpp::as<std::vector<std::string> >(full_FDF_spec["Genotype_as_numbers"]);

      std::vector<std::string> fvarsdvect =
      Rcpp::as<std::vector<std::string> > (full_FDF_spec["Genotype_as_fvarsd"]);

      std::vector<std::string> death =
        Rcpp::as<std::vector<std::string> > (full_FDF_spec["Death_as_fvars"]);


      if(fvarsdvect.size() != genotNames.size() )
        throw std::logic_error("fvarsdvect (fitnessLandscapeVariables) and "
            "genotNames (fitnessLandscape_df$Genotypes) "
            "are of different lenght. Should have been caught in R");
      
      // Fill up the map genotypes (as string of ints) to fitness (as string:
      // the expression that exprTk will evaluate).
      // Length given by number of genotypes in fitness landscape
      for(size_t i = 0; i != genotNames.size(); ++i) {
        flS.flFDDmap.insert({genotNames[i], death[i]});
      }

      // Fill up the map genotypes (as string of ints) to fVars (as string)
      // If a mismatch in $fitnessLandscape_df and
      // $fitnessLandscapeVariables, this does silly things.
      for(size_t i = 0; i != genotNames.size(); ++i) {
        flS.flfVarsDmap.insert({genotNames[i], fvarsdvect[i]});
      }
    } else {

      flS.flFDDmap.clear();//Set to 0 flFDDmap
      flS.flfVarsDmap.clear();//Set to 0 flfVarsDmap

      std::vector<std::string> genotNames =
        Rcpp::as<std::vector<std::string> >(fl_df["Genotype"]);

      Rcpp::NumericVector death = fl_df["Death"];

      for(size_t i = 0; i != genotNames.size(); ++i) {
        flS.fldmap.insert({genotNames[i], death[i]});
      }
    }

  } else {
    flS.flFDDmap.clear();//Set to 0 flFDDmap
    flS.flfVarsDmap.clear();//Set to 0 flfVarsDmap
    flS.fldmap.clear(); //Set to 0 fldmap
  }

  return flS;
}

genesWithoutInt convertNoInts(Rcpp::List nI) {
  genesWithoutInt genesNoInt;
  genesNoInt.names = Rcpp::as<std::vector<std::string> >(nI["Gene"]);
  genesNoInt.NumID = Rcpp::as<std::vector<int> >(nI["GeneNumID"]);
  genesNoInt.s = Rcpp::as<std::vector<double> >(nI["s"]);
  genesNoInt.shift = genesNoInt.NumID[0]; // we assume mutations always
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

std::vector<int> sortedAllOrder(const std::vector<epistasis>& E) {

  std::vector<int> allG;
  for(auto const &ec : E) {
    for(auto const &g : ec.NumID) {
      allG.push_back(g);
    }
  }
  sort(allG.begin(), allG.end());
  allG.erase( unique( allG.begin(), allG.end()),
		      allG.end());
  return allG;
}

std::vector<int> sortedAllPoset(const std::vector<Poset_struct>& Poset) {
  // Yes, this could be done inside rTable_to_Poset but this is cleaner
  // and will only add very little time.
  std::vector<int> allG;
  for(auto const &p : Poset) {
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
  bool rone = as<bool>(rFE["gMOneToOne"]);
  Rcpp::IntegerVector drv = rFE["drv"];
  bool fdb = as<bool>(rFE["frequencyDependentBirth"]);
  bool fdd = as<bool>(rFE["frequencyDependentDeath"]);
  std::string fType = as<std::string>(rFE["frequencyType"]);
  Rcpp::List flg = rFE["fitnessLandscape_gene_id"];
  // clang does not like this: Rcpp::DataFrame fl_df = rFE["fitnessLandscape_df"];
  Rcpp::List fl_df = rFE["fitnessLandscape_df"];

  // In the future, if we want noInt and fitnessLandscape, all
  // we need is use the fitness landscape with an index smaller than those
  // of noInt. So we can use noInt with shift being those in fitnessLandscape.
  // BEWARE: will need to modify also createNewGenotype.
  //<std::vector<std::string> > fvariables = as<std::vector<std::string> > (fvars);
  // if(fl_df.nrows()) {
  if(fl_df.size()) {
    
    Rcpp::List full_FDF_spec;
    if (fdb || fdd) {
      // Should not be used in convertFitnessLandscape if not fdb and not fdd
      full_FDF_spec = rFE["full_FDF_spec"];
    }

    fe.fitnessLandscape = convertFitnessLandscape(flg, fl_df, full_FDF_spec, fdb, fdd); 
    if (fdb) {
      fe.fVarsb = as<std::vector<std::string> > (full_FDF_spec["Genotype_as_fvarsb"]);
      fe.frequencyType = fType;
    }

    if (fdd) {
      fe.fVarsd = as<std::vector<std::string> > (full_FDF_spec["Genotype_as_fvarsd"]);
      fe.frequencyType = fType;
    }
  }

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
  // If this is null, use the nullFitnessEffects function; never
  // end up here.
  if( (rrt.size() + re.size() + ro.size() + rgi.size() + fl_df.size()) == 0) {
      throw std::logic_error("\n Nothing inside this fitnessEffects; why are you here?"
			     "  Bug in R code.");
  }

  // At least for now, if we use fitness landscape nothing else allowed
  if(fl_df.size() && ((rrt.size() + re.size() + ro.size() + rgi.size()) > 0)) {
    throw std::logic_error("\n Fitness landscape specification."
			   " There should be no other terms. "
			   " Bug in R code");
  }

  fe.Gene_Module_tabl = R_GeneModuleToGeneModule(rgm);
  fe.allOrderG = sortedAllOrder(fe.orderE);
  fe.allPosetG = sortedAllPoset(fe.Poset);
  fe.gMOneToOne = rone;
  fe.allGenes = allGenesinFitness(fe);
  fe.genomeSize =  fe.Gene_Module_tabl.size() - 1 + fe.genesNoInt.s.size() +
    fe.fitnessLandscape.NumID.size();
  fe.drv = as<std::vector<int> > (drv);
  sort(fe.drv.begin(), fe.drv.end()); //should not be needed, but just in case
  // cannot trust R gives it sorted

  fe.frequencyDependentBirth = fdb; //new line to insert frequencyDependentBirth
  fe.frequencyDependentDeath = fdd; //new line to insert frequencyDependentDeath

  if(fe.genomeSize != static_cast<int>(fe.allGenes.size())) {
    throw std::logic_error("\n genomeSize != allGenes.size(). Bug in R code.");
  }
  // At least for now
  if(fe.fitnessLandscape.NumID.size() > 0) {
    if(fe.genomeSize != static_cast<int>(fe.fitnessLandscape.NumID.size())) {
      throw std::logic_error("\n genomeSize != genes in fitness landscape."
			     "Bug in R code.");
    }
  }
  return fe;
}

// Before making allGenesinGenotype return a unique vector: we do a
// set_difference below. If we look at the help
// (http://en.cppreference.com/w/cpp/algorithm/set_difference) if we had
// more repetitions of an element in allGenes than in sortedparent we
// could have a problem. But if you look at function "allgenesinFitness",
// which is the one used to give the allgenes vector, you will see that
// that one returns only one entry per gene, as it parses the geneModule
// structure. So even if allGenesinGenotype returns multiple entries
// (which they don't), there will be no bugs as the maximum number of
// entries in the output of setdiff will be 0 or 1 as m is 1. But
// allGenesinGenotype cannot return more than one element as can be seen
// in createNewGenotype: an element only ends in the order component if it
// is not in the epistasis component.  So there are no repeated elements
// in allGenes or in sortedparent below. Also, beware that does not break
// correct fitness evaluation of fitnessEffects where the same gene is in
// epistasis and order, as can be seen in the tests and because of how we
// evaluate fitness, where genes in a genotype in the orderEffects bucket
// are placed also in the epist for fitness eval. See evalGenotypeFitness
// and createNewGenotype.
void obtainMutations(const Genotype& parent,
		     const fitnessEffectsAll& fe,
		     int& numMutablePosParent,
		     std::vector<int>& newMutations,
		     //randutils::mt19937_rng& ran_gen
		     std::mt19937& ran_gen,
		     std::vector<double> mu) {
  //Ugly: we return the mutations AND the numMutablePosParent This is
  // almost ready to accept multiple mutations. And it returns a vector,
  // newMutations.
  std::vector<int> sortedparent = allGenesinGenotype(parent);
  std::vector<int> nonmutated;
  set_difference(fe.allGenes.begin(), fe.allGenes.end(),
		 sortedparent.begin(), sortedparent.end(),
		 back_inserter(nonmutated));
  // numMutablePos is used not only for mutation but also to decide about
  // the dummy or null mutation case.
  numMutablePosParent = nonmutated.size();
  if(nonmutated.size() < 1)
    throw std::out_of_range("Trying to obtain a mutation when nonmutated.size is 0."
			    " Bug in R code; let us know.");
  if(mu.size() == 1) { // common mutation rate
    //chromothripsis would not use this, or this is the limit case with a
    // single mutant
    std::uniform_int_distribution<int> rpos(0, nonmutated.size() - 1);
    newMutations.push_back(nonmutated[rpos(ran_gen)]);
  } else { // per-gene mutation rate.
    // Remember that mutations always indexed from 1, not from 0.
    // We take an element from a discrete distribution, with probabilities
    // proportional to the rates.
    // FIXMEmaybe:varmutrate give a warning if the only mu is for mu = 0?
    std::vector<double> mu_nm;
    for(auto const &nm : nonmutated) mu_nm.push_back(mu[nm - 1]);
    std::discrete_distribution<int> rpos(mu_nm.begin(), mu_nm.end());
    newMutations.push_back(nonmutated[rpos(ran_gen)]);
  }
  // randutils
  // // Yes, the next will work, but pick is simpler!
  // // size_t rpos = ran_gen.uniform(static_cast<size_t>(0), nonmutated.size() - 1);
  // //  newMutations.push_back(nonmutated[rpos]);
  // int posmutated = ran_gen.pick(nonmutated);
  // newMutations.push_back(posmutated);
}


// std::vector<int> genesInOrderModules(const fitnessEffectsAll& fe) {
//   vector<int> genes;
//   if(fe.gMOneToOne)
//     genes = fe.alOrderG;
//   else {
//     for(auto const &m : fe.allOrderG) {
//       genes.push_back()
//     }

//   }
//   return genes;
// }


fitness_as_genes fitnessAsGenes(const fitnessEffectsAll& fe) {
  // Give the fitnessEffects in terms of genes, not modules.

  // Extract the noInt. Then those in order effects by creating a multimap
  // to go from map to genes. Then all remaining genes are those only in
  // poset. By set_difference.
  fitness_as_genes fg = zero_fitness_as_genes();

  // fitness_as_genes fg;
  fg.flGenes = fe.fitnessLandscape.NumID;
  if(fg.flGenes.size()) {
    return fg;
  }

  fg.noInt = fe.genesNoInt.NumID;
  std::multimap<int, int> MG;
  for( auto const &mt : fe.Gene_Module_tabl) {
    MG.insert({mt.ModuleNumID, mt.GeneNumID});
  }
  for (auto const &o : fe.allOrderG) {
    for(auto pos = MG.lower_bound(o); pos != MG.upper_bound(o); ++pos)
      fg.orderG.push_back(pos->second);
  }
  sort(fg.orderG.begin(), fg.orderG.end());

  std::vector<int> tmpv = fg.orderG;
  tmpv.insert(tmpv.end(),fg.noInt.begin(), fg.noInt.end());
  sort(tmpv.begin(), tmpv.end()); // should not be needed

  set_difference(fe.allGenes.begin(), fe.allGenes.end(),
		 tmpv.begin(), tmpv.end(),
		 back_inserter(fg.posetEpistG));
  return fg;
}


std::map<int, std::string> mapGenesIntToNames(const fitnessEffectsAll& fe) {
  // This is a convenience, used in the creation of output.
  // Sure, we could do this when reading the data in.
  // The noInt in convertNoInts.
  std::map<int, std::string> gg;

  for(auto const &mt : fe.Gene_Module_tabl) {
    gg.insert({mt.GeneNumID, mt.GeneName});
  }
  // this is pedantic, as what is the size_type of NumID and of names?
  // for(decltype(fe.genesNoInt.s.size()) i = 0;
  //     i != fe.genesNoInt.s.size(); ++i)

  for(size_t i = 0;
      i != fe.genesNoInt.NumID.size(); ++i){
    gg.insert({fe.genesNoInt.NumID[i], fe.genesNoInt.names[i]});
  }

  for(size_t i = 0;
      i != fe.fitnessLandscape.NumID.size(); ++i){
    gg.insert({fe.fitnessLandscape.NumID[i], fe.fitnessLandscape.names[i]});
  }


  return gg;
}

// It is simple to write specialized functions for when
// there are no restrictions or no order effects , etc.
Genotype createNewGenotype(const Genotype& parent,
			   const std::vector<int>& mutations,
			   const fitnessEffectsAll& fe,
			   std::mt19937& ran_gen,
			   //randutils::mt19937_rng& ran_gen
			   bool random = true) {
  // random: if multiple mutations, randomly shuffle the ordered ones?
  // This is the way to go if chromothripsis, but not if we give an
  // initial mutant

  Genotype newGenot = parent;
  std::vector<int> tempOrder; // holder for multiple muts if order.
  bool sort_rest = false;
  bool sort_epist = false;
  bool sort_flgenes = false;

  // Order of ifs: I suspect order effects rare. No idea about
  // non-interaction genes, but if common the action is simple.

  // A gene that is involved both in order effects and epistasis, only
  // ends up in the orderEff container, for concision (even if a gene can
  // be involved in both orderEff and epistasis and rT). But in the
  // genotype evaluation, in evalGenotypeFitness, notice that we create
  // the vector of genes to be checked against epistais and order effects
  // using also those from orderEff:
  //   std::vector<int> mutG (ge.epistRtEff);
  //   mutG.insert( mutG.end(), ge.orderEff.begin(), ge.orderEff.end());
  for(auto const &g : mutations) {
    // If we are dealing with a fitness landscape, that is as far as we go here
    // at least for now. No other genes affect fitness.
    // But this can be easily fixed in the future; like this?
    // if(g <= (fe.fitnessLandscape.NumID.size() + 1)) {
    // and restructure the else logic for the noInt
    if(fe.fitnessLandscape.NumID.size()) {
      newGenot.flGenes.push_back(g);
      sort_flgenes = true;
    } else {
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
  }

  // If there is order but multiple simultaneous mutations
  // (chromothripsis), we randomly insert them

  // initMutant cannot use this: we give the order.
  // That is why the call from initMutant uses random = false
  if( (tempOrder.size() > 1) && random)
    shuffle(tempOrder.begin(), tempOrder.end(), ran_gen);
  // The new randutils engine:
  // if(tempOrder.size() > 1)
  //   ran_gen.shuffle(tempOrder.begin(), tempOrder.end());


  for(auto const &g : tempOrder)
    newGenot.orderEff.push_back(g);

  // Sorting done at end, in case multiple mutations
  if(sort_rest)
    sort(newGenot.rest.begin(), newGenot.rest.end());
  if(sort_epist)
    sort(newGenot.epistRtEff.begin(), newGenot.epistRtEff.end());
  if(sort_flgenes)
    sort(newGenot.flGenes.begin(), newGenot.flGenes.end());
  return newGenot;

  // Prepare specialized functions:??
  // Specialized functions:
  // Never interactions: push into rest and sort. Identify by shift == 1.
  // Never no interactions: remove the if. shift == -9.
}



// A paranoid check in case R code has bugs.
void breakingGeneDiff(const vector<int>& genotype,
		      const vector<int>& fitness) {
  std::vector<int> diffg;

  set_difference(genotype.begin(), genotype.end(),
		 fitness.begin(), fitness.end(),
		 back_inserter(diffg));
  if(diffg.size()) {
    Rcpp::Rcout << "Offending genes :";
    for(auto const &gx : diffg) {
      Rcpp::Rcout << " " << gx;
    }
    Rcpp::Rcout << "\t Genotype: ";
    for(auto const &g1 : genotype) Rcpp::Rcout << " " << g1;
    Rcpp::Rcout << "\t Fitness: ";
    for(auto const &g1 : fitness) Rcpp::Rcout << " " << g1;

    Rcpp::Rcout << "\n ";
    throw std::logic_error("\n At least one gene in the genotype not in fitness effects."
			   " Bug in R code.");

  }
}

void checkNoNegZeroGene(const vector<int>& ge) {
  if( ge[0] == 0 )
    throw std::logic_error("\n Genotype cannot contain 0. Bug in R code.");
  else if(ge[0] < 0)
    throw std::logic_error("\n Genotype cannot contain negative values. Bug in R code.");
}

void checkLegitGenotype(const Genotype& ge,
			const fitnessEffectsAll& F) {
  if((ge.orderEff.size() + ge.epistRtEff.size() + ge.rest.size()) == 0) {
    // An empty genotype is always legitimate, even if silly
    return;
  }
  // vector<int> g0 = allGenesinFitness(F);
  vector<int> g0 = F.allGenes;
  vector<int> allgG = allGenesinGenotype(ge);
  checkNoNegZeroGene(allgG);
  breakingGeneDiff(allgG, g0);
}

void checkLegitGenotype(const vector<int>& ge,
			const fitnessEffectsAll& F) {
  if (ge.size() == 0) {
    // An empty genotype is always legitimate, even if silly
    return;
  }
  // std::vector<int> g0 = allGenesinFitness(F);
  vector<int> g0 = F.allGenes;
  std::vector<int> allgG (ge);
  sort(allgG.begin(), allgG.end());
  checkNoNegZeroGene(allgG);
  breakingGeneDiff(allgG, g0);
}



Genotype convertGenotypeFromInts(const std::vector<int>& gg,
				 const fitnessEffectsAll& fe) {
  // A genotype is of one kind or another depending on what genes are of
  // what type.
  Genotype newGenot;

  if(gg.size() != 0) {
    // check_disable_later
    checkLegitGenotype(gg, fe);

    // Very similar to logic in createNewGenotype for placing each gene in
    // its correct place, which needs to look at module mapping.
    for(auto const &g : gg) {
      if(fe.fitnessLandscape.NumID.size()) {
	      newGenot.flGenes.push_back(g);
      } else {
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
    }

    sort(newGenot.flGenes.begin(), newGenot.flGenes.end());
    sort(newGenot.rest.begin(), newGenot.rest.end());
    sort(newGenot.epistRtEff.begin(), newGenot.epistRtEff.end());
  } else {
    newGenot = wtGenotype(); // be explicit!!
  }
  return newGenot;
}


Genotype convertGenotypeFromR(Rcpp::IntegerVector rG,
			      const fitnessEffectsAll& fe) {
  std::vector<int> gg = Rcpp::as<std::vector<int> > (rG);
  return convertGenotypeFromInts(gg, fe);
}


bool match_order_effects(const std::vector<int>& O,
			 const std::vector<int>& G) {
  //As the name says: we check if the order effect is matched
  if(G.size() < O.size()) return false;

  std::vector<int>::const_iterator p;
  std::vector<size_t> vdist;

  auto itb = G.begin();

  for(auto const &o : O) {
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
  for(auto const &o : OE) {
    if(match_order_effects(o.NumID, mutatedM))
      s.push_back(o.s);
  }
  return s;
}


bool match_negative_epist(const std::vector<int>& E,
			  const std::vector<int>& G) {
  // When we have things like -1, 2 in epistasis. We need to check 2 is
  // present and 1 is not present. E is the vector of epistatic coeffs,
  // and G the sorted genotype.

  if(G.size() < 1) return false;

  for(auto const &e : E) {
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
    throw std::logic_error("mutatedModules not sorted in evalEpistasis."
			   " Bug in R code.");

  for(auto const &p : Epistasis ) {
    if(p.NumID[0] > 0 ) {
      if(includes(mutatedModules.begin(), mutatedModules.end(),
		   p.NumID.begin(), p.NumID.end()))
	s.push_back(p.s);
    } else {
      if(match_negative_epist(p.NumID, mutatedModules))
	s.push_back(p.s);
    }
  }
  return s;
}

std::vector<double> evalPosetConstraints(const std::vector<int>& mutatedModules,
					 const std::vector<Poset_struct>& Poset,
					 const std::vector<int>& allPosetG) {

  // check_disable_later
  if(! is_sorted(mutatedModules.begin(), mutatedModules.end()))
    throw std::logic_error("mutatedModules not sorted in evalPosetConstraints."
			   " Bug in R code.");

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

  std::vector<int> MPintersect;


  std::set_intersection(allPosetG.begin(), allPosetG.end(),
			mutatedModules.begin(), mutatedModules.end(),
			std::back_inserter(MPintersect));

  // We know MPintersect is sorted, so we can avoid an O(n*n) loop
  size_t i = 0;
  for(auto const &m : MPintersect) {
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
  return s;
}

double frequency(const int& pos, const std::vector<spParamsP>& popParams){

  double fqc;
  double numerator = popParams[pos].popSize;
  double denominator = 0.0;

  for(size_t i = 0; i < popParams.size(); i++){
      denominator += popParams[i].popSize;
  }

  fqc = numerator / denominator;
  return fqc;
}

//This function returns 0 if genotype is not in Genotypes and position+1 otherwise
int findPositionInGenotypes(const std::vector<Genotype>& Genotypes,
	const std::vector<int> genotype){

		int index;

		std::vector<std::vector<int> > flGenesInGenotypes;
	for(size_t i = 0; i < Genotypes.size(); i++){
		flGenesInGenotypes.push_back(Genotypes[i].flGenes);
	}

	int pos = std::find(flGenesInGenotypes.begin(),
	  flGenesInGenotypes.end(), genotype) - flGenesInGenotypes.begin();

	int size = flGenesInGenotypes.size();

	if(pos < size)
		index = pos + 1;
	else
		index = 0;

  return index;
}

//This function works oposite to concatIntsString #include <sstream>
std::vector<int> stringVectorToIntVector(const std::string str){

  std::vector<int> vector;
  std::stringstream ss(str);

  int i;
  while (ss >> i) {
        vector.push_back(i);
        if(ss.peek() == ','|| ss.peek() == ' '){
            ss.ignore();
        }
  }
  return vector;
}


//This function produces the map that links fVars (keys) to its frequencies (values)
std::map<std::string, double> evalFVars(const fitnessEffectsAll& F,
					const std::vector<Genotype>& Genotypes,
					const std::vector<spParamsP>& popParams,
          bool birth){
  
  // Added birth, to indicate if you want the fvars of birth if true or the fvars of death
  // if false

  std::map<std::string, double> mapFvarsValues;
  std::map<std::string, std::string> fvarsmap;
  std::string freqType = F.frequencyType;

  if(birth) {
    fvarsmap = F.fitnessLandscape.flfVarsBmap;
    
  } else {
    fvarsmap = F.fitnessLandscape.flfVarsDmap;
  }
  
  
  for(const auto& iterator : fvarsmap) {
    std::vector<int> genotype = stringVectorToIntVector(iterator.first);//genotype (as int vector)
    std::string var = iterator.second;//variable associated to genotype
    int position = findPositionInGenotypes(Genotypes, genotype);
    
    if(position != 0){
      int realPos = position - 1;
      if(freqType == "abs"){
          double freqAbs = popParams[realPos].popSize;
          mapFvarsValues.insert({var, freqAbs});
      } else {
          double freqRel = frequency(realPos, popParams);
          mapFvarsValues.insert({var, freqRel});
      }
    } else {
      double freq = 0.0;
      mapFvarsValues.insert({var, freq});
    }
  }
  return mapFvarsValues;
}

//This function produces the map that links fVars (keys) to its fitness/birth/death/mutation (values)
std::map<std::string, double> evalFVarsFitness(const fitnessEffectsAll& F,
					const std::vector<Genotype>& Genotypes,
					const std::vector<spParamsP>& popParams,
          std::string value){
  
  // Added birth, to indicate if you want the fvars of birth if true or the fvars of death
  // if false

  std::map<std::string, double> mapFvarsValues;
  std::map<std::string, std::string> fvarsmap;

  fvarsmap = F.fitnessLandscape.flfVarsBmap;  
  
  for(const auto& iterator : fvarsmap) {
    std::vector<int> genotype = stringVectorToIntVector(iterator.first);//genotype (as int vector)
    std::string var = iterator.second;//variable associated to genotype
    double val = 0.0;
    int position = findPositionInGenotypes(Genotypes, genotype);
    
    if(position != 0){
      int realPos = position - 1;
      if(value == "birth"){
        val = popParams[realPos].birth;
      }else if(value == "death"){
        val = popParams[realPos].death;
      }else if(value == "mutation"){
        val = popParams[realPos].mutation;
      }else{
        std::string errorMessage = "Invalid value to evalFVarsFitness";
        throw std::invalid_argument(errorMessage);
      }  
    } 

    var = var.substr(2);

    mapFvarsValues.insert({var, val});
  }

  return mapFvarsValues;
}


double totalPop(const std::vector<spParamsP>& popParams){
	double sum = 0.0;
	for(size_t i = 0; i < popParams.size(); i++){
      sum += popParams[i].popSize;
  }
	return sum;
}

double evalGenotypeFDFitnessEcuation(const Genotype& ge,
				     const fitnessEffectsAll& F,
				     const std::vector<Genotype>& Genotypes,
                                     const std::vector<spParamsP>& popParams,
                                     const double& currentTime,
                                     bool birth){

  double f;

  std::map<std::string, double> mapFvarsValues = evalFVars(F, Genotypes, popParams, birth);

  std::string gs = concatIntsString(ge.flGenes);
  std::string expr_string;

  if(birth) {
    expr_string = F.fitnessLandscape.flFDBmap.at(gs);
  }
  else {
    expr_string = F.fitnessLandscape.flFDDmap.at(gs);
  }
  

  double N = totalPop(popParams);
  double T = currentTime;
 
  typedef exprtk::symbol_table<double> symbol_table_t;
  typedef exprtk::expression<double> expression_t;
  typedef exprtk::parser<double> parser_t;

  symbol_table_t symbol_table;
  for(auto& iterator : mapFvarsValues) {
    symbol_table.add_variable(iterator.first, iterator.second);
  }
  symbol_table.add_constant("N", N);//We reserve N to total population size
  symbol_table.add_constant("T", T); //Pass current time to exprtk
  symbol_table.add_constants();

  expression_t expression;
  expression.register_symbol_table(symbol_table);
  
  parser_t parser;
  
  if (!parser.compile(expr_string, expression)){
    Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;
    for (std::size_t i = 0; i < parser.error_count(); ++i){
      typedef exprtk::parser_error::type error_t;
      error_t error = parser.get_error(i);
      // FIXMEmaybe: Use warning or error to capture it easily in tests?
      REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
	       i,
	       error.token.position,
	       exprtk::parser_error::to_str(error.mode).c_str(),
	       error.diagnostic.c_str(),
	       expr_string.c_str());
    }
    std::string errorMessage1 = "Wrong evalGenotypeFDFitnessEcuation evaluation, ";
    std::string errorMessage2 = "probably bad fitness columm especification.";
    std::string errorMessage = errorMessage1 + errorMessage2;
    throw std::invalid_argument(errorMessage);
  }
  f = expression.value();
  return f;
}

// Now that death can be specified, if death is specified it
// will return a tuple of numbers with the birth evaluation and
// the death evaluation. Otherwise, it will return a list with
// just the birth evaluation
std::vector<double> evalGenotypeFitness(const Genotype& ge,
	const fitnessEffectsAll& F,
	const std::vector<Genotype>& Genotypes,
	const std::vector<spParamsP>& popParams,
	const double& currentTime){
  
  // should we disabe this check? check_disable_later
  checkLegitGenotype(ge, F);

  std::vector<double> s;

  if( ((ge.orderEff.size() + ge.epistRtEff.size() +
	ge.rest.size() + ge.flGenes.size() ) == 0) && !F.frequencyDependentBirth &&
  !F.frequencyDependentDeath ) {
    s.push_back(0.0);
    return s; //return 0.0 explicitly?? FIXME
  }

  // If we are dealing with a fitness landscape, that is as far as we go here
  // at least for now. No other genes affect fitness.
  // But this can be easily fixed in the future; do not return
  // s below, but keep adding, maybe the noIntGenes.
  // Recall also  prodFitness uses, well, the prod of 1 + s
  // so we want an s s.t. 1 + s = birth rate passed,
  // which is the value in the fitness landscape as interpreted now.
  // i.e., s = birth rate - 1;
  if(F.fitnessLandscape.NumID.size()) {
    std::string gs = concatIntsString(ge.flGenes);
    if(F.frequencyDependentBirth){ //possible also with Genotype.size()==0 and popParams.size==0 ?
      if(F.fitnessLandscape.flFDBmap.find(gs) == F.fitnessLandscape.flFDBmap.end()) {
	      s.push_back(-1.0);
      } else {
	      s.push_back(evalGenotypeFDFitnessEcuation(ge, F, Genotypes, popParams, currentTime, true) - 1);
      }
    } else {
      if(F.fitnessLandscape.flbmap.find(gs) == F.fitnessLandscape.flbmap.end()) {
	      s.push_back(-1.0);
      } else {
	      s.push_back(F.fitnessLandscape.flbmap.at(gs) - 1);
      }
    }

    if(F.frequencyDependentDeath){ 
      if(F.fitnessLandscape.flFDDmap.find(gs) == F.fitnessLandscape.flFDDmap.end()) {
	      s.push_back(-1.0);
      } else {
	      s.push_back(evalGenotypeFDFitnessEcuation(ge, F, Genotypes, popParams, currentTime, false) - 1);
      }
    } else {
      if (F.fitnessLandscape.fldmap.size()) { // Check if death was specified. If not, dont do anything
        if(F.fitnessLandscape.fldmap.find(gs) == F.fitnessLandscape.fldmap.end()) {
                s.push_back(-1.0);
        } else {
          s.push_back(F.fitnessLandscape.fldmap.at(gs) - 1);
        }
      }
    }
  }

  // Genes without any restriction or epistasis are just genes. No modules.
  // So simple we do it here.
  if(F.genesNoInt.shift > 0) {
    int shift = F.genesNoInt.shift;
    for(auto const & r : ge.rest ) {
      s.push_back(F.genesNoInt.s[r - shift]);
    }
  }

  // For the rest, there might be modules. Three different effects on
  // fitness possible: as encoded in Poset, general epistasis, order effects.

  // Epistatis and poset are checked against all mutations. Create single
  // sorted vector with all mutations and map to modules, if needed. Then
  // eval.

  // Why not use a modified genotypeSingleVector without the no ints? We
  // could, but not necessary. And you can place genes in any order you
  // want, since this is not for order restrictions. That goes below.
  // Why do I put the epist first? See previous answer.
  // Why do I sort if one to one? binary searches. Not done below for order.
  std::vector<int> mutG (ge.epistRtEff);
  // A gene can be involved in epistasis and order. This gene would only
  // be in the orderEff vector, as seen in "createNewGenotype" or
  // "convertGenotypeFromInts"
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

vector<int> getGenotypeDrivers(const Genotype& ge, const vector<int>& drv) {
  // Returns the actual mutated drivers in a genotype.
  // drv comes from R, and it is the vector with the
  // numbers of the genes, not modules.
  vector<int> presentDrv;
  vector<int> og = allGenesinGenotype(ge);
  set_intersection(og.begin(), og.end(),
		   drv.begin(), drv.end(),
		   back_inserter(presentDrv));
  return presentDrv;
}


// FIXME: refactor?
// Why does this function take as arguments popParams and
// Genotypes?  It shouldn't, since its behavior does not depend on those.
// It does because we call evalGenotypeFitness which calls in turn
// evalGenotypeFDFitnessEcuation, etc. But we probably want to do this
// differently. From creating a new function, using overloading, to do s =
// evalGenotypeFitness, that would be double newfun(A, B) =
// evalGenotypeFitness(A, B); and maybe having two different
// evalGenotypeFitness.

// This same problem happens with currentTime; it is used by evalGenotypeFDFitnessEcuation,
// etc

double evalMutator(const Genotype& fullge,
		   const std::vector<int>& full2mutator,
		   const fitnessEffectsAll& muEF,
		   const std::vector<Genotype>& Genotypes,
		   const std::vector<spParamsP>& popParams,
		   const double& currentTime, // because we use it in evalGenotypeFitness
		   bool verbose = false) {
  // In contrast to nr_fitness, that sets birth and death, this simply
  // returns the multiplication factor for the mutation rate. This is used
  // by mutationFromParent and mutationFromScratch

  // Remember that the fitnessEffectsAll struct for mutator does not use
  // the same mapping from gene names to gene numerical IDs as for
  // fitness. fitnessEffectsAll and its associated algorithms expects the
  // present genes to be indexed as successive integers. That is not
  // necessarily the case if only some of the genes are in the mutator
  // fitnessEffectsAll. So we need to remap the gene numerical IDs.  We
  // could try remapping by gene inside the struct, but painful and
  // error-prone. Much simpler at least for now to do:

  // full genotype -> vector of ints (preserving order)
  //     -> convert the ints to the ints for the genotype of mutator
  //     -> genotype in terms of mutator

  // the "genotype in terms of mutator" is never preserved. It is just a
  // transient mapping.

  // This will NOT work if we ever have order effects for mutator as we do
  // not record order for those that matter for mutator if they do not matter for
  // fitness.

  vector<int> g1 = genotypeSingleVector(fullge);
  vector<int> g2;
  int tmp;
  for (auto const & i : g1) {
    tmp = full2mutator[i - 1]; //gives a -9 if no correspondence
    if( tmp > 0 ) g2.push_back(tmp);
  }
  if(g2.size() == 0) {
    return 1.0;
  } else {
    Genotype newg = convertGenotypeFromInts(g2, muEF);
    vector<double> s = evalGenotypeFitness(newg, muEF, Genotypes, popParams, currentTime);

    // just for checking
    if(verbose) {
      std::string sprod = "mutator product";
      Rcpp::Rcout << "\n Individual " << sprod << " terms are :";
      for(auto const &i : s) Rcpp::Rcout << " " << i;
      Rcpp::Rcout << std::endl;
    }
    return prodMuts(s);
  }
}

std::vector<Genotype> genotypesFromScratch(const std::vector<std::string>& genotNames){

  std::vector<Genotype> genotypesVector;

	for(size_t i = 0; i != genotNames.size(); ++i){
		//Genotype genotype = wtGenotype();
		Genotype genotype;
		genotype.flGenes = stringVectorToIntVector(genotNames[i]);
		genotypesVector.push_back(genotype);
	}

	return genotypesVector;
}

std::vector<spParamsP> popParamsFromScratch(const std::vector<int>& spPopSizes){

	std::vector<spParamsP> popParamsVector;

	for(size_t i = 0; i != spPopSizes.size(); ++i){
		spParamsP spparams;
		spparams.popSize = spPopSizes[i];
		popParamsVector.push_back(spparams);
	}

  return popParamsVector;
}

// [[Rcpp::export]]
Rcpp::NumericVector evalRGenotype(Rcpp::IntegerVector rG,
	Rcpp::List rFE,
	Rcpp::IntegerVector spPop,
	bool verbose,
	bool prodNeg,
	Rcpp::CharacterVector calledBy_,
	double currentTime) {
  // Can evaluate both ONLY fitness or ONLY mutator. Not both at the same
  // time. Use evalRGenotypeAndMut for that.

  const std::string calledBy = Rcpp::as<std::string>(calledBy_);
  const bool fdb = as<bool>(rFE["frequencyDependentBirth"]);
  const bool fdd = as<bool>(rFE["frequencyDependentDeath"]);
  const bool deathSpec = as<bool>(rFE["deathSpec"]);

  if(rG.size() == 0 && fdb == false && fdd == false) {
    Rcpp::Rcout << "NOTE: you have evaluated fitness/mutator status of a genotype of length zero  (WT?) in non fdf fitness. It is 1 by decree. \n";
    if(deathSpec) {
      return Rcpp::NumericVector::create(1.0, 1.0);
    }
    else {
      return Rcpp::NumericVector::create(1.0);
    }
  }

  std::vector<Genotype> Genotypes;
  std::vector<spParamsP> popParams;
  if(fdb || fdd){
    std::vector<int> spPopSizes = as<std::vector<int> > (spPop);
    // fitnessLandscape_df is the same as columns "Genotype_as_numbers" and "Fitness_as_fvars"
    // in full_FDF_spec
    Rcpp::List fl_df = rFE["fitnessLandscape_df"];
    std::vector<std::string> genotNames = Rcpp::as<std::vector<std::string> >(fl_df["Genotype"]);
    Genotypes = genotypesFromScratch(genotNames);
    popParams = popParamsFromScratch(spPopSizes);
  }
  
  fitnessEffectsAll F = convertFitnessEffects(rFE);
  Genotype g = convertGenotypeFromR(rG, F);
  vector<double> s = evalGenotypeFitness(g, F, Genotypes, popParams, currentTime);
  
  if(verbose) {
    std::string sprod;
    if(calledBy == "evalGenotype") {
      sprod = "s";
    } else { // if (calledBy == "evalGenotypeMut") 
      sprod = "mutator product";
    }
    Rcpp::Rcout << "\n Individual " << sprod << " terms are :";
    for(auto const &i : s) Rcpp::Rcout << " " << i;
    Rcpp::Rcout << std::endl;
  }
  if(calledBy == "evalGenotype") {
    if(!prodNeg) //for Bozic models. Comes from R.
        if (deathSpec) {
          // When death is specified, s is a vector with two elements:
          // the first element is the evaluation of the birth and the
          // second one is the evaluation of the death
          vector<double> birth(s.begin(), s.begin()+1);
          vector<double> death(s.begin()+1, s.end());

          return Rcpp::NumericVector::create(prodFitness(birth), prodFitness(death));

        } else {
          return Rcpp::NumericVector::create(prodFitness(s));
        }
         
    else {
      return Rcpp::NumericVector::create(prodDeathFitness(s));
    }
      
  } else { //if (calledBy == "evalGenotypeMut") 
    return Rcpp::NumericVector::create(prodMuts(s));
  }

}


// [[Rcpp::export]]
Rcpp::NumericVector evalRGenotypeAndMut(Rcpp::IntegerVector rG,
					Rcpp::List rFE,
					Rcpp::List muEF,
					Rcpp::IntegerVector spPop,
					Rcpp::IntegerVector full2mutator_,
					bool verbose,
					bool prodNeg,
					double currentTime) {
  // Basically to test evalMutator. We repeat the conversion to genotype,
  // but that is unavoidable here.

  const bool fdb = as<bool>(rFE["frequencyDependentBirth"]);
  const bool fdd = as<bool>(rFE["frequencyDependentDeath"]);
  
  NumericVector out(2);
  std::vector<Genotype> Genotypes;
  std::vector<spParamsP> popParams;
  
  if(fdb || fdd){
    std::vector<int> spPopSizes = as<std::vector<int> > (spPop);
    // fitnessLandscape_df is the same as columns "Genotype_as_numbers" and "Fitness_as_fvars"
    // in full_FDF_spec
    Rcpp::List fl_df = rFE["fitnessLandscape_df"];
    std::vector<std::string> genotNames = Rcpp::as<std::vector<std::string> >(fl_df["Genotype"]);
    Genotypes = genotypesFromScratch(genotNames);
    popParams = popParamsFromScratch(spPopSizes);
  }

  // For fitness. Except for "evalGenotypeFromR", all is done as in the
  // rest of the internal code for evaluating a genotype.
  fitnessEffectsAll F = convertFitnessEffects(rFE);
  fitnessEffectsAll muef = convertFitnessEffects(muEF);
  Genotype g = convertGenotypeFromR(rG, F);
  vector<double> s = evalGenotypeFitness(g, F, Genotypes, popParams, currentTime);
  
  if(!prodNeg)
    out[0] = prodFitness(s);
  else
    out[0] = prodDeathFitness(s);
  if(verbose) {
    std::string sprod = "s";
    Rcpp::Rcout << "\n Individual " << sprod << " terms are :";
    for(auto const &i : s) Rcpp::Rcout << " " << i;
    Rcpp::Rcout << std::endl;
  }

  const std::vector<int> full2mutator = Rcpp::as<std::vector<int> >(full2mutator_);
  out[1] = evalMutator(g, full2mutator, muef, Genotypes, popParams, currentTime, verbose);

  return out;
}

// FIXME refactor
// as in evalMutator, that this takes Genotypes and popParams is arguably
// bad design and very confusing (those arguments have nothing to do with mutations)
// Same problem happens with currentTime. Yes, this sucks.

// return mutation rate for a new genotype; set to dummyMutationRate in specific cases
//                                          do not use max w.r.t dummyMutationRate
//                                          indiscriminately: we want to see
//                                          dummyMutationRate being set
double mutationFromScratch(const std::vector<double>& mu,
			   const spParamsP& spP,
			   const Genotype& g,
			   const fitnessEffectsAll& fe,
			   const int mutationPropGrowth,
			   const std::vector<int> full2mutator,
			   const fitnessEffectsAll& muEF,
			   const std::vector<Genotype>& Genotypes,
			   const std::vector<spParamsP>& popParams,
			   const double& currentTime,
			   const double& dummyMutationRate) {

  double tmp;
  double mumult;
  if(full2mutator.size() > 0) { // so there are mutator effects
    mumult = evalMutator(g, full2mutator, muEF, Genotypes, popParams, currentTime);
  } else mumult = 1.0;
  
  
  // Here the code for altering mutation rate
  // with a procedure like ExprTk for fitness?
  // Nope: alter directly spParams.mutation.
  // where the updateRatesFDF... are called.
  // In BNB_nr.cpp, in nr_innerBNB function
  // when we are sampling.
  if(spP.numMutablePos == 0) {
    return(dummyMutationRate);
    }
  if(mu.size() == 1) {
    if(mutationPropGrowth) {
      tmp = mumult * mu[0] * spP.numMutablePos * spP.birth;
      if(tmp <= dummyMutationRate) {
        Rcpp::Rcout << "mFS-messagesMPL: Mutable positions left: ";
        if(mumult == 1.0) {
          // letters match codes for varmutrate
          Rcpp::Rcout << "mFS-message-2-B:  constant mut rate "
                << "no mutator and mutationPropGrowth "
                << ". birth rate = " << spP.birth << "\n";
        } else {
          Rcpp::Rcout << "mFS-message-2-C:  constant mut rate "
                << " mutator and mutationPropGrowth "
                << ". birth rate = " << spP.birth
                << ". mumult = " << mumult << "\n";
        }
        Rcpp::Rcout << "\n mutation rate = " << tmp << " < dummyMutationRate "
              << dummyMutationRate
              << ". Expect numerical problems.\n";
        // Remember genotypes of birth rate 0 never created
      }
      return(tmp);
    } else {
      tmp = mumult * mu[0] * spP.numMutablePos;
      if(tmp <= dummyMutationRate) {
	Rcpp::Rcout << "mFS-messagesMPL: Mutable positions left: ";
	if(mumult == 1.0) {
	  // letters match codes for varmutrate
	  Rcpp::Rcout << "mFS-message-2-A:  constant mut rate "
		      << "no mutator and no mutationPropGrowth ";
	} else {
	  Rcpp::Rcout << "mFS-message-2-D:  constant mut rate"
		      << " mutator and no mutationPropGrowth "
		      << ". mumult = " << mumult << "\n";
	}
	Rcpp::Rcout << "\n mutation rate = " << tmp << " < dummyMutationRate "
		    << dummyMutationRate
		    << ". Expect numerical problems.\n";
	  }
      return(tmp);
    }
  } else { // variable mutation rate
    std::vector<int> sortedG = allGenesinGenotype(g);
    std::vector<int> nonmutated;
    set_difference(fe.allGenes.begin(), fe.allGenes.end(),
		   sortedG.begin(), sortedG.end(),
		   back_inserter(nonmutated));
    double mutrate = 0.0;
    for(auto const &nm : nonmutated) {
      mutrate += mu[nm - 1];
    }
    if(mutrate == 0) {
      if(nonmutated.size() > 0)
	throw std::logic_error("\n This case should not exist: mFS-except-4-pre\n");
      if(nonmutated.size() == 0)
	throw std::logic_error("\n This case should not exist: mFS-except-4-post\n");
      // if(verbosity >= 1)
      // 	Rcpp::Rcout << "mFS-message-4 . No mutable positions. Mutation set to dummyMutationRate "
      // 		    << dummyMutationRate << "\n";
      // return(dummyMutationRate);
    }
   
    if(mutationPropGrowth)
      mutrate *= spP.birth;
    tmp = mumult * mutrate;

    if(tmp <= dummyMutationRate) {
      Rcpp::Rcout << "mFS-messagesMPL: Mutable positions left: ";
      if( (mumult == 1.0) && (!mutationPropGrowth) ) {
	Rcpp::Rcout << "mFS-message-5-A: variable mut rate "
		    << "no mutator and no mutationPropGrowth\n ";
      } else if ((mumult == 1.0) && (mutationPropGrowth) ) {
	Rcpp::Rcout << "mFS-message-5-B:  variable mut rate "
		    << "no mutator and mutationPropGrowth "
		    << ". birth rate = " << spP.birth << "\n";
      } else if ( (mumult != 1.0) && (mutationPropGrowth) ) {
	Rcpp::Rcout << "mFS-message-5-C:  variable mut rate "
		    << "mutator and mutationPropGrowth "
		    << ". birth rate = " << spP.birth
		    << ". mumult = " << mumult << "\n";
      } else if ( (mumult != 1.0) && (!mutationPropGrowth) ) {
	Rcpp::Rcout << "mFS-message-5-D:  variable mut rate "
		    << "mutator and no mutationPropGrowth "
		    << ". mumult = " << mumult << "\n";
      } else {
	throw std::logic_error("\n This case should not exist: mFS-except\n");
	  }
      
      Rcpp::Rcout << "\n mutation rate = " << tmp << " < dummyMutationRate "
		  << dummyMutationRate << ". Expect numerical problems.\n";
	}
    return(tmp);
  }
}


// list -> vector of int vectors.
//    If input is of length zero, return object of length 0
//    No checks are made of whether this makes sense (undefined behav
//    if this is a list of floats, strings, etc)
//    Check for order is necessary for fixation genotypes but not for init mutant.
std::vector < std::vector<int> > list_to_vector_of_int_vectors(Rcpp::List vlist,
							       bool check_ordered = true) {
  // As it says. We check each vector is sorted!
  std::vector < std::vector<int> > vv(vlist.size());
  if (vlist.size() != 0) {
    for(int i = 0; i != vlist.size(); ++i) {
      vv[i] = Rcpp::as<std::vector<int> >(vlist[i]);
      if( check_ordered && (! is_sorted(vv[i].begin(), vv[i].end())) )
	throw std::logic_error("Fixation genotypes not sorted. Bug in R code.");
    }
  } else {
    vv.resize(0);
  }
  return vv;
}














// // [[Rcpp::export]]
// void wrap_list_to_vector_of_int_vectors(Rcpp::List vlist) {
//   std::vector < std::vector<int> > vo(vlist.size());
//   vo = list_to_vector_of_int_vectors(vlist);
//   for(int ii = 0; ii != vo.size(); ++ii) {
//     Rcpp::Rcout << "\n";
//     Rcpp::Rcout << " list position " << ii + 1 << ": ";
//     for(int jj = 0; jj != vo[ii].size(); ++jj ) {
//       Rcpp::Rcout << vo[ii][jj] << " ";
//     }
//   }
//   Rcpp::Rcout << "\n";
// }



// Wrong when/if there are mutator effects. For suppose there wre, and
// they affected the parent, but no new mutator gene affects the child.
// We will, however, multiply twice by the mutator effect.  Therefore, we
// disable this for now.  We could fix this, checking if there are new
// mutation effects, or mutliplying/subtracting only new, etc. But too
// much of a mess.
// double mutationFromParent(const std::vector<double>& mu,
// 			  const spParamsP& newP,
// 			  const spParamsP& parentP,
// 			  const std::vector<int>& newMutations,
// 			  // const std::vector<int>& nonmutated,
// 			  const int mutationPropGrowth,
// 			  const Genotype& fullge,
// 			  const std::vector<int> full2mutator,
// 			  const fitnessEffectsAll& muEF) {
//   double mumult;
//   if(full2mutator.size() > 0) { // so there are mutator effects
//     mumult = evalMutator(fullge, full2mutator, muEF);
//   } else mumult = 1.0;

//   if(mu.size() == 1) {
//     if(mutationPropGrowth)
//       return(mumult * mu[0] * newP.numMutablePos * newP.birth);
//     else
//       return(mumult * mu[0] * newP.numMutablePos);
//   } else {
//     double mutrate = parentP.mutation;
//     for(auto const mutated : newMutations) {
//       mutrate -= mu[mutated - 1];
//     }
//     if(mutationPropGrowth)
//       mutrate *= newP.birth;
//     return(mumult * mutrate);
//   }
// }



// About order of genes and their names, etc

// We first read the R gene module. The $geneModule. The function is
// R_GeneModuleToGeneModule

// We also read the no interaction. They have their own number-name
// correspondence, within the noInt genes part. See the struct
// genesWithoutInt.  But that already comes order from R with numbers
// starting after the last gene with interaction. See the R function
// allFitnessEffects.


