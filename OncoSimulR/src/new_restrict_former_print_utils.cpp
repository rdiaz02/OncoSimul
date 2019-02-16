//     Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

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


// This code is not being used in production. It was used during
// development, to see what C++ was thinking about what R was passing to
// it. It is left here, but is now all commented out. 



// // [[Rcpp::export]]
// void readFitnessEffects(Rcpp::List rFE,
// 			bool echo) {
//   // fitnessEffectsAll fitnessEffects;
//   // convertFitnessEffects(rFE, fitnessEffects);
//   fitnessEffectsAll fitnessEffects = convertFitnessEffects(rFE);
//   if(echo) {
//      printFitnessEffects(fitnessEffects);
//   }
// }





// void printPoset(const std::vector<Poset_struct>& Poset) {

//   int counterInfs = 0;
//   int counterNegInfs = 0;
//   Rcpp::Rcout << "\n **********  Poset or Restriction table (internal) *******" 
// 	      << std::endl;
//   if(!Poset.size()) {
//     Rcpp::Rcout << "No posets: restriction table of size 0"<< std::endl;
//   } else {
//     Rcpp::Rcout << "Size = " << (Poset.size() - 1) << std::endl;
//     for(size_t i = 1; i != Poset.size(); ++i) {
//       // We do not show the Poset[0]
//       Rcpp::Rcout <<"\t Dependent Module or gene (child) " << i 
// 		  << ". childNumID: " << Poset[i].childNumID 
// 		  << ". child full name: " << Poset[i].child
// 		  << std::endl;
//       Rcpp::Rcout <<"\t\t typeDep = " << depToString(Poset[i].typeDep) << ' ' ;
//       Rcpp::Rcout <<"\t s = " << Poset[i].s << " ";
//       Rcpp::Rcout <<"\t sh = " << Poset[i].sh << std::endl;
//       if(std::isinf(Poset[i].sh))
// 	++counterInfs;
//       if(std::isinf(Poset[i].sh) && (Poset[i].sh < 0))
// 	++counterNegInfs;
//       Rcpp::Rcout << "\t\t Number of parent modules or genes = " << 
// 	Poset[i].parents.size() << std::endl;
//       Rcpp::Rcout << "\t\t\t Parents IDs: ";
//       for(auto const &c : Poset[i].parentsNumID)
// 	Rcpp::Rcout << c << "; ";
//       Rcpp::Rcout << std::endl;
//       Rcpp::Rcout << "\t\t\t Parents names: ";
//       for(auto const &c : Poset[i].parents)
// 	Rcpp::Rcout << c << "; ";
//       Rcpp::Rcout << std::endl;
    
//       // for(size_t j = 0; j != Poset[i].deps.size(); ++j) {
//       //   Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": " 
//       // 		  << Poset[i].deps[j] << std::endl;
//     }
//     Rcpp::Rcout << std::endl;

//     if(counterInfs) {
//       Rcpp::Rcout << "In sh there were " << counterNegInfs 
// 		  << " negative infinites and "
// 		  << (counterInfs - counterNegInfs) 
// 		  << " positive infinites" << std::endl;
//     }
//   }
// }


// void printGene_Module_table(const 
// 		       std::vector<Gene_Module_struct>& Gene_Module_tabl,
// 		       const bool gMOneToOne) {
//   // Rcpp::Rcout << 
//   //   "\n\n******** geneModule table (internal) *******:\nGene name\t Gene NumID\t Module name\t Module NumID\n";
//   // for(auto it = Gene_Module_tabl.begin(); it != Gene_Module_tabl.end(); ++it) {
//   //   Rcpp::Rcout << '\t' << it->GeneName << '\t' << it->GeneNumID << '\t' 
//   // 		<< it->ModuleName << '\t' << it->ModuleNumID << std::endl;
//   // }

//   Rcpp::Rcout << 
//     "\n\n******** geneModule table (internal) *******:\n" <<
//     std::setw(14) << std::left << "Gene name" << std::setw(14) << "Gene NumID" << std::setw(14)
// 	      << "Module name" << std::setw(14) << "Module NumID" << "\n";
//   for(auto it = Gene_Module_tabl.begin(); it != Gene_Module_tabl.end(); ++it) {
//     Rcpp::Rcout << std::setw(14) << std::left << it->GeneName << std::setw(14)
// 		<< it->GeneNumID << std::setw(14) << it->ModuleName
// 		<< std::setw(14) << it->ModuleNumID << std::endl;
//   }


//   if(gMOneToOne)
//     Rcpp::Rcout << "This is a dummy module table: each module is one gene."
// 		<< std::endl;
// }



// void printOtherEpistasis(const std::vector<epistasis>& Epistasis,
// 			 const std::string effectName,
// 			 const std::string sepstr) {
//   Rcpp::Rcout << "\n **********  General " << effectName << "s (internal) *******"
// 	      << std::endl;
//   if(!Epistasis.size()) {
//     Rcpp::Rcout << "No general " << effectName << std::endl;
//   } else {
//     Rcpp::Rcout << " Number of " << effectName <<"s = " << Epistasis.size();
//     for(size_t i = 0; i != Epistasis.size(); ++i) {
//       Rcpp::Rcout << "\n\t " << effectName << " " << i + 1 << ": " <<
// 	". Modules or Genes (names) = " << Epistasis[i].names[0];
//       for(size_t j = 1; j != Epistasis[i].NumID.size(); ++j) {
// 	Rcpp::Rcout << sepstr << Epistasis[i].names[j] ;
//       }
//       Rcpp::Rcout << ".\t Modules or Genes (NumID) = " << Epistasis[i].NumID[0];
//       for(size_t j = 1; j != Epistasis[i].NumID.size(); ++j) {
// 	Rcpp::Rcout << sepstr << Epistasis[i].NumID[j] ;
//       }
//       Rcpp::Rcout << ".\t s = " << Epistasis[i].s;
//     }
//   }
//   Rcpp::Rcout << std::endl;
// }

// void printNoInteractionGenes(const genesWithoutInt& genesNoInt) {
//   Rcpp::Rcout << "\n **********  All remaining genes without interactions (internal) *******"
// 	      << std::endl;
  
//   if(genesNoInt.shift <= 0) {
//     Rcpp::Rcout << "No other genes without interactions" << std::endl;
//   } else {
//     Rcpp::Rcout << std::setw(14) << std::left << "Gene name" << std::setw(14)
// 		<< "Gene NumID" << std::setw(14) << "s" << std::endl;
//     for(size_t i = 0; i != genesNoInt.NumID.size(); ++i) {
//       Rcpp::Rcout << std::setw(14) << std::left << genesNoInt.names[i]
// 		  << std::setw(14) << genesNoInt.NumID[i]	
// 		  << std::setw(14) << genesNoInt.s[i] << '\n';
//     }
//   }
// }

// void printAllOrderG(const std::vector<int> ge) {
//   Rcpp::Rcout << "\n **********  NumID of genes/modules in the order restrict. (internal) *******"
// 	      << std::endl;
//   for(auto const &g : ge)
//     Rcpp::Rcout << g << " ";
//   Rcpp::Rcout << std::endl;
// }


// void printFitnessEffects(const fitnessEffectsAll& fe) {
//   printGene_Module_table(fe.Gene_Module_tabl, fe.gMOneToOne);
//   printPoset(fe.Poset);
//   printOtherEpistasis(fe.orderE, "order effect", " > ");
//   printOtherEpistasis(fe.Epistasis, "epistatic interaction", ", ");
//   printNoInteractionGenes(fe.genesNoInt);
//   printAllOrderG(fe.allOrderG);
// }
