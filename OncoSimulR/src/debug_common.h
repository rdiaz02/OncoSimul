#ifndef _DEBUG_COMMON_H__
#define _DEBUG_COMMON_H__

#include<Rcpp.h>
#define DEBUGW

// #define DEBUGV
// DEBUGV2 leads to very verbose debugging output
//#define DEBUGV2
#define DP {Rcpp::Rcout <<"\nHola";		\
    Rcpp::Rcout << "\n A test of DPtest";	\
    Rcpp::Rcout << "     tSample  = " << tSample << "\n";} 
#define DP1(x) {Rcpp::Rcout << "\n DEBUG2: I am at " << x << std::endl;}
#define DP2(x) {Rcpp::Rcout << "\n DEBUG2: Value of " << #x << " = " << x << std::endl;}
#define DP3(x, t){                                                    \
    Rcpp::Rcout <<"\n DEBUG2:" ;                                      \
    for(int xut = 0; xut < t; ++xut) Rcpp::Rcout << "\t ";            \
    Rcpp::Rcout << "  I am at " << x << std::endl;}
#define DP4(x, t){                                                    \
    Rcpp::Rcout <<"\n DEBUG2:" ;                                      \
    for(int xut = 0; xut < t; ++xut) Rcpp::Rcout << "\t ";            \
    Rcpp::Rcout << "  Value of " << #x << " = " << x << std::endl; }

// // Windows compiler in BioC is pre 4.8.0, so no to_string
// // From http://stackoverflow.com/a/5590404
// #define SSTR( x ) dynamic_cast< std::ostringstream & >( 
//        ( std::ostringstream() << std::dec << x ) ).str()

#define DP_ti_notfinite {if(std::abs(ti - ti2) > 1e-5) { \
	  DP2(ti); \
	  DP2(ti2); \
	  DP2(numerator); \
	  DP2(numerator2); \
	  DP2(denominator); \
	  DP2(denominator2); \
	  DP2(invspr); \
	  DP2(r); \
	  DP2(r1); \
	  DP2(tmp); \
	  DP2(tmp2); \
	  DP2(oneminusr); \
	} \
	Rcpp::Rcout << "\n r1 = " << r1; \
	Rcpp::Rcout << "\n r = " << r; \
	Rcpp::Rcout << "\n numerator = " << numerator; \
	Rcpp::Rcout << "\n denominator = " << denominator; \
	Rcpp::Rcout << "\n ti2 = " << ti2; \
	Rcpp::Rcout << "\n numerator2 = " << numerator2; \
	Rcpp::Rcout << "\n denominator2 = " << denominator2; \
	Rcpp::Rcout << "\n is r > 1? " << (r > 1.0) << "\n"; \
	Rcpp::Rcout << "\n is r < 0? " << (r < 0.0) << "\n"; \
	Rcpp::Rcout << "\n is eq12 < r? " << (eq12 < r) << "\n"; \
	Rcpp::Rcout << "\n tmp = " << tmp << "\n"; \
	Rcpp::Rcout << "\n tmp2 = " << tmp2 << "\n"; \
	Rcpp::Rcout << "\n eq12 = " << eq12 << "\n"; \
	print_spP(spP); \
}


#ifdef DEBUGW
#define ASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << "\n\n";		\
    }									\
  }
#else
#define ASSERT(x) do {} while(0)
#endif


#ifdef DEBUGW
#define STOPASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << std::endl;		\
	throw std::out_of_range("STOPASSERT");				\
    }									\
  }
#else
#define STOPASSERT(x) do {} while(0)
#endif


#endif

// If you want VERY verbose debugging output



#ifdef DEBUGV2
#define DEBUGfs {Rcpp::Rcout << " DEBUGV\n";	       \
    Rcpp::Rcout << "\n ForceSample? " << forceSample   \
 		<< "  tSample " << tSample	       \
 		<< "  currentTime " << currentTime; }
// We always warn about this, since interaction with ti==0
#define DEBUGfs2 {
	Rcpp::Rcout << "\n Forced sampling triggered for next loop: \n    " \
	<< " popParams[nextMutant].popSize = " \
	<< popParams[nextMutant].popSize \
	<< " > ratioForce * detectionSize \n"; \
	Rcpp::Rcout << " when nextMutant = " \
	<< nextMutant \
	<< " at iteration " << iter << "\n"; \
}
#define DEBUG_detect_duplicates(x, y) {detect_ti_duplicates(mapTimes, x, y); \
  }
#define DEBUG_52(y, x, z) {Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, " \
				       << z				\
				       << "\n     tSample  = " << tSample \
				       << "\n\n**   Species  = " << x	\
				       << "\n       genotype =  ";	\
    print_Genotype(Genotypes[x]);					\
    Rcpp::Rcout << "\n       popSize = " << popParams[x].popSize	\
		<< "\n       currentTime = " << currentTime		\
		<< "\n       popParams[i].nextMutationTime = "		\
		<< y							\
		<< " \n     species R " << popParams[x].R		\
		<< " \n     species W " << popParams[x].W		\
		<< " \n     species death " << popParams[x].death	\
		<< " \n     species birth " << popParams[x].birth;	\
  }
#define DEBUGfsnl {Rcpp::Rcout << "\n Forced sampling triggered for next loop: \n    " \
			       <<  " popParams[nextMutant].popSize = "	\
			       <<  popParams[nextMutant].popSize	\
			       << " > ratioForce * detectionSize \n";	\
    Rcpp::Rcout << " when nextMutant = " << nextMutant			\
		<<" at iteration " << iter << "\n";			\
}
#define DEBUG_1456 {if( (currentTime - popParams[nextMutant].timeLastUpdate) < 0.0) { \
  DP1("currentTime was set to minNextMutationTime above");	  \
  DP2(currentTime);						  \
  DP2(minNextMutationTime);					  \
  DP2(tSample);							  \
  DP2(popParams[nextMutant].timeLastUpdate);				\
  DP2( (currentTime -  popParams[nextMutant].timeLastUpdate) );		\
  DP2( (currentTime <  popParams[nextMutant].timeLastUpdate) );		\
  DP2( (currentTime ==  popParams[nextMutant].timeLastUpdate) );	\
  DP2(nextMutant);							\
  DP2(u_1);								\
  DP2(u_2);								\
  DP2(tmpdouble1);							\
  DP2(tmpdouble2);							\
  DP2(popParams[nextMutant].timeLastUpdate);				\
  DP2(popParams[u_1].timeLastUpdate);					\
  DP2(popParams[u_2].timeLastUpdate);					\
  DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) ); \
  DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) ); \
  DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) ); \
  print_spP(popParams[nextMutant]);					\
  throw std::out_of_range("new species: currentTime - timeLastUpdate[sp] out of range. ***###!!!Serious bug!!!###***"); \
 }									\
}
 // Yes, the difference could be 0 if two next mutation times are identical.
// You enable detect_ti_duplicates and use trigger-duplicated-ti.R
// to see it.
// Often the involved culprits (nextMutant and the other, say sp)
// were lastUpdated with tiny difference and they were, when updated
// given an identical ti, each in its own run.
// Key is not timeLastUpdate. This is a possible sequence of events:
//    - at time t0, species that will become nextMutant is updated and gets ti = tinm
//    - t1: species u1 gets ti = tinm
//    - t2: species u2 gets some ti > tinm
//    - tinm becomes minimal, so we mutate u1, and it mutates to u2
//    - (so now the timeLastUpdate of u1 = u2 = tinm)
//    - nextMutant is now mutated, and it mutates to u2, which becomes sp
//    - tinm = timeLastUpdate of u1 and u2.
//    - You will also see that number of mutations, or genotypes are such
//      that, in this case, u2 is the most mutated, etc.
//    - If you enable the detect_ti_duplicates, you would have seen duplicated ti
//      for nextMutant and u1

//   Even simpler is if above, nextMutant will mutate to u1 (not u2) so u1 becomes sp.
//next is set to minNextMutationTime above
#define DEBUG_1536 { if( (currentTime - popParams[sp].timeLastUpdate) < 0.0) { \
      DP2(currentTime);							\
      DP2(minNextMutationTime);						\
      DP2(tSample);							\
      DP2(popParams[sp].timeLastUpdate);				\
      DP2( (currentTime -  popParams[sp].timeLastUpdate) );		\
      DP2( (currentTime <  popParams[sp].timeLastUpdate) );		\
      DP2( (currentTime ==  popParams[sp].timeLastUpdate) );		\
      DP2(sp);								\
      DP2(nextMutant);							\
      DP2(u_1);								\
      DP2(u_2);								\
      DP2(tmpdouble1);							\
      DP2(tmpdouble2);							\
      DP2(popParams[sp].timeLastUpdate);				\
      DP2(popParams[nextMutant].timeLastUpdate);			\
      DP2(popParams[u_1].timeLastUpdate);				\
      DP2(popParams[u_2].timeLastUpdate);				\
      DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) ); \
      DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) ); \
      DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) ); \
      print_spP(popParams[sp]);						\
      print_spP(popParams[nextMutant]);					\
      throw std::out_of_range("currentTime - timeLastUpdate[sp] out of range.  ***###!!!Serious bug!!!###***"); \
    }									\
    if( (currentTime - popParams[nextMutant].timeLastUpdate) < 0.0) {	\
      DP2(currentTime);							\
      DP2(minNextMutationTime);						\
      DP2(tSample);							\
      DP2(popParams[nextMutant].timeLastUpdate);			\
      DP2( (currentTime -  popParams[nextMutant].timeLastUpdate) );	\
      DP2( (currentTime <  popParams[nextMutant].timeLastUpdate) );	\
      DP2( (currentTime ==  popParams[nextMutant].timeLastUpdate) );	\
      DP2(sp);								\
      DP2(nextMutant);							\
      DP2(u_1);								\
      DP2(u_2);								\
      DP2(tmpdouble1);							\
      DP2(tmpdouble2);							\
      DP2(popParams[sp].timeLastUpdate);				\
      DP2(popParams[nextMutant].timeLastUpdate);			\
      DP2(popParams[u_1].timeLastUpdate);				\
      DP2(popParams[u_2].timeLastUpdate);				\
      DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) ); \
      DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) ); \
      DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) ); \
      print_spP(popParams[sp]);						\
      print_spP(popParams[nextMutant]);					\
      throw std::out_of_range("currentTime - timeLastUpdate[nextMutant] out of range. ***###!!!Serious bug!!!###***"); \
    }									\
    Rcpp::Rcout <<"\n     Mutated to existing species " << sp		\
		<< " (Genotype = ";					\
    print_Genotype(Genotypes[sp]);					\
    Rcpp::Rcout <<  ")"                                                 \
		<< "\n from species "  <<   nextMutant			\
		<< " (Genotypes = ";					\
    print_Genotype(Genotypes[nextMutant]);				\
    Rcpp::Rcout	<< ")";							\
  }
#define DEBUGrrr { Rcpp::Rcout << "\n reachDetection = " << reachDetection; \
    Rcpp::Rcout << "\n forceRerun =  " << forceRerun  << "\n";		\
  }
#define DEBUG_nr {Rcpp::Rcout << "\n\n     ********* 5.9 ******\n " \
			      << "     Species  = " << i	    \
			      << "\n      Genotype = ";		    \
    print_Genotype(Genotypes[i]);				    \
    Rcpp::Rcout << "\n      pre-update popSize = "		    \
		<< popParams[i].popSize				    \
		<< "\n      time of sample = " << tSample	    \
		<< "\n      popParams[i].timeLastUpdate = "	    \
		<< popParams[i].timeLastUpdate			    \
		<< ";\n     t for Algo2 = "			    \
		<< tSample - popParams[i].timeLastUpdate	    \
		<< " \n     species R " << popParams[i].R	    \
		<< " \n     species W " << popParams[i].W	    \
		<< " \n     species death " << popParams[i].death   \
		<< " \n     species birth " << popParams[i].birth;  \
}
#define DEBUG_nr2 {Rcpp::Rcout << "\n\n     Removing species i = " << i \
			       << " with genotype = ";			\
    print_Genotype(Genotypes[i]);					\
  }
#define DEBUG_nr3 {Rcpp::Rcout << "\n\n   post-update popSize = "	\
		<< popParams[i].popSize << "\n";			\
}
#define DEBUG_ti {Rcpp::Rcout << "\n\n\n WARNING: ti == 0. Setting it to DBL_MIN \n"; \
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) /  (spP.R + spP.W - 2.0 * spP.birth) , spP.p \opSize);
	Rcpp::Rcout << "\n tmp2 = " << tmp2; \
	Rcpp::Rcout << "\n tmp = " << tmp; \
	Rcpp::Rcout << "\n invspr = " << invspr; \
	Rcpp::Rcout << "\n invpop = " << invpop; \
	Rcpp::Rcout << "\n numerator = " << numerator; \
	Rcpp::Rcout << "\n denominator = " << denominator; \
	Rcpp::Rcout << "\n numerator == denominator? " << (numerator == denominator); \
	Rcpp::Rcout << "\n is r > 1? " << (r > 1.0); \
	Rcpp::Rcout << "\n is r < 0? " << (r < 0.0); \
	Rcpp::Rcout << "\n r = " << r; \
	Rcpp::Rcout << "\n is r == 1? " << (r == 1.0L); \
	Rcpp::Rcout << "\n oneminusr = " << oneminusr; \
	Rcpp::Rcout << "\n is oneminusr == 0? " << (oneminusr == 0.0L); \
	Rcpp::Rcout << "\n r1 = " << r1; \
	Rcpp::Rcout << "\n is r1 == 1? " << (r1 == 1.0); \
	Rcpp::Rcout << "\n is eq12 < r? " << (eq12 < r); \
}
#else
#define DEBUGfs {} do {} while(0)
#define DEBUGfs2 {} do {} while(0)
#define DEBUG_detect_duplicates(x, y) do {} while(0)
#define DEBUG_52(y, x, z) do {} while(0)
#define DEBUGfsnl {} do {} while(0)
#define DEBUG_1456 {} do {} while(0)
#define DEBUG_1536 {} do {} while(0)
#define DEBUG_rrr {} do {} while(0)
#define DEBUG_nr {} do {} while(0)
#define DEBUG_nr2 {} do {} while(0)
#define DEBUG_nr3 {} do {} while(0)
#define DEBUG_ti {} do {} while(0)
#endif

