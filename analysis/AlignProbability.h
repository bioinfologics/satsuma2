#ifndef _ALIGNPROBABILITY_H_
#define _ALIGNPROBABILITY_H_



class DNAVector;



double GetMatchProbability(const DNAVector & target, 
			   const DNAVector & query, 
			   int startTarget,
			   int startQuery,
			   int length,
			   double targetSize);


double GetMatchProbabilityEx(double & ident,
			     const DNAVector & target, 
			     const DNAVector & query, 
			     int startTarget,
			     int startQuery,
			     int length,
			     double targetSize);


double GetRawScore(double ident,
		   int len,
		   int tLen,
		   double gcQuery,
		   double gcTarget);

double GetPValueForColaScore(double ident, 
			     double length, 
			     int qLen, 
			     int tLen);



#endif 

