
#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "analysis/AlignProbability.h"
#include "analysis/DNAVector.h"

#include <math.h>

double CDF(double x, double m, double s)
{
  
  double r = 0.5 * (1. + erf((x-m)/s/1.414213562));

  //cout << "x=" << x << " m=" << m << " s=" << s << " cdf=" << r << endl;
  return r;
}

double Sigma(double p, int N)
{
  return sqrt(p * (1. - p) * (double)N);
}


double GCAdjustExpect(double gc, int N, double gc_target)
{
  double at_target = 1. - gc_target;
  
  double r = gc * gc_target;
  
  //cout << "Expect gc matches: " << r << endl;
  //cout << "Expect at matches: " << (double)(N - gc) * at_target << endl; 

  r += (double)(N - gc) * at_target;


  return r / (double)N / 2.;
  
}

double GetMatchProbability(const DNAVector & target, 
			   const DNAVector & query, 
			   int startTarget,
			   int startQuery,
			   int length,
			   double targetSize)
{

  double ident;
  return GetMatchProbabilityEx(ident,
			       target, 
			       query, 
			       startTarget,
			       startQuery,
			       length,
			       targetSize);


}

double GetMatchProbabilityEx(double & ident,
			     const DNAVector & target, 
			     const DNAVector & query, 
			     int startTarget,
			     int startQuery,
			     int length,
			     double targetSize)
{

  int i;
  double gcCountTarget = 0;
  double gcCountQuery = 0;
  double matches = 0;

  int trimLeft = -1;
  int lastMatch = 0;

  for (i=0; i<length; i++) {
    if (target[i+startTarget] == query[i+startQuery]) {
      if (trimLeft == -1)
	trimLeft = i;
      lastMatch = i;
    }
 
    matches += DNA_Equal(target[i+startTarget], query[i+startQuery]);

    gcCountTarget += DNA_C(target[i+startTarget]) + DNA_G(target[i+startTarget]); 
    gcCountQuery += DNA_C(query[i+startQuery]) +  DNA_G(query[i+startQuery]);

    //if (target[i+startTarget] == 'G' || target[i+startTarget] == 'C')
    //  gcCountTarget++;
    //if (query[i+startQuery] == 'G' || query[i+startQuery] == 'C')
    //gcCountQuery++;
  }
  
  int trimRight = (length - 1 - lastMatch);
  //cout << "leftTrim=" << trimLeft << " rightTrim=" << trimRight << endl; 
  //cout << "Adjusting length from " << length;
  

  //length -= trimLeft;
  //length -= trimRight;
  //cout << " to " << length << endl;


  ident = matches/(double)length;
  
  double gc = gcCountQuery;
      
  double p_match = GCAdjustExpect(gc, length, gcCountTarget/(double)length);
     
  //cout << "ident=" << ident << " p_match=" << p_match << " gc_target=" << (double)gcCountTarget/(double)length << endl;


  double s = Sigma(p_match, length);
      
  double m = p_match * (double)length;

  double x = (double)length * ident;
    
  double cdf = CDF(m, x, s);
  double expect = cdf * targetSize; // Twice the genome size if we do RC!!
  double p_val = exp(-expect);
     
  return p_val;
}

double GetRawScore(double ident,
		   int length,
		   int targetSize,
		   double gcQuery,
		   double gcTarget)
{

  int gcCountTarget = (int)(0.5 + (double)length * gcTarget);
  int gcCountQuery = (int)(0.5 + (double)length * gcQuery);

  double gc = (double)gcCountQuery;
      
  double p_match = GCAdjustExpect((int)gc, length, (double)gcCountTarget/(double)length);
     
  //cout << "ident=" << ident << " p_match=" << p_match << " gc_target=" << (double)gcCountTarget/(double)length << endl;


  double s = Sigma(p_match, length);
      
  double m = p_match * (double)length;

  double x = (double)length * ident;
    
  double cdf = CDF(m, x, s);

  //cout << "cdf=" << cdf << endl;

  double expect = cdf * targetSize; // Twice the genome size if we do RC!!
  double p_val = exp(-expect);
     
  //cout << "expect=" << expect << endl;
  return p_val;
}

double GetExpectCola(double ident,
		     int length,
		     double expect_ident,
		     int tries)
{
  double p_match = expect_ident;
    
  double scale = 4.0; // Adjust standard deviation for normal distribution 
  int scaledLength = (double)length/scale;
  if(length<2) { scaledLength = 1; } // Added to prevent 0 sigma for length<2 TODO - review
  double s = expect_ident * scale * Sigma(p_match, scaledLength);

  //cout << "Sigma=" << s << endl; 
      
  double m = p_match * (double)length;
  double x = (double)length * ident;
    
  //cout << "len=" << length << " ident=" << ident << " s=" << s << endl;

  double cdf = CDF(m, x, s); // Cumulative distribution function

  double targetSize = tries;

  double expect = cdf * targetSize; // Adjust for sequence lengths
  
  double p_val = -expm1(-expect); // Poisson distribution

  //cout << "s=" << s << " cdf=" << cdf << " size=" << targetSize << " expect=" << expect << " p=" << p_val << endl; 
  //cout << "expect=" << expect << " p-val=" << cdf << endl;
  return p_val;
}


double GetPValueForColaScore(double ident, double len, int qLen, int tLen) 
{
/*
  double tries = qLen * tLen; // How many times can we try?
  tries *= 2.; // Exeggerate a bit
  if (tries < 0)
    tries = 1;

  double dd = GetExpectCola(ident, len, 0.43, tries);
  //cout << ident << " " << len << " " << dd << " " << tries << endl;
  return dd;
*/
  double lambda = 0.28;
  double mu     = 8.5;
  double x      = ident;
  double p      = 1 - exp(-exp(-lambda*(x-mu)));
  return p;
}

