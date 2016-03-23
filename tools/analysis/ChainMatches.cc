#ifndef FORCE_DEBUG
#define NDEBUG
#endif



#include <string>

#include "base/CommandLineParser.h"
#include "analysis/MatchDynProg.h"
#include "analysis/SequenceMatch.h"

//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{
  //commandArg<string> aStringCmmd("-q","query fasta sequence");
  //commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<string> iStringCmmd("-i","XCorr input file (binary)");
  commandArg<string> oStringCmmd("-o","output file (binary)");
  commandArg<int> mLenCmmd("-min_len","minimum length to force keep", 220);
  commandArg<double> mIdCmmd("-min_id","minimum identity to force keep", 0.57);
  commandArg<bool> dupCmmd("-dups","allow for duplications", false);
  

  commandLineParser P(argc,argv);
  P.SetDescription("Chains matches found by HomologyByXCorr via dynamic programming");
  //P.registerArg(aStringCmmd);
  //P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(iStringCmmd);
  P.registerArg(mLenCmmd);
  P.registerArg(mIdCmmd);
  P.registerArg(dupCmmd);
  //P.registerArg(lQoffCmmd);
  //P.registerArg(lToffCmmd);

  P.parse();

  // String sQuery = P.GetStringValueFor(aStringCmmd);
  //String sTarget = P.GetStringValueFor(bStringCmmd);
  string input = P.GetStringValueFor(iStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);

  int minLenForce = P.GetIntValueFor(mLenCmmd);
  double minIdentForce = P.GetDoubleValueFor(mIdCmmd);
  bool bDup = P.GetBoolValueFor(dupCmmd);

  //int t_offset = P.GetIntValueFor(lToffCmmd);
  //int q_offset = P.GetIntValueFor(lQoffCmmd);
  


  int i, j, k;

  MultiMatches in, out;

  in.Read(input);
  cout << "Done loading." << endl;

  cout << "Total matches: " << in.GetMatchCount() << endl;

  in.Sort();
  in.Collapse();

  if (bDup)
    RunMatchDynProgMult(out, in);
  else
    RunMatchDynProg(out, in);
  
  //out.Collapse();
  cout << "Saving result..." << endl;
  out.Write(output);
  return 0;
}

