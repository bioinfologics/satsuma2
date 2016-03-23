#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "src/Papaya/Papaya.h"


int main(int argc,char** argv)
{
  commandArg<string> targetCmmd("-t","target fasta file");
  commandArg<string> queryCmmd("-q","query fasta file");
  commandArg<string> outCmmd("-o","output file");
  commandArg<int> nCmmd("-n","number of processes", 1);
  commandArg<int> numCmmd("-k","number of 12-mers used for seeding", 2);
  commandArg<int> wordCmmd("-w","size of a k-mer", 12);
  commandArg<int> laCmmd("-la","number of lookahead k-mers", 12);
  commandArg<int> skipCmmd("-la_skip","skip k-mers in lookahead", 2);
  commandArg<bool> selfCmmd("-self","suppress matches to itself (query=target)", false);
  commandArg<double> identCmmd("-min","min identity", 0.999);
  commandArg<int> maxGapCmmd("-max_gap","maximum gap", 500);
  commandArg<int> minCmmd("-min_len","minimum alignment length", 30);
  commandArg<bool> repCmmd("-noreps","do not align to repeats (soft-masked)", false);
  commandLineParser P(argc,argv);

  P.SetDescription("Welcome to Papaya/Cola, alignments with seeds (and more)!");

  P.registerArg(targetCmmd);
  P.registerArg(queryCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(numCmmd);
  P.registerArg(laCmmd);
  P.registerArg(skipCmmd);
  P.registerArg(selfCmmd);
  P.registerArg(repCmmd);
  P.registerArg(wordCmmd);
  P.registerArg(maxGapCmmd);
  P.registerArg(minCmmd);
  P.registerArg(identCmmd);

  P.parse();

  string targetFile = P.GetStringValueFor(targetCmmd);
  string queryFile  = P.GetStringValueFor(queryCmmd);
  string outFile    = P.GetStringValueFor(outCmmd);
  int numKmers      = P.GetIntValueFor(numCmmd);
  int lookAhead     = P.GetIntValueFor(laCmmd);
  int lookAheadSkip = P.GetIntValueFor(skipCmmd);
  int nProc         = P.GetIntValueFor(nCmmd);
  int wordSize      = P.GetIntValueFor(wordCmmd);
  bool bSelf        = P.GetBoolValueFor(selfCmmd);
  bool bNoReps      = P.GetBoolValueFor(selfCmmd);
  int maxGap        = P.GetIntValueFor(maxGapCmmd);
  int minLen        = P.GetIntValueFor(minCmmd);
  double minIdent   = P.GetDoubleValueFor(identCmmd);


  PapayaSeedingParams seedingParams(numKmers, lookAhead, lookAheadSkip, 
                                    wordSize, bSelf, bNoReps);
  PapayaAlignerParams alignerParams(maxGap, minLen, minIdent);
  PapayaAligner pAligner;
  pAligner.initialize(targetFile, seedingParams, nProc);
  pAligner.processQuery(queryFile, outFile, alignerParams);

  return 0;

}
  
