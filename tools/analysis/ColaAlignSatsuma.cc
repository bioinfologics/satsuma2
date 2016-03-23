#include <string>
#include <unistd.h>
#include "../base/CommandLineParser.h"
#include "../src/Cola/NSGAaligner.h"
#include "RefineSatsuma.h"


int main(int argc,char** argv)
{
  
  commandArg<string> aStringCmmd("-q","query sequence");
  commandArg<string> bStringCmmd("-t","target sequence");
  commandArg<string> satsumaCmmd("-s","satsuma summary file");
  commandArg<int> alignerTypeCmd("-a","Aligner type - Choose 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW ", 1);
  commandArg<int> outputModeCmmd("-m","Output mode - 0:Full  1:Summary CSV  2:MultiAlignFormat", 2);
  commandArg<string> realignOutCmmd("-r","realignment output file", "reAlign.out");
  commandArg<string> gapOutCmmd("-g","gap output file", "gapAlign.out");
  commandArg<string> bothOutCmmd("-o","combined output file", "");
  commandArg<int> mBoundCmmd("-b","Merge Boundary for merging consecutive blocks", 5);
  commandArg<int> threadCmmd("-T","Number of cores to run with", 2);

  commandLineParser P(argc,argv);
  P.SetDescription("Realigns global alignment in satsuma format with Cola.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(satsumaCmmd);
  P.registerArg(alignerTypeCmd);
  P.registerArg(outputModeCmmd);
  P.registerArg(bothOutCmmd);
  P.registerArg(realignOutCmmd);
  P.registerArg(gapOutCmmd);
  P.registerArg(mBoundCmmd);
  P.registerArg(threadCmmd);

  P.parse();

  string aString        = P.GetStringValueFor(aStringCmmd);
  string bString        = P.GetStringValueFor(bStringCmmd);
  string satsumaFile    = P.GetStringValueFor(satsumaCmmd);
  AlignerType aType     = AlignerType(P.GetIntValueFor(alignerTypeCmd));
  int outputMode        = P.GetIntValueFor(outputModeCmmd); 
  string realignOutFile = P.GetStringValueFor(realignOutCmmd);
  string gapOutFile     = P.GetStringValueFor(gapOutCmmd);
  string bothOutFile    = P.GetStringValueFor(bothOutCmmd);
  int  mergeBoundary    = P.GetIntValueFor(mBoundCmmd); 
  int  numOfThreads     = P.GetIntValueFor(threadCmmd);

  vecDNAVector query, target;
  query.Read(aString);
  target.Read(bString);

  if (bothOutFile == "") {
    RefineSatsuma refinementAligner(target, query,
				    satsumaFile, realignOutFile,
				    gapOutFile, outputMode);
    refinementAligner.setAligner(AlignerParams(-1, aType));
    refinementAligner.alignAll(mergeBoundary, numOfThreads);
    string catCmd1 = "cat " + gapOutFile + ".* >"  + gapOutFile;
    string catCmd2 = "cat "  + realignOutFile + ".* > " + realignOutFile;
    string rmCmd = "rm " + gapOutFile + ".* " + realignOutFile + ".*";
    int ret = system(catCmd1.c_str());
    ret = system(catCmd2.c_str());
    ret = system(rmCmd.c_str());
 
  } else {
    RefineSatsuma refinementAligner(target, query,
				    satsumaFile, bothOutFile,
				    bothOutFile, outputMode);

    refinementAligner.setAligner(AlignerParams(-1, aType));
    refinementAligner.alignAll(mergeBoundary, numOfThreads);
    string catCmd1 = "cat " + bothOutFile + ".* >"  + bothOutFile;
    string rmCmd = "rm " + bothOutFile + ".*";
    int ret = system(catCmd1.c_str());
    ret = system(rmCmd.c_str());
    
  }
  return 0;
}
