#ifndef FORCE_DEBUG
#define NDEBUG
#endif



#include <string.h>
#include <unistd.h>
#include "base/CommandLineParser.h"
#include "analysis/SequenceMatch.h"
#include "analysis/GridSearch.h"
#include "analysis/SeqChunk.h"
#include "analysis/MatchDynProg.h"
#include "util/SComm.h"
#include "util/SysTime.h"
#include <math.h>
#include "util/FindProcess.h"
#include <pthread.h>


//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{
  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<string> oStringCmmd("-o","output directory");
  commandArg<string> distanceCmmd("-seeddist","distance between pre-filter seeds (increase for close genomes)", "1");
  commandArg<string> seedWidthCmmd("-filterwidth","width of the seed filter", "2");


  commandLineParser P(argc,argv);
  P.SetDescription("Wrapper around high-sensitivity alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);

  P.registerArg(distanceCmmd);
  P.registerArg(seedWidthCmmd);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string theDistance = P.GetStringValueFor(distanceCmmd);
  string seedwidth = P.GetStringValueFor(seedWidthCmmd);
   
 
 
  int i, j;

  cout << "SATSUMA: Welcome to SatsumaSyntenyLight! Current date and time: " << GetTimeStatic() << endl;

  string probe = output + "/satsuma.log";

  FILE * pProbe = fopen(probe.c_str(), "w");
  if (pProbe == NULL) {
    string mk_cmd = "mkdir " + output;
    cout << "Creating directory " << output << endl;
    system(mk_cmd.c_str());
    pProbe = fopen(probe.c_str(), "w");
    fprintf(pProbe, "Starting\n");
    fclose(pProbe);
   } else {
    fclose(pProbe);
  }


  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;


  if (output == "") {
    cout << "You MUST specify a valid output directory (option '-o')!" << endl;
    return -1;
  }


  char exec_dir[8192];
  strcpy(exec_dir, pExec);
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      break;
    }
  }


  //==================================================================
  string seedFile = output + "/xcorr_aligns.seeds.out";

  string cmd = exec_dir;
  cmd += "FilterGridSeeds -t ";
  cmd += sTarget;
  cmd += " -q ";
  cmd += sQuery;
  cmd += " -o ";
  cmd += seedFile;
  cmd += " -distance ";
  cmd += theDistance;
  cmd += " -w ";
  cmd += seedwidth;
  
  cout << "Running seed pre-filter " << cmd << endl;
  system(cmd.c_str());
  
 //==================================================================

  

  string mergeCmd = exec_dir;
  
  mergeCmd += "/ChainMatches -i " + seedFile;
  mergeCmd += " -o " + output + "/xcorr_aligns.chained.out";
  
  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());


  mergeCmd = exec_dir;

  mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.chained.out";
  mergeCmd += " -q " + sQuery + " -t " + sTarget;
  mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
  mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";
  
  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());
  
 
 

  cout << "SATSUMA: all done, date and time: " << GetTimeStatic() << endl;
  cout << "all done!" << endl;
  

  return 0;
}

