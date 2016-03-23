#include <string>
#include "base/CommandLineParser.h"
//#include "base/FileParser.h"


void Run(const string & cmd) {
  cout << "Running: " << cmd << endl;
  int ret = system(cmd.c_str());
  if (ret < 0) {
    cout << "ERROR: program " << endl;
    cout << cmd << endl;
    cout << "died with exit code " << ret << endl;
    exit(-1);
  }

}


int main( int argc, char** argv )
{

  commandArg<string> tCmmd("-t","target fasta file");
  commandArg<string> qCmmd("-q","query fasta file");
  commandArg<string> oCmmd("-o","output directory");
  commandArg<int> nCmmd("-n","number of processes");
  commandArg<bool> localCmmd("-local","prepare file for local run");
  commandLineParser P(argc,argv);
  P.SetDescription("Alignment of closely related genomes.");
  P.registerArg(tCmmd);
  P.registerArg(qCmmd);
  P.registerArg(oCmmd);
  P.registerArg(nCmmd);
  P.registerArg(localCmmd);
  
  P.parse();
  
  string tName = P.GetStringValueFor(tCmmd);
  string qName = P.GetStringValueFor(qCmmd);
  string oName = P.GetStringValueFor(oCmmd);
  int n = P.GetIntValueFor(nCmmd);
  bool local = P.GetBoolValueFor(localCmmd);
  
  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;

  int i;

  char exec_dir[8192];
  strcpy(exec_dir, pExec);
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      break;
    }
  }
  string exec = exec_dir;

  string mk = "mkdir ";
  mk += oName;
  int s = system(mk.c_str());
  
  string qt = " -t ";
  qt += tName;
  qt += " -q ";
  qt += qName;
  qt += " ";

  string cmmd = exec + "FilterGridSeeds ";
  cmmd += qt;
  cmmd += " -o ";
  cmmd += oName;
  cmmd += "/xcorr_aligns.seeds.out";
  
  Run(cmmd);

  cmmd = exec + "ChainMatches -i ";
  cmmd += oName;
  cmmd += "/xcorr_aligns.seeds.out -o ";
  cmmd += oName;
  cmmd += "/xcorr_aligns.chained.out";
  
  Run(cmmd);

  cmmd = exec + "MergeXCorrMatches -i ";
  cmmd += oName;
  cmmd += "/xcorr_aligns.chained.out -o ";
  cmmd += oName;
  cmmd += "/chained_summary.out ";
  cmmd += qt;
  cmmd += " > ";
  cmmd += oName;
  cmmd += "/MergeXCorrMatches.out";
  
  Run(cmmd);


  cmmd = exec + "PrepareForDetAligns -i ";
  cmmd += oName;
  cmmd += "/chained_summary.out ";
  cmmd += " > ";
  cmmd += oName;
  cmmd += "/chained_merged.out ";
  
  Run(cmmd);

	 
  cmmd = exec + "SplitFileChunk -i ";
  cmmd += oName;
  cmmd += "/chained_merged.out -n ";
  char tmp[256];
  sprintf(tmp, "%d", n);
  cmmd += tmp;
  
  Run(cmmd);

  string out = oName;
  out += "/detailed_commands";
  FILE * pOut = fopen(out.c_str(), "w");

  for (i=0; i<n; i++) {    
    cmmd = exec + "DPAlignSatsuma";
    cmmd += qt;
    cmmd += " -s ";
    cmmd += oName;
    sprintf(tmp, "%d", i);
    cmmd += "/chained_merged.out";
    cmmd += ".";
    cmmd += tmp;
    cmmd += " -o ";
    cmmd += oName;
    cmmd += "/satsuma_summary.";
    cmmd += tmp;
    cmmd += ".out";
    cmmd += " > ";
    cmmd += oName;
    cmmd += "/satsuma_aligns.";
    cmmd += tmp;
    cmmd += ".out";

    if (local)
      cmmd += " &";

    fprintf(pOut, "%s\n", cmmd.c_str());
  }
  fclose(pOut);

  cout << "All done. The list of commands is found in " << oName << "/detailed_commands" << endl;
 
  return 0;
}
