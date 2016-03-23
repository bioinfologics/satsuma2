#include <string>
#include "../base/CommandLineParser.h"

bool Run(const string & exec, const string & cmd) 
{

  string ex = exec;
  //if (ex != "")
  //ex += "/";
  ex += cmd;

  cout << "Executing: " << ex << endl;
  int ret = system(ex.c_str());   // TODO: comment 3 lines back in
  if (ret == 0)
    return true;
  cout << "********************************************************" << endl;
  cout << "ERROR!!! See log file for more error messages. ABORTING!" << endl;
  cout << "********************************************************" << endl;
  exit(-1);
  return false;
}


int main( int argc, char** argv )
{

  commandArg<string> targetCmmd("-t","target fasta file (in chromosome coordinates)");
  commandArg<string> queryCmmd("-q","query fasta file (the assembly)");
  commandArg<string> outCmmd("-o","output directory");
  commandArg<int> nCmmd("-n","number of CPUs (for full Satsuma run)", 25);
  commandArg<bool> thoroughCmmd("-thorough","runs a full Satsuma alignment (slow!!)", false);
  commandArg<bool> chrCmmd("-pseudochr","maps scaffolds into chromosomes", false);
  commandArg<bool> satCmmd("-s","run SatsumaSynteny at the end", false);

  commandLineParser P(argc,argv);

  P.SetDescription("Runs a pipeline that scaffolds an assembly by synteny.");
  P.registerArg(targetCmmd);
  P.registerArg(queryCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(thoroughCmmd);
  P.registerArg(chrCmmd);
  P.registerArg(satCmmd);
  
  P.parse();
  
  string target = P.GetStringValueFor(targetCmmd);
  string query = P.GetStringValueFor(queryCmmd);
  string outName = P.GetStringValueFor(outCmmd);

  int nth = P.GetIntValueFor(nCmmd);
  bool bThorough = P.GetBoolValueFor(thoroughCmmd);
  bool bChr = P.GetBoolValueFor(chrCmmd);
  bool bFinal = P.GetBoolValueFor(satCmmd);

  int i;


  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;

  char exec_dir[8192];
  strcpy(exec_dir, pExec);
  bool bTrunc = false;
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      bTrunc = true;
      break;
    }
  }
  if (!bTrunc)
    strcpy(exec_dir, "");
   
  string cmd;
  cmd = "mkdir " + outName;
  system(cmd.c_str());

  if (bThorough) {
    cmd = "SatsumaSynteny -q ";
    cmd += target;
    cmd += " -t ";
    cmd += query;
    cmd += " -n ";
    cmd += Stringify(nth);
    cmd += " -o ";
    cmd += outName;
    cmd += "/xcorr_first";
    Run(exec_dir, cmd);

  } else {
    cmd = "mkdir ";
    cmd += outName;
    cmd += "/xcorr_first";
    system(cmd.c_str());

    // FilterGridSeeds
    cmd = "FilterGridSeeds -q ";
    cmd += target;
    cmd += " -t ";
    cmd += query;
    cmd += " -o ";
    cmd += outName;
    cmd += "/xcorr_first/filterseeds.dat";

    Run(exec_dir, cmd);

    // ChainMatches
    cmd = "ChainMatches -i ";
    cmd += outName;
    cmd += "/xcorr_first/filterseeds.dat -o ";
    cmd += outName;
    cmd += "/xcorr_first/filterseeds_chained.dat";
     
    Run(exec_dir, cmd);
    cmd = "MergeXCorrMatches -q ";
    cmd += target;
    cmd += " -t ";
    cmd += query;
    cmd += " -i ";
    cmd += outName;
    cmd += "/xcorr_first/filterseeds_chained.dat";
    cmd += " -o ";
    cmd += outName;
    cmd += "/xcorr_first/satsuma_summary.chained.out";

    Run(exec_dir, cmd);    

  }
  
  cmd = "SortSatsuma -i ";
  cmd += outName;
  cmd += "/xcorr_first/satsuma_summary.chained.out > ";
  cmd += outName;
  cmd += "/xcorr_first/satsuma_summary.chained.out_sorted";
  Run(exec_dir, cmd);    


  cmd = "MergeScaffoldsBySynteny -i ";
  cmd += outName;
  cmd += "/xcorr_first/satsuma_summary.chained.out_sorted -f ";
  cmd += query;
  cmd += " -o ";
  cmd += outName;
  cmd += "/superscaffolds.fasta";
  Run(exec_dir, cmd);    
    
  // Change query name
  query = outName;
  query += "/superscaffolds.fasta";
 

  //=======================================================
  if (bChr) {

    if (bThorough) {
      cmd = "SatsumaSynteny -q ";
      cmd += target;
      cmd += " -t ";
      cmd += query;
      cmd += " -n ";
      cmd += Stringify(nth);
      cmd += " -o ";
      cmd += outName;
      cmd += "/xcorr_second";
      Run(exec_dir, cmd);
    
    } else {
      
      cmd = "mkdir ";
      cmd += outName;
      cmd += "/xcorr_second";
      //Run("", cmd);
      system(cmd.c_str());

      // FilterGridSeeds
      cmd = "FilterGridSeeds -q ";
      cmd += target;
      cmd += " -t ";
      cmd += query;
      cmd += " -o ";
      cmd += outName;
      cmd += "/xcorr_second/filterseeds.dat";
      
      Run(exec_dir, cmd);
      
      // ChainMatches
      cmd = "ChainMatches -i ";
      cmd += outName;
      cmd += "/xcorr_second/filterseeds.dat -o ";
      cmd += outName;
      cmd += "/xcorr_second/filterseeds_chained.dat";
      
      Run(exec_dir, cmd);
      cmd = "MergeXCorrMatches -q ";
      cmd += target;
      cmd += " -t ";
      cmd += query;
      cmd += " -i ";
      cmd += outName;
      cmd += "/xcorr_second/filterseeds_chained.dat";
      cmd += " -o ";
      cmd += outName;
      cmd += "/xcorr_second/satsuma_summary.chained.out";
      
      Run(exec_dir, cmd);    
      
    }
    
    cmd = "SortSatsuma -i ";
    cmd += outName;
    cmd += "/xcorr_second/satsuma_summary.chained.out > ";
    cmd += outName;
    cmd += "/xcorr_second/satsuma_summary.chained.out_sorted";
    Run(exec_dir, cmd);    
    
    
    cmd = "OrderOrientBySynteny -i ";
    cmd += outName;
    cmd += "/xcorr_second/satsuma_summary.chained.out_sorted -f ";
    cmd += query;
    cmd += " -o ";
    cmd += outName;
    cmd += "/pseudochromosomes.fasta";
    Run(exec_dir, cmd);    
    
    // Change query name
    query = outName;
    query += "/pseudochromosomes.fasta";
    

  }

  //=======================================================
  
  if (bFinal) {
    cmd = "SatsumaSynteny -t ";
    cmd += target;
    cmd += " -q ";
    cmd += query;
    cmd += " -n ";
    cmd += Stringify(nth);
    cmd += " -o ";
    cmd += outName;
    cmd += "/xcorr_final";
    Run(exec_dir, cmd);
  }

  return 0;
}

