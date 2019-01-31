#include <string>
#include "../base/CommandLineParser.h"
#include <stdlib.h>

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

std::string getEnvVar( std::string const & key )
{
  char * val = getenv( key.c_str() );
  return val == NULL ? std::string("") : std::string(val);
}

int main( int argc, char** argv )
{

  commandArg<string> targetCmmd("-t","target fasta file (in chromosome coordinates)");
  commandArg<string> queryCmmd("-q","query fasta file (the assembly)");
  commandArg<string> outCmmd("-o","output directory");
  commandArg<int> nCmmd("-n","number of CPUs (for full Satsuma run)", 25);
  commandArg<bool> thoroughCmmd("-thorough","runs a full Satsuma alignment (slow!!)", false);
  commandArg<bool> chrCmmd("-nopseudochr","skips mapping scaffolds into chromosomes", false);
  commandArg<bool> mergeCmmd("-nomerge","skips merging scaffolds", false);
  commandArg<bool> satCmmd("-s","run SatsumaSynteny at the end", false);

  commandLineParser P(argc,argv);

  P.SetDescription("Runs a pipeline that scaffolds an assembly by synteny.");
  P.registerArg(targetCmmd);
  P.registerArg(queryCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(thoroughCmmd);
  P.registerArg(chrCmmd);
  P.registerArg(mergeCmmd);
  P.registerArg(satCmmd);
  
  P.parse();
  
  string target = P.GetStringValueFor(targetCmmd);
  string query = P.GetStringValueFor(queryCmmd);
  string outName = P.GetStringValueFor(outCmmd);

  int nth = P.GetIntValueFor(nCmmd);
  bool bThorough = P.GetBoolValueFor(thoroughCmmd);
  bool bNoChr = P.GetBoolValueFor(chrCmmd);
  bool bNoMerge = P.GetBoolValueFor(mergeCmmd);
  bool bFinal = P.GetBoolValueFor(satCmmd);

  int i;

  string satsuma2_path  = getEnvVar("SATSUMA2_PATH");  
  string current_path = getEnvVar("PWD");

  if (satsuma2_path==""){
    cout << "ERROR: SATSUMA2_PATH variable not set, please set it to the binary path." <<endl;
    return -1;
  }
  cout<< "Path for Satsuma2: '"<<satsuma2_path<<"'"<<endl;


  string exec_dir = satsuma2_path + "/";

   
  string cmd;
  cmd = "mkdir " + outName;
  system(cmd.c_str());

  if (!bNoMerge) {
    if (bThorough) {
      cmd = "SatsumaSynteny2 -q ";
      cmd += target;
      cmd += " -t ";
      cmd += query;
      cmd += " -slaves ";
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
  }

  //=======================================================
  if (!bNoChr) {

    if (bThorough) {
      cmd = "SatsumaSynteny2 -q ";
      cmd += target;
      cmd += " -t ";
      cmd += query;
      cmd += " -slaves ";
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
    cmd = "SatsumaSynteny2 -t ";
    cmd += target;
    cmd += " -q ";
    cmd += query;
    cmd += " -slaves ";
    cmd += Stringify(nth);
    cmd += " -o ";
    cmd += outName;
    cmd += "/xcorr_final";
    Run(exec_dir, cmd);
  }

  return 0;
}

