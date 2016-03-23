#include <string>
#include <stdio.h>
#include "base/CommandLineParser.h"
#include "analysis/SatsumaAlign.h"
#include <unistd.h>


class SubmitSeeds
{
public:
  SubmitSeeds(const string & e) {
    m_exec = e;
  }


  void Submit(SAlignVec & res, const string & params, const string & outdir, int n, bool bNoSubmit);

private:
  string m_exec;
  vec<string> m_results;
};


void SubmitSeeds::Submit(SAlignVec & res, const string & params, const string & outdir, int n, bool bNoSubmit)
{
  int i;
  
  m_results.clear();

  for (i=0; i<n; i++) {
    string cmd = m_exec + "PapayaSeed " + params + " ";
    string out = outdir + "/" + "seeds."; 
    string log = outdir + "/" + "PapayaSeed."; 
    char tmp[128];
    sprintf(tmp, "%d", i);
    out += tmp;
    out += ".out";

    log += tmp;
    log += ".log";


    sprintf(tmp, " -block %d -n_block %d ", i, n); 

    m_results.push_back(out);
    cmd += tmp;
    cmd += " -o " + out;
    cmd += " > " + log + " &";

    if (!bNoSubmit) {
      cout << "Submitting " << cmd << endl;
      system(cmd.c_str());
    } else {
      cout << "NOT Submitting " << cmd << endl;
    }
  }

  while (m_results.isize() > 0) {
    sleep(1);
    for (i=0; i<m_results.isize(); i++) {
      string doneFile = m_results[i] + ".done";
      FILE * pTest = fopen(doneFile.c_str(), "r");
      if (pTest != NULL) {
	fclose(pTest);
	sleep(1);
	cout << "Reading result: " << m_results[i] << endl;
	res.Read(m_results[i]);
	cout << "done!" << endl;
	//if (i < m_results.isize()-1)
	m_results[i] = m_results[m_results.isize()-1];
	m_results.resize(m_results.isize()-1);
      }
    }
  }

  //cout << "done here." << endl;
}

int main(int argc,char** argv)
{

  char execPath[512];
  strcpy(execPath, argv[0]);
  cout << argv[0] << endl;
  execPath[strlen(execPath)-6] = 0;
  string exec = execPath;
  cout << "Path to executables: " << exec << endl;

  commandArg<string> targetCmmd("-t","target fasta file");
  commandArg<string> queryCmmd("-q","query fasta file");
  commandArg<string> outCmmd("-o","output directory");
  commandArg<int> nCmmd("-n","number of processes", 1);
  commandArg<int> numCmmd("-k","number of 12-mers used for seeding", 2);
  commandArg<int> wordCmmd("-w","size of a k-mer", 12);
  commandArg<int> laCmmd("-la","number of lookahead k-mers", 12);
  commandArg<int> skipCmmd("-la_skip","skip k-mers in lookahead", 2);
  commandArg<bool> seedCmmd("-no_seed","run PapayaSeed", false);
  commandArg<bool> verboseCmmd("-no_verbose","DO NOT print detailed alignments", false);
  commandArg<bool> selfCmmd("-self","suppress matches to itself (query=target)", false);
  commandArg<bool> repCmmd("-noreps","do not align to repeats (soft-masked)", false);
  commandLineParser P(argc,argv);

  P.SetDescription("Welcome to Papaya, alignments with seeds!");

  P.registerArg(targetCmmd);
  P.registerArg(queryCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(numCmmd);
  P.registerArg(laCmmd);
  P.registerArg(skipCmmd);
  P.registerArg(seedCmmd);
  P.registerArg(verboseCmmd);
  P.registerArg(selfCmmd);
  P.registerArg(repCmmd);
  P.registerArg(wordCmmd);

  P.parse();

  string targetName = P.GetStringValueFor(targetCmmd);
  string queryName = P.GetStringValueFor(queryCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  int num_kmers = P.GetIntValueFor(numCmmd);
  int look_ahead = P.GetIntValueFor(laCmmd);
  int look_ahead_skip = P.GetIntValueFor(skipCmmd);
  int nProc = P.GetIntValueFor(nCmmd);
  int wordSize = P.GetIntValueFor(wordCmmd);
  bool bRunSeed = P.GetBoolValueFor(seedCmmd);
  bool bVerbose = P.GetBoolValueFor(verboseCmmd);
  bool bSelf = P.GetBoolValueFor(selfCmmd);
  bool bNoReps = P.GetBoolValueFor(selfCmmd);
    
  string cmd = "mkdir " + outName;
  system(cmd.c_str());

  cmd = " -t " + targetName + " -q " + queryName + " -k ";


  char tmp[512];
  sprintf(tmp, "%d -la %d -la_skip %d -w %d ", num_kmers, look_ahead, look_ahead_skip, wordSize);
  cmd += tmp;

  if (bSelf)
    cmd += " -self ";
  if (bNoReps)
    cmd += " -noreps ";

  SubmitSeeds ss(exec);

  SAlignVec result;

  ss.Submit(result, cmd, outName, nProc, bRunSeed);

  string outFile = outName + "/" + "allseeds.out";

  cout << "Writing result to " << outFile << endl;
  
  result.Write(outFile);

  string outFileDetails = outName + "/" + "papaya.aligns";
  string outFileSummary = outName + "/" + "papaya_summary.out";
 
  cmd = exec + "PapayaDetails -t " + targetName + " -q " + queryName;
  cmd += " -s " + outFile + " -o " + outFileDetails;
  cmd += " -summary " + outFileSummary;

  if (!bVerbose) {
    cmd += " -verbose";
    cout << "Running " << cmd << endl;
    system(cmd.c_str());
  }

  return 0;

}
  
