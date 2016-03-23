//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif



#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "base/CommandLineParser.h"
#include "analysis/SequenceMatch.h"


void BuildList(svec<string> & xcCommands, 
	       svec<string> & xcOutput, 
	       const string & simple, 
	       const string & outDir,
	       int nBlocks, 
	       bool bLSF,
	       const string & exec_dir)
{

  int i;
  char blocks[64];
  sprintf(blocks, "%d", nBlocks);


  //cout << "outdir: " << outDir << endl;

  for (i=0; i<nBlocks; i++) {
    char num[64];
    sprintf(num, "%d", i);
    string submit;
    if (bLSF)
      submit = "bsub -o " + outDir + "/HomologyByXCorr." + num + ".out " + exec_dir;
    else
      submit = exec_dir;

    string out = outDir + "/xcorr_aligns." + num + ".data";
    string cmd = submit + simple + " -o " + out;
    string chunk = " -nblocks ";
    chunk += blocks;
    chunk += " -block ";
    chunk += num;
    cmd += chunk;
    if (!bLSF)
      cmd += " > " + outDir + "/HomologyByXCorr." + num + ".out &";
    xcOutput.push_back(out);
    xcCommands.push_back(cmd);
  }


}



//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{
  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<string> oStringCmmd("-o","output directory");
  commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
  commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
  commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);
  commandArg<int> blockCmmd("-n","number of blocks", 1);
  commandArg<bool> lsfCmmd("-lsf","submit jobs to LSF", false);
  commandArg<bool> nogoCmmd("-nosubmit","do not run jobs", false);
  commandArg<bool> nowaitCmmd("-nowait","do not wait for jobs", false);
  commandArg<bool> chainCmmd("-chain_only","only chain the matches", false);
  commandArg<bool> refineCmmd("-refine_only","only refine the matches", false);
  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.99999);
  commandArg<bool> protCmmd("-proteins","align in protein space", false);
  commandArg<bool> sameCmmd("-same_only","only align sequences that have the same name.", false);
  commandArg<double> cutoffCmmd("-cutoff","signal cutoff", 1.8);
  commandArg<bool> selfCmmd("-self","ignore self-matches.", false);


  commandLineParser P(argc,argv);
  P.SetDescription("Wrapper around high-sensitivity alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(lIntCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(qChunkCmmd);
  P.registerArg(blockCmmd);
  P.registerArg(lsfCmmd);
  P.registerArg(nogoCmmd);
  P.registerArg(nowaitCmmd);
  P.registerArg(chainCmmd);
  P.registerArg(refineCmmd);
  P.registerArg(probCmmd);
  P.registerArg(protCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(selfCmmd);


  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  int nBlocks = P.GetIntValueFor(blockCmmd);
  bool bLSF = P.GetBoolValueFor(lsfCmmd);
  bool bNogo = P.GetBoolValueFor(nogoCmmd);
  bool bNowait = P.GetBoolValueFor(nowaitCmmd);
  bool bChainOnly = P.GetBoolValueFor(chainCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);
  bool bProteins = P.GetBoolValueFor(protCmmd);
  bool bRefOnly = P.GetBoolValueFor(refineCmmd);
  bool bSameOnly = P.GetBoolValueFor(sameCmmd);
  bool bSelf = P.GetBoolValueFor(selfCmmd);
  double sigCutoff = P.GetDoubleValueFor(cutoffCmmd);

  int i, j;


  string probe = output + "/satsuma.log";
  FILE * pProbe = fopen(probe.c_str(), "r");
  if (pProbe != NULL) {
    cout << "ERROR: Please remove all files from the output directory!!" << endl;
    return -1;
  }
  pProbe = fopen(probe.c_str(), "w");
  if (pProbe == NULL) {
    string mk_cmd = "mkdir " + output;
    cout << "Creating firectory " << output << endl;
    system(mk_cmd.c_str());
    pProbe = fopen(probe.c_str(), "w");
    fprintf(pProbe, "Starting\n");
    fclose(pProbe);
  } else {
    fclose(pProbe);
  }



  //cout << "outdir " << output << endl;

  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;

  if (bLSF && strstr(pExec, "/") == NULL) {
    cout << "When submitting to LSF, Satsuma must be run with the full path (e.g. '/usr/home/dummy/Satsuma...'" << endl;
    return -1;
  }

  if (output == "") {
    cout << "You MUST specify a valid output directory (option '-o')!" << endl;
    return -1;
  }


  char exec_dir[1024];
  strcpy(exec_dir, pExec);
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      break;
    }
  }

  
  svec<string> xcCommands;
  svec<string> xcOutput;

  string simple = "HomologyByXCorr -q " + sQuery;
  simple += " -t " + sTarget;
  char tmp[1024];
  sprintf(tmp, " -l %d -q_chunk %d -t_chunk %d -min_prob %f -cutoff %f ", minLen, queryChunk, targetChunk, minProb, sigCutoff);

  simple += tmp;

  if (bProteins) {
    simple += " -proteins ";
  }
  if (bSameOnly) {
    simple += " -same_only ";
  }

  cout << "Will process XCorr alignments in " << nBlocks << " processes." << endl;
  //cout << "outdir " << output << endl;
  BuildList(xcCommands, xcOutput, simple, output, nBlocks, bLSF, exec_dir);


  usleep(1000);

  if (!bNogo && !bChainOnly && !bRefOnly) {
    if (!bLSF) {
      if (xcCommands.isize() > 48) {
	cout << "ERROR: As a precaution, you cannot run more than 48 jobs locally!" << endl;
	exit (-1);
      }
    }
    for (i=0; i<xcCommands.isize(); i++) {
      cout << "Running:" << endl;
      cout << "     " << xcCommands[i] << endl << endl;
      system(xcCommands[i].c_str());
      if (i > 0 && (i % 25 == 0)) {
	cout << "Sleeping before submitting more..." << endl;
	//for (int x=0; x<10; x++)
	sleep(30);
      }
    }
  } else {
    cout << "NOT submitting jobs!!" << endl;
  }

  if (bNowait) {
    cout << "NOT waiting for jobs, terminating." << endl;
    return 0;
  }

  if (strstr(exec_dir, "/") == NULL)
    strcpy(exec_dir, ".");

  if (!bChainOnly && !bRefOnly) {
    cout << "Waiting for jobs to complete." << endl;
  
    MultiMatches multi;
    
    for (i=0; i<xcOutput.isize(); i++) {
      FILE * pTest = NULL;
      while (pTest == NULL) {
	pTest = fopen(xcOutput[i].c_str(), "rb");
	if (pTest != NULL) {
	  sleep(1);
	  multi.MergeRead(xcOutput[i]);
	  cout << "Reading output file " << xcOutput[i] << endl;
	} else {
	  sleep(1);
	  //cout << "." << endl;
	}
      }
      
      fclose(pTest);
      
    }
    cout << "Collected all data!!" << endl;
  
    cout << "Sorting output." << endl;
    multi.Sort();
    
    cout << "Saving alignments to " << output << "/xcorr_aligns.all.out" << endl;
    multi.Write(output + "/xcorr_aligns.all.out");
    multi.Clear();

    //cout << "Running merge..." << endl;
    string mergeCmd = exec_dir;
    
    mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.all.out";
    mergeCmd += " -q " + sQuery + " -t " + sTarget;
    mergeCmd += " -o " + output + "/satsuma_summary.out";
    mergeCmd += " > " + output + "/MergeXCorrMatches.out";
    
    if (bSelf) 
      mergeCmd += " -self";

    cout << "Running " << mergeCmd << endl;
    system(mergeCmd.c_str());
  } else {
    cout << "Skipping waiting, merging and merging." << endl;
  }



  if (!bRefOnly) {
    string mergeCmd = exec_dir;
    mergeCmd += "/ChainMatches -i " + output + "/xcorr_aligns.all.out";
    mergeCmd += " -o " + output + "/xcorr_aligns.chained.out";
    cout << "Running " << mergeCmd << endl;
    system(mergeCmd.c_str());



    mergeCmd = exec_dir;
    
    mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.chained.out";
    mergeCmd += " -q " + sQuery + " -t " + sTarget;
    mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
    mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";

    if (bSelf) 
      mergeCmd += " -self";
   
    cout << "Running " << mergeCmd << endl;
    system(mergeCmd.c_str());
  }
 
  cout << "Running refined matches (fill in the blanks)." << endl;
  string refCmd = exec_dir;
  refCmd += "/HomologyByXCorr -cutoff 1.2  -q " + sQuery;
  refCmd += " -t " + sTarget;
  refCmd += " -o " + output + "/xcorr_aligns.refined.out";
  refCmd += " -guide " + output + "/xcorr_aligns.chained.out";
  refCmd += " > " + output + "/HomologyByXCorr.refined.out";

  cout << "Running " << refCmd << endl;
  system(refCmd.c_str());
  
  
  refCmd = exec_dir;
  
  refCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.refined.out";
  refCmd += " -q " + sQuery + " -t " + sTarget;
  refCmd += " -o " + output + "/satsuma_summary.refined.out";
  refCmd += " > " + output + "/MergeXCorrMatches.refined.out";
  
  if (bSelf) 
    refCmd += " -self";


  cout << "Running " << refCmd << endl;
  system(refCmd.c_str());
 


  cout << "all done!" << endl;
  

  return 0;
}

