#ifndef FORCE_DEBUG
#define NDEBUG
#endif


//#define FAKE

#include <string.h>
#include <unistd.h>
#include "base/CommandLineParser.h"
#include "analysis/SequenceMatch.h"
#include "analysis/GridSearch.h"
#include "analysis/SeqChunk.h"
#include "analysis/MatchDynProg.h"
#include "util/SysTime.h"
#include <math.h>
#include <pthread.h>
#include <netinet/in.h> 
#include <sys/socket.h>
#include "analysis/WorkQueue.h"


class SyntenyInterpolator
{
public:
  SyntenyInterpolator() {
    m_min = 40;
    m_num = 2;
    m_space = 100000;
  }

  int Interpolate(MultiMatches & chained); 

private:
  double Slope(double & sx, double & sy, const MultiMatches & chained, int first, int last) {
    int n = chained.GetMatchCount();
    int i;
    double eps = 0;
    double div = 0;
    // Idiotic...
    //cout << "Matches: " << last - first << endl;
    for (i=first+1; i<=last; i++) {
      const SingleMatch & m = chained.GetMatch(i);
      //cout << m.GetQueryID() << " " << m.GetStartQuery();
  /*    if (m.IsRC())
	cout << " -" << endl;
      else
	cout << " +" << endl;*/
      const SingleMatch & l = chained.GetMatch(i-1);

      double x = m.GetStartTarget() - l.GetStartTarget();
      double y = m.GetStartQuery() - l.GetStartQuery();
      double s = y/x;
      s = atan(s);
      eps += s;
      div += 1.;
    }
    
    eps /= div;
    sy = sin(eps);
    sx = cos(eps);
    //cout << "Slope: eps=" << eps << " x=" << sx << " y=" << sy << endl; 
    return eps;
  }

  int m_min;
  int m_num;
  int m_space;
};





int SyntenyInterpolator::Interpolate(MultiMatches & chained) 
{
  int n = chained.GetMatchCount();
  int i, j;

  //MultiMatches added;

  int last = 0;
  int count = 0;
  for (i=1; i<n; i++) {
    const SingleMatch & m = chained.GetMatch(i);
    const SingleMatch & l = chained.GetMatch(i-1);

    int x1 = m.GetStartTarget();
    int y1 = m.GetStartQuery();
    int x2 = l.GetStartTarget();
    int y2 = l.GetStartQuery();
    
    
	
    if (m.GetTargetID() == l.GetTargetID() &&
	m.GetQueryID() == l.GetQueryID() &&
	m.IsRC() == l.IsRC() &&
	x1 - x2 < m_space/2 &&
	y1 - y2 < m_space/2 && y1 > y2) { 
      if (x1 - x2 > l.GetLength() &&
	  y1 - y2 > l.GetLength()) {
	count++;
      }

    } else {
      if (count >= m_min) {
	double sx = 0.; 
	double sy = 0;
	Slope(sx, sy, chained, last, i-1);
	SingleMatch first = chained.GetMatch(last);
	SingleMatch end = chained.GetMatch(i-1);
	int x = first.GetStartTarget();
	int y = first.GetStartQuery();
	for (j=1; j<=m_num; j++) {
	  int xx = x - (int)((double)m_space * (double)j * sx); 
	  int yy = y - (int)((double)m_space * (double)j * sy); 
	  SingleMatch dd;
	  dd.SetProbability(0.999);
	  dd.SetIdentity(0.999);

	  dd.SetQueryTargetID(first.GetQueryID(), first.GetTargetID(), 
			      chained.GetQuerySize(first.GetQueryID()));
	  dd.SetPos(yy, xx, 50+i-last, first.IsRC());
	  chained.AddMatch(dd);
	  //cout << "Adding extra (left): " << xx << " / " << yy << endl;
	  

	}
	x = end.GetStartTarget();
	y = end.GetStartQuery();
	for (j=1; j<=m_num; j++) {
	  int xx = x + (int)((double)m_space * (double)j * sx); 
	  int yy = y + (int)((double)m_space * (double)j * sy);
 
	  SingleMatch dd;
	  dd.SetProbability(0.999);
	  dd.SetIdentity(0.999);

	  dd.SetQueryTargetID(end.GetQueryID(), end.GetTargetID(), 
			      chained.GetQuerySize(end.GetQueryID()));
	  dd.SetPos(yy, xx, 50+i-last, end.IsRC());
	  chained.AddMatch(dd);
	  cout << "Adding extra (right): " << xx << " / " << yy << endl;

	}
	
      }
      last = i;
      count = 0;
    }
    
  }
  return 0;
}



//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{
  //ALG: parse arguments (TODO: review arguments)
  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<string> oStringCmmd("-o","output directory");
  commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
  commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
  commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);
  commandArg<int> slavesCmmd("-slaves","number of processing slaves", 1);
  commandArg<bool> lsfCmmd("-lsf","submit jobs to LSF", false);
  commandArg<int> perCmmd("-m","number of jobs per block", 4);
  commandArg<bool> refineNotCmmd("-do_refine","refinment steps", false);
  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.99999);
  commandArg<double> cutoffCmmd("-cutoff","signal cutoff", 1.8);
  commandArg<string> seedCmmd("-seed","loads seeds and runs from there (xcorr*data)", "");
  commandArg<int> blockPixelCmmd("-pixel","number of blocks per pixel", 24);
  commandArg<bool> filterCmmd("-nofilter","do not pre-filter seeds (slower runtime)", false);
  commandArg<bool> dupCmmd("-dups","allow for duplications in the query sequence", false);


  commandLineParser P(argc,argv);
  P.SetDescription("Wrapper around high-sensitivity alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(lIntCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(qChunkCmmd);
  P.registerArg(lsfCmmd);
  P.registerArg(refineNotCmmd);
  P.registerArg(probCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(perCmmd);
  P.registerArg(slavesCmmd);
  P.registerArg(seedCmmd);
  P.registerArg(blockPixelCmmd);
  P.registerArg(filterCmmd);
  P.registerArg(dupCmmd);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string seedFile = P.GetStringValueFor(seedCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  int perBlock = P.GetIntValueFor(perCmmd);
  int slave_count = P.GetIntValueFor(slavesCmmd);
  bool bLSF = P.GetBoolValueFor(lsfCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);
  bool bNoRef = P.GetBoolValueFor(refineNotCmmd);
  double sigCutoff = P.GetDoubleValueFor(cutoffCmmd);
  int blocksPerPixel = P.GetIntValueFor(blockPixelCmmd);
  bool bFilter = P.GetBoolValueFor(filterCmmd);
  bool bDup = P.GetBoolValueFor(dupCmmd);
  MultiMatches matches;
  
  int i, j;


  cout << "SATSUMA: Welcome to SatsumaSynteny! Current date and time: " << GetTimeStatic() << endl;

  string satsuma2_path(std::getenv("SATSUMA2_PATH"));
  string current_path(std::getenv("PWD"));
  if (satsuma2_path==""){
    cout << "ERROR: SATSUMA2_PATH variable not set, please set it to the binary path." <<endl;
    return -1;
  }
  cout<< "Path for Satsuma2: '"<<satsuma2_path<<"'"<<endl;
  //TODO: test for the binaries to be there!

  //ALG: create output dir
  string probe = output + "/satsuma.log";
  FILE * pProbe = fopen(probe.c_str(), "r");
  if (pProbe != NULL) {
    cout << "ERROR: Please remove all files from the output directory!!" << endl;
    return -1;
  }
  pProbe = fopen(probe.c_str(), "w");
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


  //==================================================================


  //ALG: setting up GRID
  vecDNAVector target, query;

  query.Read(sQuery);
  target.Read(sTarget);

  cout << "Setting up grid." << endl;
  GridSearch grid(targetChunk, blocksPerPixel);
  cout << "Preparing..." << endl;
  grid.SetUp(target, query);
  
  cout << "Initializing multimatches..." << endl;
  svec<string> queryNames, targetNames;

  for (i=0; i<target.isize(); i++)
    targetNames.push_back(target.NameClean(i));
  for (i=0; i<query.isize(); i++)
    queryNames.push_back(query.NameClean(i));
  matches.SetNames(queryNames, targetNames);


  for (i=0; i<target.isize(); i++) {
    matches.SetTargetSize(i, target[i].isize());
  }
  for (i=0; i<query.isize(); i++)
    matches.SetQuerySize(i, query[i].isize());
  cout << "Done." << endl;
  target.clear();
  query.clear();
  
  
  bool bNoChain = false;


  WorkQueue wq(minLen, sQuery, queryChunk, sTarget, targetChunk, minProb, sigCutoff, slave_count);
  //==================================================================
  //ALG: create filtered seeds
  if (!bFilter && seedFile == "") {
    for (long long i=11; i<32; i+=2){
      seedFile = output + "/kmatch_results.k"+ to_string(i);
      string cmd;
      cmd = "echo \"cd " + current_path + ";"+ satsuma2_path + "/KMatch " + sQuery + " " + sTarget;
      cmd += " " + to_string(i) + " " + seedFile + " " + to_string(i) + " " + to_string(i-1);
      cmd += "; touch " + seedFile + ".finished\"|qsub -l ncpus=2,mem=100G";
      cout << "Running seed pre-filter " << cmd << endl;
      system(cmd.c_str());
    }
    //TODO:wait for the kmatches to finish
    cout << "Waiting for seed pre-filters..." << endl;
    seedFile = output + "/kmatch_results.k";
    int finished=0;
    while (finished<11){
      sleep(3);
      for (long long i=11; i<32; i+=2){
        FILE * pProbe = fopen((seedFile+to_string(i)+".finished").c_str(), "r");
        if (pProbe != NULL) {
          fclose(pProbe);
          cout<<"loading results for k="<<i<<endl;
          wq.results_from_file((seedFile+to_string(i)).c_str());
          string cmd="rm "+seedFile+to_string(i)+".finished";
          system(cmd.c_str());
          finished++;
        }
      }
    }
    cout << "Seed pre-filters finished" << endl;
  }
  
  //ALG: processing the seeds. Why is this different to just process matches? (besides the filter)
  
  unsigned long int new_matches_count;
  wq.collect_new_matches(matches); //matches.Read(seedFile);
  matches.Sort();
  matches.Collapse();
  matches.LengthFilter(24);//XXX:hardcoded lenght filter-->terrible!!!

  //=============================================================

  cout << "SATSUMA: Launching slaves, date and time: " << GetTimeStatic() << endl;
  wq.setup_queue();//XXX totally wrong to name the slaves count like that!
  
  //MultiMatches workMatches;
  cout << "SATSUMA: Entering main search loop, date and time: " << GetTimeStatic() << endl;
  int exitCounter = 0;
  int superExitCounter = 0;
  int superExitCounterThresh = 20;
  int superExitCounterNot = 0;
  int targets_queue_size=slave_count*perBlock*3;
  int main_iteration=0;
  //ALG: main loop
  while (true) {
    //TODO: ALG: collect matches (get count of tasks and results)
    main_iteration++;
    cout<<"MAIN: starting iteration "<<main_iteration<<endl;
    t_collect_status collect_status=wq.collect_new_matches(matches);
    cout<<"MAIN: "<<collect_status.matches<<" new matches collected"<<endl;
    //ALG: if new matches
    if (collect_status.matches > 0) {
      //ALG: chain matches
      cout << "MAIN: Running the chaining step..." << endl;
      MultiMatches chained;
      matches.Sort();
      if (bDup)
        RunMatchDynProgMult(chained, matches);
      else
        RunMatchDynProg(chained, matches);
      cout << "MAIN: Total matches: " << chained.GetMatchCount() << endl;
      matches = chained;
      SyntenyInterpolator inter;
      inter.Interpolate(chained); 
      cout << "MAIN: Total matches (interpolated): " << chained.GetMatchCount() << endl;
      cout << "MAIN: Updating grid's Target Weights" << endl;
      grid.UpdateTargetWeights(chained);
    }
    //TODO: ALG: collect targets to fill in the queue
    int targets_to_collect=targets_queue_size - wq.pending_pair_count();
    if (!collect_status.matches && !targets_to_collect) {
      cout<< "MAIN: nothing changed, skipping cycle and waiting 3 seconds"<<endl;
      sleep(3);
      continue;
    }
    cout << "MAIN: Collecting " << targets_to_collect << " new targets from the grid " << endl;
    svec<GridTarget> newTargets;
    int realTargets = grid.CollectTargets(newTargets, targets_to_collect);
    cout << "MAIN: Targets retrieved: " << newTargets.isize() << endl;
    for (i=0; i<newTargets.isize(); i++) {
        wq.add_pair(newTargets[i].TargetFirst(),
            newTargets[i].TargetLast(),
            newTargets[i].QueryFirst(),
            newTargets[i].QueryLast(),
            newTargets[i].IsFast());
    }

    //TODO: ALG: check CONVERGENCE condition
    cout<< "MAIN: STATUS: slaves_processed, new_matches, slave_to_match_relation, targets_needed, targets_collected"<<endl;
    cout<< "MAIN: STATUS: "<< collect_status.slaves << ", " \
                           << collect_status.matches << ", " \
                           << (collect_status.slaves ? ((double) collect_status.matches)/collect_status.slaves : 0 ) << ", " \
                           << targets_to_collect << ", " \
                           << newTargets.isize() << endl;
    if (newTargets.isize()+wq.pending_pair_count()==0){//XXX: extreme convergence, only useful to test functions
      break;
    }



  }

  wq.close_queue(); //XXX: TODO: wait till all targets have been processed!

  //ALG: write final output
  cout << "Writing final output!" << endl;
  string finalOut = output + "/xcorr_aligns.almost.out"; 
  matches.Write(finalOut);


  //ALG: run ChainMatches
  string mergeCmd = satsuma2_path;

  mergeCmd += "/ChainMatches -i " + output + "/xcorr_aligns.almost.out";
  mergeCmd += " -o " + output + "/xcorr_aligns.final.out";

  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());


  //ALG: run /MergeXCorrMatches
  mergeCmd = satsuma2_path;

  mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.final.out";
  mergeCmd += " -q " + sQuery + " -t " + sTarget;
  mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
  mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";

  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());

  //ALG: if bNoRef, refine matches
  if (bNoRef) {
    cout << "Running refined matches (fill in the blanks)." << endl;
    string refCmd = satsuma2_path;
    refCmd += "/HomologyByXCorr -cutoff 1.2  -q " + sQuery;
    refCmd += " -t " + sTarget;
    refCmd += " -o " + output + "/xcorr_aligns.refined.out";
    refCmd += " -guide " + output + "/xcorr_aligns.final.out";
    refCmd += " > " + output + "/HomologyByXCorr.refined.out";

    cout << "Running " << refCmd << endl;
    system(refCmd.c_str());


    refCmd = satsuma2_path;

    refCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.refined.out";
    refCmd += " -q " + sQuery + " -t " + sTarget;
    refCmd += " -o " + output + "/satsuma_summary.refined.out";
    refCmd += " > " + output + "/MergeXCorrMatches.refined.out";

    cout << "Running " << refCmd << endl;
    system(refCmd.c_str());
  }


  cout << "all done!" << endl;



  cout << "SATSUMA: all done, date and time: " << GetTimeStatic() << endl;
  cout << "all done!" << endl;


  return 0;
}

