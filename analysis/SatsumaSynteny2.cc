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
  commandArg<int> qSeedChunkCmmd("-q_chunk_seed","query chunk size (seed)", 8192);
  commandArg<int> tSeedChunkCmmd("-t_chunk_seed","target chunk size (seed)", 8192);
  commandArg<int> blockCmmd("-n","number of processes", 1);
  commandArg<int> blockIniCmmd("-ni","number of initial search blocks", -1);
  commandArg<bool> lsfCmmd("-lsf","submit jobs to LSF", false);
  commandArg<bool> lsfCmmdIni("-lsf_ini","submit jobs to LSF", false);
  commandArg<int> perCmmd("-m","number of jobs per block", 32);
  commandArg<int> slavepathCmmd("-slave_path","full path to the slave binary");
  commandArg<int> slavesCmmd("-slaves","number of processing slaves", 1);
  commandArg<int> threadsCmmd("-threads","number of threads in each slave", 4);
  commandArg<bool> nogoCmmd("-nosubmit","do not run jobs", false);
  commandArg<bool> nowaitCmmd("-nowait","do not wait for jobs", false);
  commandArg<bool> chainCmmd("-chain_only","only chain the matches", false);
  commandArg<bool> refineCmmd("-refine_only","only refine the matches", false);
  commandArg<bool> refineNotCmmd("-do_refine","refinment steps", false);
  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.99999);
  commandArg<bool> protCmmd("-proteins","align in protein space", false);
  commandArg<double> cutoffCmmd("-cutoff","signal cutoff", 1.8);
  commandArg<double> cutoffSeedCmmd("-cutoff_seed","signal cutoff (seed)", 2.0);
  commandArg<string> resumeCmmd("-resume","resumes w/ the output of a previous run (xcorr*data)", "");
  commandArg<string> seedCmmd("-seed","loads seeds and runs from there (xcorr*data)", "");
  commandArg<int> blockPixelCmmd("-pixel","number of blocks per pixel", 24);
  commandArg<bool> filterCmmd("-nofilter","do not pre-filter seeds (slower runtime)", false);
  commandArg<string> distanceCmmd("-seeddist","distance between pre-filter seeds (increase for close genomes)", "1");
  commandArg<string> seedWidthCmmd("-filterwidth","width of the seed filter", "2");
  commandArg<bool> dupCmmd("-dups","allow for duplications in the query sequence", false);


  commandLineParser P(argc,argv);
  P.SetDescription("Wrapper around high-sensitivity alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(lIntCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(qChunkCmmd);
  P.registerArg(tSeedChunkCmmd);
  P.registerArg(qSeedChunkCmmd);
  P.registerArg(blockCmmd);
  P.registerArg(blockIniCmmd);
  P.registerArg(lsfCmmd);
  P.registerArg(lsfCmmdIni);
  P.registerArg(nogoCmmd);
  P.registerArg(nowaitCmmd);
  P.registerArg(chainCmmd);
  P.registerArg(refineCmmd);
  P.registerArg(refineNotCmmd);
  P.registerArg(probCmmd);
  P.registerArg(protCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(cutoffSeedCmmd);
  P.registerArg(perCmmd);
  P.registerArg(slavepathCmmd);
  P.registerArg(slavesCmmd);
  P.registerArg(threadsCmmd);
  P.registerArg(resumeCmmd);
  P.registerArg(seedCmmd);
  P.registerArg(blockPixelCmmd);
  P.registerArg(filterCmmd);
  P.registerArg(distanceCmmd);
  P.registerArg(dupCmmd);
  P.registerArg(seedWidthCmmd);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string prevRun = P.GetStringValueFor(resumeCmmd);
  string seedFile = P.GetStringValueFor(seedCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  int targetChunkSeed = P.GetIntValueFor(tSeedChunkCmmd);
  int queryChunkSeed = P.GetIntValueFor(qSeedChunkCmmd);
  int nBlocks = P.GetIntValueFor(blockCmmd);
  int nBlocksIni = P.GetIntValueFor(blockIniCmmd);
  int perBlock = P.GetIntValueFor(perCmmd);
  string slave_path = P.GetStringValueFor(slavepathCmmd);
  int slave_count = P.GetIntValueFor(slavesCmmd);
  int slave_threads = P.GetIntValueFor(threadsCmmd);
  bool bLSF = P.GetBoolValueFor(lsfCmmd);
  bool bLSFIni = P.GetBoolValueFor(lsfCmmdIni);
  bool bNogo = P.GetBoolValueFor(nogoCmmd);
  bool bNowait = P.GetBoolValueFor(nowaitCmmd);
  bool bChainOnly = P.GetBoolValueFor(chainCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);
  bool bProteins = P.GetBoolValueFor(protCmmd);
  bool bRefOnly = P.GetBoolValueFor(refineCmmd);
  bool bNoRef = P.GetBoolValueFor(refineNotCmmd);
  double sigCutoff = P.GetDoubleValueFor(cutoffCmmd);
  double sigCutoffSeed = P.GetDoubleValueFor(cutoffSeedCmmd);
  int blocksPerPixel = P.GetIntValueFor(blockPixelCmmd);
  bool bFilter = P.GetBoolValueFor(filterCmmd);
  string theDistance = P.GetStringValueFor(distanceCmmd);
  string seedwidth = P.GetStringValueFor(seedWidthCmmd);
  bool bDup = P.GetBoolValueFor(dupCmmd);
  MultiMatches matches;
  
  if (nBlocksIni == -1)
    nBlocksIni = nBlocks;
  int i, j;

  //Terribly bad fix!




  cout << "SATSUMA: Welcome to SatsumaSynteny! Current date and time: " << GetTimeStatic() << endl;
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


  char exec_dir[8192];
  strcpy(exec_dir, pExec);
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      break;
    }
  }


  //==================================================================
  //ALG: create filtered seeds
  if (!bFilter && seedFile == "") {
    seedFile = output + "/xcorr_aligns.seeds.out";

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
  }
  

  //==================================================================

  //ALG: spawn slaves

  WorkQueue wq(minLen, sQuery, queryChunk, sTarget, targetChunk, minProb, sigCutoff, slave_path, slave_count, slave_threads);
  wq.setup_queue();//XXX totally wrong to name the slaves count like that!

  


  //ALG: setting up GRID
  vecDNAVector target, query;

  query.Read(sQuery);
  target.Read(sTarget);

  cout << "Setting up grid." << endl;
  GridSearch grid(targetChunk, blocksPerPixel);
  cout << "Preparing..." << endl;
  grid.SetUp(target, query);

  cout << "Done." << endl;
  target.clear();
  query.clear();

  //ALG: adding the seeds to the workqueue
  
  int collectCounter = 0;
  bool bNoChain = false;
  cout << "Loading seeds from " << seedFile << endl;
  matches.Read(seedFile);
  matches.Sort();
  matches.Collapse();

  MultiMatches tmp;
  if (bDup)
    RunMatchDynProgMult(tmp, matches);
  else
    RunMatchDynProg(tmp, matches);

  bNoChain = true;
  cout << "Done!" << endl;
  matches = tmp;



  //=============================================================


  
  //MultiMatches workMatches;
  cout << "SATSUMA: Entering main search loop, date and time: " << GetTimeStatic() << endl;
  int last = 0;
  int exitCounter = 0;
  int superExitCounter = 0;
  int superExitCounterThresh = 20;
  int superExitCounterNot = 0;
  int counter = 0;
  int plus_counter = 0;

  int last_pair_in_current_batch=0;
  unsigned long int pending_matches=1;//XXX change this to a meaningfull value
  //ALG: main loop
  while (true) {
    //ALG: wait for slaves to finish, and update the matches count
    while (pending_matches == 0  || wq.pending_pair_count() >= slave_count * slave_threads * 4 ){
        pending_matches+=wq.collect_new_matches(matches);
        cout<<"Main loop waiting (pending_matches="<<pending_matches<<" pending_pairs="<<wq.pending_pair_count()<<") ..."<<endl;
        sleep(1); //XXX: needs to have something to do if no new matches because none was good
    }
    pending_matches=0;
    

    bool bSetUp = false;

    //ALG: if new matches
    if (matches.GetMatchCount() - last > 0) {
      bSetUp = true;
      
      //ALG: chain matches
      cout << "Running the chaining step..." << endl;
      
      MultiMatches chained;
      matches.Sort();

      if (!bNoChain) {
        if (bDup)
          RunMatchDynProgMult(chained, matches);
        else
          RunMatchDynProg(chained, matches);
      } else {
        chained = matches;
        bNoChain = false;
      }

      counter++;

      int startMatches = chained.GetMatchCount();
      cout << "Total matches: " << chained.GetMatchCount() << endl;
      matches = chained;

      SyntenyInterpolator inter;
      inter.Interpolate(chained); 
      cout << "Total matches (interpolated): " << chained.GetMatchCount() << endl;

      grid.ClearTargetWeights();

      //ALG: for each chained match
      int lastY1 = -1;
      int lastY2 = -1;
      cout<<"Considering "<<chained.GetMatchCount()<<" matches..."<<endl;
      for (i=0; i<chained.GetMatchCount(); i++) {
        const SingleMatch & m = chained.GetMatch(i);

        int x1 = m.GetStartTarget();
        int y1 = m.GetStartQuery();
        int x2 = m.GetStartTarget() + m.GetLength();
        int y2 = m.GetStartQuery() + m.GetLength();

        if (m.IsRC()) {
          int s = chained.GetQuerySize(m.GetQueryID());
          int tmp = y1;
          y1 = s - y2;
          y2 = s - tmp;
        }

        if (x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0)
          continue;
        if (x2 >= chained.GetTargetSize(m.GetTargetID()) ||
            y1 >= chained.GetQuerySize(m.GetQueryID()) ||
            y2 >= chained.GetQuerySize(m.GetQueryID())) {
          continue;
        }
        //ALG: consider 
        //cout << "Considering." << endl;
        grid.ConsiderTargets(m.GetTargetID(), x1, x2, m.GetQueryID(), y1, y2, m.GetIdentity());
        lastY1 = y1;
        lastY2 = y2;
      }
      cout<<"Matches considered."<<endl;


      plus_counter++;


    }
    //ALG: collect targets
    last = matches.GetMatchCount();
    svec<GridTarget> newTargets;
    collectCounter++;
    int seedPlus = 0;
    cout << "Collecting." << endl;
    //int realTargets = grid.CollectTargets(newTargets, nBlocks * perBlock, seedPlus);
    //TODO: ckech the size of the queue here, implement the proposed changes
    nBlocks=4*slave_count*slave_threads;
    int realTargets = grid.CollectTargets(newTargets, 4*slave_count*slave_threads, seedPlus);
    cout << "Targets retrieved: " << newTargets.isize() << endl;
    cout << "Ready to re-feed." << endl;



    int giveUpThresh = (nBlocks * perBlock)/4+1;

    //ALG: check end condition
    if (newTargets.isize() >= nBlocks * perBlock)
      superExitCounterNot++;

    if (newTargets.isize() < nBlocks * perBlock) {
      if (superExitCounterNot > 4 * superExitCounterThresh)
        superExitCounter++;
    }


    if (realTargets < giveUpThresh 
        || superExitCounter >= superExitCounterThresh
        || newTargets.isize() == 0) {
      cout << "Counting down: " << exitCounter << endl;
      exitCounter++;

      if (exitCounter > 10 || newTargets.isize() == 0) {
        //ALG: final collect
        cout << "Final collect!" << endl;
        string die = output + "/slaves.directive";
        FILE * pDie = fopen(die.c_str(), "w");
        fprintf(pDie, "exit\n");
        fclose(pDie);

        //ALG: wait for slaves to die
        //wq.close_queue();
        cout << "All done, exiting!!" << endl;
        //ALG: exit main loop
        break;
      }
    } else {

      cout << "Resetting exit counter." << endl;
      exitCounter = 0;
    }


    // Check for alive processesss.
    //ALG: distribute the work in blocks 

    for (i=0; i<newTargets.isize(); i++) {
      int to = i + perBlock;
      if (to >= newTargets.isize())
        to = newTargets.isize();
      for (j=i; j<to; j++) {
        wq.add_pair(newTargets[j].TargetFirst(),
            newTargets[j].TargetLast(),
            newTargets[j].QueryFirst(),
            newTargets[j].QueryLast(),
            newTargets[j].IsFast());

      }
      i = j-1;

    }
  }



  wq.close_queue();
  //ALG: write final output
  cout << "Writing final output!" << endl;
  string finalOut = output + "/xcorr_aligns.almost.out"; 
  matches.Write(finalOut);


  //ALG: run ChainMatches
  string mergeCmd = exec_dir;

  mergeCmd += "/ChainMatches -i " + output + "/xcorr_aligns.almost.out";
  mergeCmd += " -o " + output + "/xcorr_aligns.final.out";

  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());


  //ALG: run /MergeXCorrMatches
  mergeCmd = exec_dir;

  mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.final.out";
  mergeCmd += " -q " + sQuery + " -t " + sTarget;
  mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
  mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";

  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());

  //ALG: if bNoRef, refine matches
  if (bNoRef) {
    cout << "Running refined matches (fill in the blanks)." << endl;
    string refCmd = exec_dir;
    refCmd += "/HomologyByXCorr -cutoff 1.2  -q " + sQuery;
    refCmd += " -t " + sTarget;
    refCmd += " -o " + output + "/xcorr_aligns.refined.out";
    refCmd += " -guide " + output + "/xcorr_aligns.final.out";
    refCmd += " > " + output + "/HomologyByXCorr.refined.out";

    cout << "Running " << refCmd << endl;
    system(refCmd.c_str());


    refCmd = exec_dir;

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

