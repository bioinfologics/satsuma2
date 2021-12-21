#ifndef FORCE_DEBUG
#define NDEBUG
#endif


//#define FAKE

#include <string.h>
#include <unistd.h>
#include <memory>
#include "../base/CommandLineParser.h"
#include "SequenceMatch.h"
#include "GridSearch.h"
#include "MatchDynProg.h"
#include "../util/SysTime.h"
#include "WorkQueue.h"
#define POSITION_CHR_CNST 10000000000L

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
    for (i=first+1; i<=last; i++) {
      const SingleMatch & m = chained.GetMatch(i);
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
	  //cout << "Adding extra (right): " << xx << " / " << yy << endl;

	}
	
      }
      last = i;
      count = 0;
    }
    
  }
  return 0;
}

typedef struct match_segments_s {
  uint64_t start;
  uint64_t end;
  bool operator <(const match_segments_s & rhs) const{
        return (start < rhs.start);
          }
} match_segments;

class TargetCoverageTracker{
  public:
    TargetCoverageTracker(){
      match_blocks= std::shared_ptr<std::vector<match_segments>>(new vector<match_segments>());
      cout<<"TargetCoverageTracker initialized a match_blocks of size"<<match_blocks->size()<<endl;
    }
    std::pair<uint64_t, uint64_t> update_and_report(const MultiMatches & m){
      //Create a std::shared_ptr<std::vector<match_segments>> from
      std::shared_ptr<std::vector<match_segments>> new_match_blocks = std::shared_ptr<std::vector<match_segments>>(new vector<match_segments>());
      new_match_blocks->resize(m.GetMatchCount());
      for (unsigned long int i=0;i<m.GetMatchCount();i++){
        (*new_match_blocks)[i].start=m.GetMatch(i).GetTargetID()*POSITION_CHR_CNST+m.GetMatch(i).GetStartTarget();
        (*new_match_blocks)[i].end=(*new_match_blocks)[i].start+m.GetMatch(i).GetLength();
        //if (i>0 && (*new_match_blocks)[i].start < (*new_match_blocks)[i-1].start) cout<<"ERROR: unsorted Matches into coverageTracker"<<endl;
      };
      std::sort(new_match_blocks->begin(),new_match_blocks->end());
      //Double list comparisson
      unsigned long int i=0,j=0,lastpoint=0,nextpoint=0, last_status=0;
      unsigned long int coverage[4]={0,0,0,0};
      while (i<match_blocks->size() || j<new_match_blocks->size()){

        //find next point after last_point;
        
        if (i<match_blocks->size()) nextpoint=(*match_blocks)[i].end;
        if (j<new_match_blocks->size() && (*new_match_blocks)[j].end > nextpoint) nextpoint=(*new_match_blocks)[j].end;

        if (i<match_blocks->size()){
          if ((*match_blocks)[i].start < nextpoint && (*match_blocks)[i].start > lastpoint) nextpoint=(*match_blocks)[i].start;
          else if ((*match_blocks)[i].end < nextpoint && (*match_blocks)[i].end > lastpoint) nextpoint=(*match_blocks)[i].end;
        }
        if (j<new_match_blocks->size()){
          if ((*new_match_blocks)[j].start < nextpoint && (*new_match_blocks)[j].start > lastpoint) nextpoint=(*new_match_blocks)[j].start;
          else if ((*new_match_blocks)[j].end < nextpoint && (*new_match_blocks)[j].end > lastpoint) nextpoint=(*new_match_blocks)[j].end;
        }

        //status? (who's covering this?)
        unsigned long int status=0;
        if (i<match_blocks->size() && lastpoint>=(*match_blocks)[i].start && nextpoint<=(*match_blocks)[i].end) status +=1;
        if (j<new_match_blocks->size() && lastpoint>=(*new_match_blocks)[j].start && nextpoint<=(*new_match_blocks)[j].end) status +=2;
        coverage[status]+=nextpoint-lastpoint;

        /*if (j<10){
          if (i<match_blocks->size()) cout<<"CoverageTracker Detail: OLD Match ["<<(*match_blocks)[i].start<<"-"<<(*match_blocks)[i].end<<"]"<<endl;
          else cout<<"CoverageTracker Detail: NO OLD Match"<<endl;
          if (j<new_match_blocks->size()) cout<<"CoverageTracker Detail: NEW Match ["<<(*new_match_blocks)[j].start<<"-"<<(*new_match_blocks)[j].end<<"]"<<endl;
          else cout<<"CoverageTracker Detail: NO NEW Match"<<endl;
          cout<<"CoverageTracker Detail: lastpoint="<<lastpoint<<" nextpoint="<<nextpoint<<" status="<<status<<endl;
        }*/
        
        //advance 
        lastpoint=nextpoint;
        while (i<match_blocks->size() && lastpoint>=(*match_blocks)[i].end) i++;
        while (j<new_match_blocks->size() && lastpoint>=(*new_match_blocks)[j].end) j++;
        last_status=status;
      }
      match_blocks=new_match_blocks;
      return std::pair<uint64_t, uint64_t>(coverage[2],coverage[1]);  




    }

  private:
    std::shared_ptr<std::vector<match_segments>> match_blocks;
};

std::string getEnvVar( std::string const & key )
{
  char * val = getenv( key.c_str() );
  return val == NULL ? std::string("") : std::string(val);
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
  commandArg<int> threadsCmmd("-threads","number of working threads per processing slave", 1);
  commandArg<int> kmmemCmmd("-km_mem","memory required for kmatch (Gb)",100);  // TODO: default to mammal genome requirement
  commandArg<bool> kmsyncCmmd("-km_sync","run kmatch jobs synchronously",true);
  commandArg<int> slmemCmmd("-sl_mem","memory requirement for slaves (Gb)", 100); // TODO: default to mammal genome requirement
  commandArg<int> perCmmd("-m","number of jobs per block", 4);
  commandArg<bool> refineNotCmmd("-do_refine","refinement steps", false);
  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.99999);
  commandArg<double> cutoffCmmd("-cutoff","signal cutoff", 1.8);
  commandArg<bool> probtableCmmd("-prob_table","approximate match prob using a table lookup in slaves", false);
  commandArg<int> minmatchesCmmd("-min_matches","minimum matches per target to keep iterating", 20);
  commandArg<string> seedCmmd("-seed","loads seeds and runs from there (kmatch files prefix)", "");
  commandArg<int> max_kmatch_freqCmmd("-max_seed_kmer_freq","maximum frequency for kmatch seed kmers", 1);
  commandArg<int> min_seedCmmd("-min_seed_length","minimum length for kmatch seeds (after collapsing)", 24);
  commandArg<string> old_seedCmmd("-old_seed","loads seeds and runs from there (xcorr*data)", "");
  commandArg<int> blockPixelCmmd("-pixel","number of blocks per pixel", 24);
  commandArg<bool> filterCmmd("-nofilter","do not pre-filter seeds (slower runtime)", false);
  commandArg<bool> dupCmmd("-dups","allow for duplications in the query sequence", false);
  commandArg<bool> dumpCycleMatchesCmmd("-dump_cycle_matches","dump matches on each cycle (for debug/testing)", false);


  commandLineParser P(argc,argv);
  P.SetDescription("Wrapper around high-sensitivity alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(lIntCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(qChunkCmmd);
  P.registerArg(slmemCmmd);
  P.registerArg(refineNotCmmd);
  P.registerArg(probCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(probtableCmmd);
  P.registerArg(minmatchesCmmd);
  P.registerArg(perCmmd);
  P.registerArg(slavesCmmd);
  P.registerArg(threadsCmmd);
  P.registerArg(kmmemCmmd);
  P.registerArg(kmsyncCmmd);
  P.registerArg(seedCmmd);
  P.registerArg(min_seedCmmd);
  P.registerArg(max_kmatch_freqCmmd);
  P.registerArg(old_seedCmmd);
  P.registerArg(blockPixelCmmd);
  P.registerArg(filterCmmd);
  P.registerArg(dupCmmd);
  P.registerArg(dumpCycleMatchesCmmd);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string seedFile = P.GetStringValueFor(seedCmmd);
  int min_seed_length = P.GetIntValueFor(min_seedCmmd);
  int max_kmatch_freq = P.GetIntValueFor(max_kmatch_freqCmmd);
  string old_seedFile = P.GetStringValueFor(old_seedCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  int perBlock = P.GetIntValueFor(perCmmd);
  int slave_count = P.GetIntValueFor(slavesCmmd);
  int threads_per_slave = P.GetIntValueFor(threadsCmmd);
  int kmatch_mem = P.GetIntValueFor(kmmemCmmd);
  int slave_mem = P.GetIntValueFor(slmemCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);
  bool bNoRef = P.GetBoolValueFor(refineNotCmmd);
  double sigCutoff = P.GetDoubleValueFor(cutoffCmmd);
  bool probtable=P.GetBoolValueFor(probtableCmmd);
  int min_matches_per_target=P.GetIntValueFor(minmatchesCmmd);
  int blocksPerPixel = P.GetIntValueFor(blockPixelCmmd);
  bool bFilter = P.GetBoolValueFor(filterCmmd);
  bool bDup = P.GetBoolValueFor(dupCmmd);
  bool bDumpCycleMatches = P.GetBoolValueFor(dumpCycleMatchesCmmd);
  bool bKm_sync = P.GetBoolValueFor(kmsyncCmmd);
  MultiMatches matches;
  
  int i, j;


  cout << "SATSUMA: Welcome to SatsumaSynteny! Current date and time: " << GetTimeStatic() << endl;

  //string satsuma2_path = std::getenv("SATSUMA2_PATH");
  //string current_path = std::getenv("PWD");

  string satsuma2_path  = getEnvVar("SATSUMA2_PATH");
  string current_path = getEnvVar("PWD");

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

  //if (bLSF && strstr(pExec, "/") == NULL) {
  //  cout << "When submitting to LSF, Satsuma must be run with the full path (e.g. '/usr/home/dummy/Satsuma...'" << endl;
  //  return -1;
  //}

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
  std::vector<string> queryNames, targetNames;

  for (i=0; i<target.size(); i++)
    targetNames.push_back(target.NameClean(i));
  for (i=0; i<query.size(); i++)
    queryNames.push_back(query.NameClean(i));
  matches.SetNames(queryNames, targetNames);


  for (i=0; i<target.size(); i++) {
    matches.SetTargetSize(i, target[i].size());
  }
  for (i=0; i<query.size(); i++)
    matches.SetQuerySize(i, query[i].size());
  cout << "Done." << endl;
  target.clear();
  query.clear();
  
  
  bool bNoChain = false;


  WorkQueue wq(minLen, sQuery, queryChunk, sTarget, targetChunk, minProb, sigCutoff, probtable, slave_count, threads_per_slave);
  //==================================================================
  //ALG: create filtered seeds
  cout << "SATSUMA: Acquiring seeds, date and time: " << GetTimeStatic() << endl;
  if (!bFilter) {
    if (old_seedFile != "") {
      cout<<"KMATCH not run because old seeds passed on."<<endl;
    }
    else if ( seedFile == "") {
      //Launch KMatch
      for (long long i=11; i<32; i+=2){
        seedFile = output + "/kmatch_results.k"+ to_string(i);
        string cmd;
        cmd = satsuma2_path + "/KMatch " + sQuery + " " + sTarget;
        cmd += " " + to_string(i) + " " + seedFile + " " + to_string(i) + " " + to_string(i-1) + " " + to_string((long long)max_kmatch_freq);
        cmd += "; touch " + seedFile + ".finished";

        cout << "Running seed pre-filter: " << endl << "  " << cmd << endl;

        // now build the command to call the shell script
        string sh_cmd;
        int ncpus = 2;  // kmatch uses max 2 threads
        sh_cmd = "sh " + satsuma2_path + "/satsuma_run.sh " + current_path + " \"" + cmd + "\" " + to_string(ncpus) + " " + to_string(kmatch_mem) + " KM" + to_string(i) + " " + to_string(bKm_sync);

        system(sh_cmd.c_str());
      }
      //wait for the kmatches to finish
      cout << "Waiting for seed pre-filters..." << endl;
      seedFile = output + "/kmatch_results.k";
      int finished=0;
      while (finished<11){
        sleep(1);
        for (long long i=11; i<32; i+=2){
          FILE * pProbe = fopen((seedFile+to_string(i)+".finished").c_str(), "r");
          if (pProbe != NULL) {
            fclose(pProbe);
            cout<<"loading results for k="<<i<<endl;
            wq.results_from_file((seedFile+to_string(i)).c_str(),grid.m_queryChunks);
            string cmd="rm "+seedFile+to_string(i)+".finished";
            system(cmd.c_str());
            finished++;
          }
        }
      }
      cout << "Seed pre-filters finished" << endl;
    } else {
      cout << "Loading pre-calculated seeds from "<<seedFile<<"*..."<<endl;
      for (long long i=11; i<32; i+=2){
        cout<<"loading results for k="<<i<<endl;
        wq.results_from_file((seedFile+to_string(i)).c_str(),grid.m_queryChunks);
      }
    }
  }
  
  
  cout << "SATSUMA: Launching slaves, date and time: " << GetTimeStatic() << endl;
  wq.setup_queue(slave_mem);   //XXX totally wrong to name the slaves count like that!
  //ALG: processing the seeds. Why is this different to just process matches? (besides the filter)
  
  if (old_seedFile != "") {
    cout<<"Reading old seeds from "<<old_seedFile<<endl;
    matches.Read(old_seedFile);
    matches.Sort();
    matches.Collapse();
  }
  else {
    cout << "SATSUMA: Processing KMATCH results, date and time: " << GetTimeStatic() << endl;
    unsigned long int new_matches_count;
    wq.collect_new_matches(matches); //matches.Read(seedFile);
    matches.Sort();
    matches.Collapse();
    matches.LengthFilter(min_seed_length);//XXX:hardcoded lenght filter-->terrible!!!
  }
  //=============================================================

  
  //MultiMatches workMatches;
  cout << "SATSUMA: Entering main search loop, date and time: " << GetTimeStatic() << endl;
  int targets_queue_size=slave_count*threads_per_slave*8;
  unsigned int main_iteration=0;
  bool first_pass=true;
  int extra_iterations=10;
  int exit_counter=extra_iterations;
  int min_slave_matches=perBlock*min_matches_per_target;//XXX: make this a parameter
  //ALG: main loop
  bool first_match_seen=false;
  TargetCoverageTracker coverage_tracker;
  while (true) {
    //TODO: ALG: collect matches (get count of tasks and results)
    main_iteration++;
    cout<<"MAIN: starting iteration "<<main_iteration<<endl;
    t_collect_status collect_status=wq.collect_new_matches(matches);
    cout<<"MAIN: "<<collect_status.matches<<" new matches collected"<<endl;
    
    //ALG: if new matches
    if (first_pass || collect_status.matches > 0) {
      first_pass=false;
      //ALG: chain matches
      cout << "MAIN: Running the chaining step..." << endl;
      MultiMatches chained;
      matches.Sort();
      matches.Collapse();
      cout << "MAIN: running coverage tracker" <<endl;
      std::pair<uint64_t,uint64_t> ctr=coverage_tracker.update_and_report(matches);
      cout << "MAIN: coverage tracker: +"<<ctr.first<<" -"<<ctr.second <<endl;


      if (bDumpCycleMatches){
        string cycle_out = output + "/cycle_" + to_string((unsigned long long) main_iteration) + ".matches"; 
        cout << "MAIN: dumping cycle matches..." << endl;
        matches.Write(cycle_out);
      }
      cout << "MAIN: DynProg..." << endl;
      if (bDup)
        RunMatchDynProgMult(chained, matches);
      else
        RunMatchDynProg(chained, matches);
      cout << "MAIN: Total matches: " << chained.GetMatchCount() << endl;
      //TODO:
      matches = chained;
      cout << "MAIN: Interpolating matches... " << chained.GetMatchCount() << endl;
      SyntenyInterpolator inter;
      inter.Interpolate(chained); 
      cout << "MAIN: Total matches (interpolated): " << chained.GetMatchCount() << endl;
      cout << "MAIN: Updating grid's Target Weights" << endl;
      grid.UpdateTargetWeights(chained);//Clears and updates weights
      
      
    }
    //TODO: ALG: collect targets to fill in the queue
    int targets_to_collect=(wq.pending_pair_count()<targets_queue_size ? targets_queue_size - wq.pending_pair_count() : 0 );
    if (targets_to_collect>targets_queue_size * .75){
      targets_queue_size+=threads_per_slave*slave_count;
      cout<<"MAIN: Increasing queue size to "<<targets_queue_size<<endl;
      targets_to_collect=(wq.pending_pair_count()<targets_queue_size ? targets_queue_size - wq.pending_pair_count() : 0 );
    }
    if (targets_to_collect<targets_queue_size * .50 && targets_queue_size>8*slave_count*threads_per_slave){
      targets_queue_size-=threads_per_slave*slave_count;
      cout<<"MAIN: Decreasing queue size to "<<targets_queue_size<<endl;
      targets_to_collect=(wq.pending_pair_count()<targets_queue_size ? targets_queue_size - wq.pending_pair_count() : 0 );
    }
    if (!collect_status.matches && !targets_to_collect) {
      cout<< "MAIN: nothing changed, skipping cycle and waiting 3 seconds"<<endl;
      sleep(3);
      main_iteration--;
      if (slave_count==0) {
        //This is only fo testing purposes (optimization and benchmarking)
        cout << "SATSUMA: forcing collect and termination (slave_count==0), date and time: " << GetTimeStatic() << endl;
        wq.close_queue();
        std::vector<GridTarget> newTargets;
        int realTargets = grid.CollectTargets(newTargets, 5000);
        wq.join();
        return 0;
      }
      continue;
    }
    cout << "MAIN: Collecting " << targets_to_collect << " new targets from the grid " << endl;
    std::vector<GridTarget> newTargets;
    int realTargets = grid.CollectTargets(newTargets, targets_to_collect);
    cout << "MAIN: Targets retrieved: " << newTargets.size() << endl;
    for (i=0; i<newTargets.size(); i++) {
        wq.add_pair(newTargets[i].TargetFirst(),
            newTargets[i].TargetLast(),
            newTargets[i].QueryFirst(),
            newTargets[i].QueryLast(),
            newTargets[i].IsFast());
    }
    //TODO: ALG: check CONVERGENCE condition
    cout<< "MAIN: STATUS: iteration, slaves_processed, new_matches, slave_to_match_relation, targets_needed, targets_collected"<<endl;
    cout<< "MAIN: STATUS: "<< main_iteration << ", " \
                           << collect_status.slaves << ", " \
                           << collect_status.matches << ", " \
                           << (collect_status.slaves ? ((double) collect_status.matches)/collect_status.slaves : 0 ) << ", " \
                           << targets_to_collect << ", " \
                           << newTargets.size() << ", " << GetTimeStatic() << endl;
    if (first_match_seen==false && collect_status.matches) {
      first_match_seen=true;
    }
    if (first_match_seen && collect_status.slaves){//XXX: extreme convergence, only useful to test functions
      if (collect_status.matches/collect_status.slaves<min_slave_matches) {
        exit_counter--;
        cout <<"MAIN: exit_counter="<<exit_counter<<endl;
        if (!exit_counter){
          break;
        }
      } else if (exit_counter!=extra_iterations){
        exit_counter=extra_iterations;
        cout <<"MAIN: exit_counter="<<exit_counter<<endl;
      }
    }



  }

  wq.close_queue(); //XXX: TODO: wait till all targets have been processed!

  //ALG: write final output
  cout << "Sorting and collapsing final matches..." << endl;
  matches.Sort();
  matches.Collapse();
  cout << "Writing final output!" << endl;
  string finalOut = output + "/xcorr_aligns.final.out"; 
  matches.Write(finalOut);


  //ALG: run ChainMatches
  /*string mergeCmd = satsuma2_path;

  mergeCmd += "/ChainMatches -i " + output + "/xcorr_aligns.almost.out";
  mergeCmd += " -o " + output + "/xcorr_aligns.final.out";

  cout << "Running " << mergeCmd << endl;
  system(mergeCmd.c_str());

  */
  //ALG: run /MergeXCorrMatches
  string mergeCmd = satsuma2_path;

  mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.final.out";
  mergeCmd += " -q " + sQuery + " -t " + sTarget;
  mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
  mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";

  cout << "Merging XCorr matches: " << endl << "  " << mergeCmd << endl;
  system(mergeCmd.c_str());

  //ALG: if bNoRef, refine matches
  if (bNoRef) {
    string refCmd = satsuma2_path;
    refCmd += "/HomologyByXCorr -cutoff 1.2  -q " + sQuery;
    refCmd += " -t " + sTarget;
    refCmd += " -o " + output + "/xcorr_aligns.refined.out";
    refCmd += " -guide " + output + "/xcorr_aligns.final.out";
    refCmd += " > " + output + "/HomologyByXCorr.refined.out";

    cout << "Running refined matches (fill in the blanks): " << endl << "  " << refCmd << endl;
    system(refCmd.c_str());


    refCmd = satsuma2_path;
    refCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.refined.out";
    refCmd += " -q " + sQuery + " -t " + sTarget;
    refCmd += " -o " + output + "/satsuma_summary.refined.out";
    refCmd += " > " + output + "/MergeXCorrMatches.refined.out";

    cout << "Merging XCorr matches: " << endl << "  " << refCmd << endl;
    system(refCmd.c_str());
  }

  cout << "Joining Workqueue thread"<<endl;
  wq.join();
  cout << "SATSUMA: all done, date and time: " << GetTimeStatic() << endl;

  return 0;
}

