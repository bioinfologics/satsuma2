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
#include "util/SComm.h"
#include "util/SysTime.h"
#include <math.h>
#include "util/FindProcess.h"
#include <pthread.h>

void AddPair2String(string & s, int targetFrom, int targetTo, int queryFrom, int queryTo, bool bFast)
{
  char tmp[8192];
  

  if (bFast) {
    sprintf(tmp, "F[%d-%d,%d-%d]", targetFrom, targetTo, queryFrom, queryTo);
  } else {
    sprintf(tmp, "[%d-%d,%d-%d]", targetFrom, targetTo, queryFrom, queryTo);
  }
  s += tmp;
}


class JobSubmittor
{
public:
  JobSubmittor(const string & exec, const string & execReuse, const string & outDir, int targetBlocks, int queryBlocks, int maxJobs = 16) {
    m_exec = exec;
    m_execReuse = execReuse;
    m_outDir = outDir;
    m_maxJobs = maxJobs;
    m_num = 0;
    m_lastProcessed = 0;
    m_targetBlocks = targetBlocks;
    m_queryBlocks = queryBlocks;
    m_farmed = 0;
    //m_pTrans = GetTransmitter();
    m_pTrans = NULL;
  }

  ~JobSubmittor() {/*delete m_pTrans;*/}


  void SubmitSeeds(int nJobs, bool lsf = false, const string& advLSFArgs="", bool lsf_ini = false);
  void SubmitBlock(int firstTarget, int lastTarget, int fromQuery, int toQuery, bool lsf = false, const string& advLSFArgs="");
  void SubmitMultiBlock(const string & multi, bool lsf = false, const string& advLSFArgs="");
  void SubmitSelection(int firstTarget, int lastTarget, const SeqChunkSelect & sel, bool lsf = false, const string& advLSFArgs="");

  int Collect();

  MultiMatches & GetMulti() {return m_multi;}
  MultiMatches & GetMultiFinal() {return m_multiFinal;}
  int Last() const {return m_lastProcessed;}
  void Update() {m_lastProcessed = m_multi.GetMatchCount();}

  bool Done() {
    if (m_waitingList.isize() == 0 /*&& m_lastProcessed == m_multi.GetMatchCount()*/)
      return true;
    return false;
  }


  int HowManyRunning() const {return m_waitingList.isize();}

  void SetNumAlive(int n) {
    if (n < m_farmed) {
      cout << "Jobs died!! Will re-farm " << m_farmed - n << endl;
    }
    m_farmed = n;
  }

private:
  void Submit(const string & command, bool lsf, const string& advLSFArgs, bool bReally=true);

  string m_exec;
  string m_execReuse;
  string m_outDir;

  svec<string> m_waitingList;
  int m_maxJobs;
  int m_num;
  int m_lastProcessed;
  
  MultiMatches m_multi;
  MultiMatches m_multiFinal;

  int m_targetBlocks;
  int m_queryBlocks;

  int m_farmed;
  SCommTransmitter * m_pTrans;
  //ChunkAcceptor m_accept;
};



void JobSubmittor::SubmitSeeds(int nJobs, bool lsf, const string& advLSFArgs, bool lsfIni)
{

  int i;


  //int spacing = nTotal / nJobs;
  int tmpJobs = m_maxJobs;

  if (lsfIni) {
    m_maxJobs = nJobs;
    lsf = true;
  }

  for (i=0; i<nJobs; i++) {
    string command = m_exec;
    char jobs[8192];
    double line = (double)i/(double)nJobs;
    sprintf (jobs, " -line %f ", line);

    command += jobs;

    Submit(command, lsf, advLSFArgs);
    if (lsf) {
      sleep(2);
      if (i > 0 && i % 25 == 0)
	sleep(25);
    }
  }
  m_maxJobs = tmpJobs;
}

void JobSubmittor::SubmitMultiBlock(const string & s, bool lsf, const string& advLSFArgs)
{

  if (m_farmed < m_maxJobs) {
    string command = m_execReuse;
    
    command += " -pairs ";
    command += s;
    
    cout << "Farming out job for pairs " << s << endl;
    
    Submit(command, lsf, advLSFArgs);
    m_farmed++;
  } else {
    cout << "Re-feeding existing jobs " << s << endl;
    
    char num[8192];
    sprintf(num, "%d", m_num);
    m_num++;
      
    string output = m_outDir + "/xcorr_aligns." + num + ".data";

    string msg = output;
    msg += " ";
    msg += s;

    m_pTrans = GetTransmitter();
    m_pTrans->SendWait(msg.c_str());
    delete m_pTrans;
    m_pTrans = NULL;

    cout << "Job was being picked up!!" << endl;

    m_waitingList.push_back(output);
  }
}

void JobSubmittor::SubmitBlock(int firstTarget, int lastTarget, int firstQuery, 
                               int lastQuery, bool lsf, const string& advLSFArgs)
{
  string command = m_exec;
  char num[8192];
  sprintf(num, "%d", m_num);

  string select = m_outDir + "/select.";
  select += num;


  SeqChunkSelect sel;
  sel.Set(m_queryBlocks, m_targetBlocks);
  
  int i;

  cout << "SUBMIT Target " << firstTarget << "-" << lastTarget;
  cout << " Query " << firstQuery << "-" << lastQuery << endl;
  for (i=firstTarget; i<=lastTarget; i++)
    sel.SetTarget(i, 1);

  for (i=firstQuery; i<=lastQuery; i++)
    sel.SetQuery(i, 1);

  //cout << "Preparing selection..." << endl;
  sel.Write(select);

  command += " -select " + select;
  Submit(command, lsf, advLSFArgs);
}


void JobSubmittor::SubmitSelection(int firstTarget, int lastTarget, const SeqChunkSelect & selection, 
                                   bool lsf, const string& advLSFArgs)
{

  string command = m_exec;
  char num[8192];
  sprintf(num, "%d", m_num);

  string select = m_outDir + "/select.";
  select += num;

  SeqChunkSelect sel = selection;
  
  int i;

  cout << "SUBMIT SEED Target " << firstTarget << "-" << lastTarget << endl;
  for (i=firstTarget; i<=lastTarget; i++)
    sel.SetTarget(i, 1);

  //cout << "Preparing selection..." << endl;
  sel.Write(select);

  command += " -select " + select;
  Submit(command, lsf, advLSFArgs);

}


void RemoveFile(const string & file) {

  string cmd = "rm " + file;
  int ret = system(cmd.c_str());

}


int JobSubmittor::Collect()
{
  int i;

  int n = 0;

  for (i=0; i<m_waitingList.isize(); i++) {
    FILE * pTest = NULL;     
    pTest = fopen(m_waitingList[i].c_str(), "rb");
    if (pTest != NULL) {
      fclose(pTest);
      usleep(100000);
      //cout << "Reading output file " << m_waitingList[i] << endl;
      //int before = m_multi.GetMatchCount();
      MultiMatches candidate;
      candidate.Read(m_waitingList[i]);
      m_multiFinal.MergeRead(m_waitingList[i]);

      cout << "Reading/removing file " << m_waitingList[i] << endl;
    
      m_multi.Merge(candidate);
      cout << "Merged." << endl;
      RemoveFile(m_waitingList[i]);
      cout << "Removed." << endl;
      
      //int after = m_multi.GetMatchCount();
      //cout << "=> Read in data, matches: " << after - before << " " << m_waitingList[i] << endl;
      
      m_waitingList[i] = m_waitingList[m_waitingList.isize()-1];
      m_waitingList.resize(m_waitingList.isize()-1);
      i--;
      n++;
    }      

  }
  return n;
}



void JobSubmittor::Submit(const string & command, bool bLSF, const string& advLSFArgs, bool bReally)
{
  char num[8192];
  sprintf(num, "%d", m_num);
  m_num++;
  string submit;
  
  if (bLSF)
    submit = "bsub " + advLSFArgs + " -o " + m_outDir + "/HomologyByXCorr." + num + ".out ";
  
  submit += command;
  
  string output = m_outDir + "/xcorr_aligns." + num + ".data";

  submit += " -o " + output;

  if (!bLSF)
    submit += " > " + m_outDir + "/HomologyByXCorr." + num + ".out &";


  if (m_waitingList.isize() >= m_maxJobs) {
    cout << "Waiting before submitting jobs...." << endl;
    while (m_waitingList.isize() >= m_maxJobs) {
      usleep(100000);
      Collect();
    }
  }

  //cout << "Running: " << submit << endl;
 
  if (bReally) {
    cout << "Submitting " << submit << endl;
    int ret = system(submit.c_str());
  } else {
    cout << "NOT submitted!" << endl;
  }
  m_waitingList.push_back(output);

}


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
    cout << "Matches: " << last - first << endl;
    for (i=first+1; i<=last; i++) {
      const SingleMatch & m = chained.GetMatch(i);
      cout << m.GetQueryID() << " " << m.GetStartQuery();
      if (m.IsRC())
	cout << " -" << endl;
      else
	cout << " +" << endl;
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
    cout << "Slope: eps=" << eps << " x=" << sx << " y=" << sy << endl; 
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
	  cout << "Adding extra (left): " << xx << " / " << yy << endl;
	  

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
  commandArg<string> lsfAdvArgsCmmd("-lsfArgs","String containing arguments for running LSF - Ensure given string is correct", "");
  commandArg<bool> lsfCmmdIni("-lsf_ini","submit jobs to LSF", false);
  commandArg<int> perCmmd("-m","number of jobs per block", 32);
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
  commandArg<bool> exRefineCmmd("-noExhaustRefine","Do not run exhaustive alignment on synteny blocks and boundaries", false);
  commandArg<bool> selfCmmd("-self","self alignment", false);  
  commandArg<int> exRefineCoresCmmd("-exhaustRefineCores","Number of cores to use in running exhaustive alignment", 0);


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
  P.registerArg(lsfAdvArgsCmmd);
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
  P.registerArg(resumeCmmd);
  P.registerArg(seedCmmd);
  P.registerArg(blockPixelCmmd);
  P.registerArg(filterCmmd);
  P.registerArg(distanceCmmd);
  P.registerArg(dupCmmd);
  P.registerArg(seedWidthCmmd);
  P.registerArg(exRefineCmmd);
  P.registerArg(exRefineCoresCmmd);
  P.registerArg(selfCmmd);

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
  bool bLSF = P.GetBoolValueFor(lsfCmmd);
  string advLSFArgs = P.GetStringValueFor(lsfAdvArgsCmmd);
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
  bool noExRefine = P.GetBoolValueFor(exRefineCmmd);
  bool bSelf = P.GetBoolValueFor(selfCmmd);
  int exRefineCores = P.GetIntValueFor(exRefineCoresCmmd);
  
  if (nBlocksIni == -1)
    nBlocksIni = nBlocks;

 
  int i, j;

  cout << "SATSUMA: Welcome to SatsumaSynteny! Current date and time: " << GetTimeStatic() << endl;

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
    int ret = system(mk_cmd.c_str());
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
    if (bSelf) {
      cmd += " -self ";
    }

    cout << "Running seed pre-filter " << cmd << endl;
    int ret = system(cmd.c_str());
  }
 //==================================================================

  
  cout << "Checking for other running instances..." << endl;
  if (getProcessCount("SatsumaSynteny") > 1) {
    cout << "ERROR: An instance of SatsumaSynteny is already running!!!" << endl;
    exit(-1);
  }
  
  int nrun = getProcessCount("HomologyByXCorr");
  if (nrun > 0) {
    cout << "ERROR: " << nrun << " instance(s) of HomologyByXCorrSlave is already running!!!" << endl;
    exit(-1);
  }



  int exitCounter = 0;
  int superExitCounter = 0;
  int superExitCounterThresh = 20;
  int superExitCounterNot = 0;
  

  string simple = exec_dir;
  simple += "HomologyByXCorr -q " + sQuery;
  simple += " -t " + sTarget;
  char tmp[8192];
  sprintf(tmp, " -l %d -q_chunk %d -t_chunk %d -min_prob %f -cutoff %f ", minLen, queryChunkSeed, targetChunkSeed, minProb, sigCutoffSeed);

  simple += tmp;

  string simpleReuse = exec_dir;
  simpleReuse += "HomologyByXCorrSlave -q " + sQuery;
  simpleReuse += " -master ";
  GetHostName(tmp, sizeof(tmp));
  cout << "I am host " << tmp << endl;
  simpleReuse += tmp;
  simpleReuse += " -t " + sTarget;
  sprintf(tmp, " -l %d -q_chunk %d -t_chunk %d -min_prob %f -cutoff %f ", minLen, queryChunk, targetChunk, minProb, sigCutoff);
  simpleReuse += tmp;

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

  int nTargets = grid.NTargetChunks();
  
  cout << "Preparing job submission." << endl;

  JobSubmittor js(simple, simpleReuse, output, grid.NTargetChunks(), grid.NQueryChunks(), nBlocks);
  

  int collectCounter = 0;

  bool bNoChain = false;
  //=============================================================
  if (prevRun != "" || seedFile != "") {

    if (prevRun != "") {
      cout << "RESUMING from " << prevRun << endl;
    } else {
      cout << "Loading seeds from " << seedFile << endl;
      prevRun = seedFile;
    }

    MultiMatches & multi = js.GetMulti();
    MultiMatches & multiFin = js.GetMultiFinal();
    
    multi.Read(prevRun);
    //multiFin.Read(prevRun);

    MultiMatches tmp;
    multi.Sort();
    multi.Collapse();
    multiFin = multi;
    
    if (bDup)
      RunMatchDynProgMult(tmp, multi);
    else
      RunMatchDynProg(tmp, multi, bSelf);
    
    bNoChain = true;
    cout << "Done!" << endl;
    multi = tmp;
    multiFin = tmp;
    

      //else
      //chained = multi;



    
    if (seedFile == "") {

      /*
      for (i=0; i<multi.GetMatchCount(); i++) {
	const SingleMatch & m = multi.GetMatch(i);
	
	int x1 = m.GetStartTarget();
	int y1 = m.GetStartQuery();
	int x2 = m.GetStartTarget() + m.GetLength();
	int y2 = m.GetStartQuery() + m.GetLength();
	
	if (m.IsRC()) {
	  int s = multi.GetQuerySize(m.GetQueryID());
	  int tmp = y1;
	  y1 = s - y1;
	  y2 = s - y2;
	}
	
	grid.SetUsed(m.GetTargetID(), x1, x2, m.GetQueryID(), y1, y2);
      }
      */
    }

  } else {
    //=============================================================
    // Twice the # of seeds
    cout << "SATSUMA: Finding starting points, date and time: " << GetTimeStatic() << endl;
    js.SubmitSeeds(nBlocksIni, bLSF, advLSFArgs, bLSFIni);
  }

  int submitted = 0;

  int tmpRun = 0;
  int tmpRunFull = 0;

  int last = 0;




  
  //MultiMatches workMatches;
  cout << "SATSUMA: Entering main search loop, date and time: " << GetTimeStatic() << endl;

  cout << "Writing initial output..." << endl;
  string tmpOut1 = output + "/xcorr_aligns.init.out"; 
  js.GetMultiFinal().Write(tmpOut1);
  cout << "done!" << endl;

  MultiMatches debug;
  debug = js.GetMultiFinal();
  debug.ClearMatches();

  int counter = 0;
  int plus_counter = 0;


  while (true) {

    int l = 0;
    while (js.Collect() == 0) {
      usleep(100000);
      
      l++;
      if (l > 2)
  	break;
    }



    MultiMatches & multi = js.GetMulti();
    MultiMatches & multiPlus = js.GetMultiFinal();


    bool bSetUp = false;


    if (multi.GetMatchCount() - last > 0) {
      bSetUp = true;
      cout << "Running the chaining step..." << endl;
      //cout << "First, we are sorting the matches." << endl;
      MultiMatches chained;
      multi.Sort();

      //if (counter > 2)
      
      if (!bNoChain) {
	if (bDup)
	  RunMatchDynProgMult(chained, multi);
	else
	  RunMatchDynProg(chained, multi);
      } else {
	chained = multi;
	bNoChain = false;
      }
      //else
      //chained = multi;

      counter++;
      
      int startMatches = chained.GetMatchCount();
      cout << "Total matches: " << chained.GetMatchCount() << endl;
      multi = chained;

      SyntenyInterpolator inter;
      inter.Interpolate(chained); 
      cout << "Total matches (interpolated): " << chained.GetMatchCount() << endl;

      grid.ClearTargetWeights();

      
      int lastY1 = -1;
      int lastY2 = -1;
      for (i=0; i<chained.GetMatchCount(); i++) {
      //for (i=startMatches; i<chained.GetMatchCount(); i++) {
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

	//cout << "Considering." << endl;
	grid.ConsiderTargets(m.GetTargetID(), x1, x2, m.GetQueryID(), y1, y2, m.GetIdentity());

	/*
	if (i > 0) {
	  const SingleMatch & minus = chained.GetMatch(i-1);
	  if (minus.GetTargetID() == m.GetTargetID() &&
	      minus.GetQueryID() == m.GetQueryID() &&
	      minus.IsRC() == m.IsRC()) {

	    int mx1 = (x1 + minus.GetStartTarget())/2;
	    int my1 = (y1 + lastY1)/2;
	   
	    grid.ConsiderTargets(m.GetTargetID(), mx1, mx1 + m.GetLength()/2, m.GetQueryID(), my1, my1 + m.GetLength()/2);
	  }
	      
	  }*/
	lastY1 = y1;
	lastY2 = y2;
      }
      //??????????????????
    

      
      /*
      //===============================================================================
      if (plus_counter % 10 == 0) {
	
	cout <<  "Finding additional loners from the full set." << endl;
	multiPlus.Sort();

	for (i=1; i<chained.GetMatchCount()-1; i++) {
	  const SingleMatch & m = chained.GetMatch(i);
	  const SingleMatch & m_l = chained.GetMatch(i-1);
	  const SingleMatch & m_r = chained.GetMatch(i+1);
	  int dLeft = m.GetStartTarget() - (m_l.GetStartTarget() + m_l.GetLength());
	  int dRight = m_r.GetStartTarget() - (m.GetStartTarget() + m.GetLength());
	  
	  if (dLeft > 1000 && dRight > 1000) {

	    int x1 = m.GetStartTarget();
	    int y1 = m.GetStartQuery();
	    int x2 = m.GetStartTarget() + m.GetLength();
	    int y2 = m.GetStartQuery() + m.GetLength();
	    
	    if (m.IsRC()) {
	      int s = chained.GetQuerySize(m.GetQueryID());
	      int tmp = y1;
	      y1 = s - y1;
	      y2 = s - y2;
	    }
	    
	    grid.ConsiderTargets(m.GetTargetID(), x1, x2, m.GetQueryID(), y1, y2, m.GetIdentity());	  
 

	  }
	}

	}*/
      
      plus_counter++;
    }

    last = multi.GetMatchCount();

    svec<GridTarget> newTargets;
    collectCounter++;
    int seedPlus = 0;
    //if (collectCounter % 40 == 0 && nBlocks > 10)
    //seedPlus = 1;
    cout << "Collecting." << endl;
    int realTargets = grid.CollectTargets(newTargets, nBlocks * perBlock, seedPlus);
    cout << "Targets retrieved: " << newTargets.isize() << endl;


    //cout << "Collecting..." << endl;
    //int l = 0;
    //while (js.Collect() == 0) {
    //  usleep(100000);
      
    //  l++;
    //  if (l > 2)
    //	break;
    //}
    cout << "Ready to re-feed." << endl;



    //------------------------------------------------
    //for (i=0; i<newTargets.isize(); i++) {      
    //  if (newTargets[i].IsFast())
    //cout << "***** Target is FAST!!" << endl;
    //}
    //-----------------------------------------------
    int giveUpThresh = (nBlocks * perBlock)/4+1;



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
	cout << "Final collect!" << endl;
	string die = output + "/slaves.directive";
	FILE * pDie = fopen(die.c_str(), "w");
	fprintf(pDie, "exit\n");
	fclose(pDie);


	int counter = 0;
	while (js.HowManyRunning() > 0) {
	  js.Collect();
	  usleep(250000);
	  counter++;
	  if (counter > 500)
	    break;
	}
	cout << "All done, exiting!!" << endl;
	break;
      }
    } else {
      cout << "Resetting exit counter." << endl;
      exitCounter = 0;
    }


    // Check for alive processesss.

    int alive = getProcessCount("HomologyByXCorr");
    js.SetNumAlive(alive);
    
    for (i=0; i<newTargets.isize(); i++) {
      int to = i + perBlock;
      string theBlocks;
      if (to >= newTargets.isize())
	to = newTargets.isize();
      for (j=i; j<to; j++) {
      
	//if (newTargets[j].IsFast())
	//cout << "Target is FAST!!" << endl;

	AddPair2String(theBlocks,
		       newTargets[j].TargetFirst(),
		       newTargets[j].TargetLast(),
		       newTargets[j].QueryFirst(),
		       newTargets[j].QueryLast(),
		       newTargets[j].IsFast());

	SingleMatch dd;
	const SeqChunk & TCStart = grid.TargetChunk(newTargets[j].TargetFirst());
	const SeqChunk & TCEnd = grid.TargetChunk(newTargets[j].TargetLast());
	const SeqChunk & QCStart = grid.QueryChunk(newTargets[j].QueryFirst());
	const SeqChunk & QCEnd = grid.QueryChunk(newTargets[j].QueryFirst());

	
	int sTStart = TCStart.GetStart();
	int sTEnd = TCEnd.GetStart() + 4096;
	int sTID = TCStart.GetID();
	int sQStart = QCStart.GetStart();
	int sQEnd = QCEnd.GetStart() + 4096;
	int sQID = QCStart.GetID();

	dd.SetProbability(0.999);
	dd.SetIdentity(0.999);

	dd.SetQueryTargetID(sQID, sTID, debug.GetQuerySize(sQID));
	dd.SetPos(sQStart, sTStart, sTEnd - sTStart, false);
	debug.AddMatch(dd);

      }
      i = j-1;
      

      cout << "Multiblock " << theBlocks << endl;
      js.SubmitMultiBlock(theBlocks, bLSF, advLSFArgs);

      submitted++;
      tmpRun++;
      tmpRunFull++;
      cout << "Submitted: " << submitted << endl;
    }


  

    if (tmpRun > 20) {
      cout << "Writing temporary output..." << endl;
      string tmpOut = output + "/xcorr_aligns.temp.out"; 
      js.GetMulti().Write(tmpOut);
      
      cout << "Writing DEBUG output..." << endl;
      string debugOut = output + "/xcorr_aligns.debug.out"; 
      debug.Write(debugOut);

      tmpRun = 0;
    }

    if (tmpRunFull > 200) {
      cout << "Writing temporary full output..." << endl;
      string tmpOut = output + "/xcorr_aligns.temp.full.out"; 
      js.GetMultiFinal().Write(tmpOut);
      tmpRunFull = 0;
    }

    
    js.Update();
  }





  cout << "Writing final output!" << endl;
  string finalOut = output + "/xcorr_aligns.almost.out"; 
  js.GetMulti().Write(finalOut);
  finalOut = output + "/xcorr_aligns.final.all.out"; 
  js.GetMultiFinal().Write(finalOut);



  string mergeCmd = exec_dir;
  
  mergeCmd += "/ChainMatches -i " + output + "/xcorr_aligns.almost.out";
  mergeCmd += " -o " + output + "/xcorr_aligns.final.out";
  
  cout << "Running " << mergeCmd << endl;
  int ret = system(mergeCmd.c_str());


  mergeCmd = exec_dir;

  mergeCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.final.out";
  mergeCmd += " -q " + sQuery + " -t " + sTarget;
  mergeCmd += " -o " + output + "/satsuma_summary.chained.out";
  mergeCmd += " > " + output + "/MergeXCorrMatches.chained.out";
  
  cout << "Running " << mergeCmd << endl;
  ret = system(mergeCmd.c_str());
  
 
  if (bNoRef) {
    cout << "Running refined matches (fill in the blanks)." << endl;
    string refCmd = exec_dir;
    refCmd += "/HomologyByXCorr -cutoff 1.2  -q " + sQuery;
    refCmd += " -t " + sTarget;
    refCmd += " -o " + output + "/xcorr_aligns.refined.out";
    refCmd += " -guide " + output + "/xcorr_aligns.final.out";
    refCmd += " > " + output + "/HomologyByXCorr.refined.out";
    
    cout << "Running " << refCmd << endl;
    int ret = system(refCmd.c_str());
    
    
    refCmd = exec_dir;
    
    refCmd += "/MergeXCorrMatches -i " + output + "/xcorr_aligns.refined.out";
    refCmd += " -q " + sQuery + " -t " + sTarget;
    refCmd += " -o " + output + "/satsuma_summary.refined.out";
    refCmd += " > " + output + "/MergeXCorrMatches.refined.out";
    
    cout << "Running " << refCmd << endl;
    ret = system(refCmd.c_str());
  }


  if(!noExRefine) {
    cout << "Running exhaustive alignment for refining synteny blocks and gap boundaries" << endl;
    if (exRefineCores == 0)
      exRefineCores = nBlocks;
    string refCmd = exec_dir;
    refCmd += "/ColaAlignSatsuma ";
    refCmd += " -s " + output + "/satsuma_summary.chained.out";
    refCmd += " -q " + sQuery + " -t " + sTarget;
    refCmd += " -o " + output + "/summary_refinedAlignments.out";
    char tmp[128];
    sprintf(tmp, " -T %d ", exRefineCores);
    refCmd += tmp;

    cout << "Running " << refCmd << endl;
    int ret = system(refCmd.c_str());
  }


  cout << "all done!" << endl;


  cout << "SATSUMA: all done, date and time: " << GetTimeStatic() << endl;
  cout << "all done!" << endl;
  

  return 0;
}

