#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "DNAVector.h"

#include <string>
#include <unistd.h>
#include "../base/CommandLineParser.h"
#include "CrossCorr.h"
#include "SequenceMatch.h"
#include "SeqChunk.h"
#include "AlignProbability.h"
#include "MatchDynProg.h"
#include "../util/mutil.h"
#include "../util/SComm.h"
#include "../base/FileParser.h"



//==============================================================
int RCQuery(int offset, int start, int len, int chunkSize, int qLen) 
{
  int l = start + qLen - offset - chunkSize;
  return l;
}

void Load(vecDNAVector & out, svec<string> & names, const string & file) 
{
  out.Read(file, names);

}


class HomologyByXCorr
{
public:
  HomologyByXCorr() {
    m_minLen = 0;
    m_pMulti = NULL;
    m_offset = 0;
    m_qLen = 0;
    m_chunk = 0;
    m_tOffset = 0;
    m_qID = -1;
    m_tID = -1;

    m_targetSize = 0;
    m_minProb = 0.99;
  }

  void SetTopCutoff(double c) {
    m_sa.SetTopCutoff(c);
  }

  void SetMinProb(double p) {
    m_minProb = p;
  }

  void SetStuff(int qID,
		int tID,
		int tOffset, 
		int qOffset, 
		int chunkSize, 
		int queryLen) {
    m_qID = qID;
    m_tID = tID;
    m_offset = qOffset;
    m_tOffset = tOffset;
    m_qLen = queryLen;
    m_chunk = chunkSize;
  }

  void Align(const CCSignal & querySignal, 
	     const CCSignal & targetSignal, 
	     const DNAVector & query,
	     const DNAVector & target,
	     int id = 0,
	     bool rc = false);

  void ClearCoverage() {
    m_bestID.clear();
    m_ident.clear();
  }

  void SetMultiMatches(MultiMatches * p) {
    m_pMulti = p;
  }

  void CoverageStats() {
    int i;
    int k = 0;
    double avg = 0.;
    for (i=0; i<m_bestID.isize(); i++) {
      if (m_bestID[i] != -1) {
	k++;
	avg += m_ident[i];
      }
    }
    if (k == 0) {
      cout << "No hits." << endl;
      return;
    }
    avg /= (double)k;
    cout << "Bases covered: " << k << " (" << 100 * (double)k/(double)m_bestID.isize() << " %)" << endl; 
    cout << "Average identity: " << 100 * avg << " %" << endl;
  }

  void SetTargetSize(double i) {
    m_targetSize = i;
  }

  void SetMinimumAlignLen(int l) {
    m_minLen = l;
  }

  void SortFilter();


private:
  CrossCorrelation m_xc;
  SeqAnalyzer m_sa;
  int m_minLen;

  svec<int> m_bestID;
  svec<double> m_ident;
  MultiMatches * m_pMulti;

  int m_tID;
  int m_qID;
  int m_tOffset;
  int m_offset;
  int m_qLen;
  int m_chunk;
  double m_targetSize;
  double m_minProb;

  vecSeqMatch m_matches;

};


void HomologyByXCorr::SortFilter() 
{

}

void HomologyByXCorr::Align(const CCSignal & querySignal, 
			    const CCSignal & targetSignal,
			    const DNAVector & query,
			    const DNAVector & target,
			    int id,
			    bool rc)
{


  int i, j;
  svec<float> result;
  m_xc.CrossCorrelate(result, targetSignal, querySignal);
  m_matches.clear();
  m_sa.MatchUp(m_matches, query, target, result);

 

  for (j=0; j<m_matches.isize(); j++) {
    if (m_matches[j].GetLength() < m_minLen)
      continue;

    int len = m_matches[j].GetLength();
    int tStart = m_tOffset + m_matches[j].GetStartTarget();
    int qStart = m_offset +  m_matches[j].GetStartQuery();
    if (rc) {
      qStart = RCQuery(m_offset, m_matches[j].GetStartQuery(), len, m_chunk, m_qLen); 
    }

    double prob = GetMatchProbability(target, 
				      query, 
				      m_matches[j].GetStartTarget(), 
				      m_matches[j].GetStartQuery(), 
				      m_matches[j].GetLength(),
				      m_targetSize);

    if (prob < m_minProb)
      continue;

    //cout << "Match # " << j << " probability " << 100.* prob << " %" << endl;
    //cout << "Start target: " << tStart << " - " << tStart + len << endl;
    double ident = PrintMatch(query, target, m_matches[j], true);
    //double ident = PrintMatch(query, target, m_matches[j], false);
   


    SingleMatch m;
    m.SetQueryTargetID(m_qID, m_tID, m_qLen);
    m.SetPos(qStart, tStart, len, rc);
    m.SetProbability(prob);
    m.SetIdentity(ident);
    m.AddMatches(ident * (double)len);
    m_pMulti->AddMatch(m);

  }
}




class ValidChunkPairs
{
public:
  ValidChunkPairs() {
    m_tLo = m_tHi = m_qLo = m_qHi = -1;
  }

  bool IsValidTarget(int iTarget) {
    if (m_targetStart.isize() == 0)
      return true;
    if (iTarget < m_tLo || iTarget > m_tHi)
      return false;
    
    for (int i=0; i<m_targetStart.isize(); i++) {
      if (iTarget >= m_targetStart[i] && iTarget <= m_targetEnd[i])
	return true;
    }
    return false;
  }

  bool IsValid(int iTarget, int iQuery, bool & fast) {
    fast = false;
    if (m_targetStart.isize() == 0)
      return true;

    if (iTarget < m_tLo || iTarget > m_tHi)
      return false;
    if (iQuery < m_qLo || iQuery > m_qHi)
      return false;

    int i;
    for (i=0; i<m_targetStart.isize(); i++) {
      if (iTarget >= m_targetStart[i] && iTarget <= m_targetEnd[i] &&
	  iQuery >= m_queryStart[i] && iQuery <= m_queryEnd[i]) {
	fast = m_fast[i];
	return true;
      }
    }
    return false;

  }

  void LoadFromString(const string & s) {
    
    CMTokenizer tokenizer;
    tokenizer.AddDelimiter("-");
    tokenizer.AddDelimiter(",");
    tokenizer.AddDelimiter(";");
    tokenizer.AddDelimiter("]");
    tokenizer.AddDelimiter("[");

    CMPtrStringList l;
    tokenizer.Tokenize(l, (const char*)s.c_str());

    m_targetStart.clear();
    m_targetEnd.clear();
    m_queryStart.clear();
    m_queryEnd.clear();
    
    int i;
    for (i=0; i<l.length(); i+=4) {
      //cout << "i=" << (const char*)*l(i) << endl;
      if (*l(i) =="F") {
	i++;
	m_fast.push_back(true);
      } else {
	m_fast.push_back(false);
      }
      int tS = atol(*l(i));
      int tE = atol(*l(i+1));
      int qS = atol(*l(i+2));
      int qE = atol(*l(i+3));
      

      if (tE > m_tHi)
	m_tHi = tE;
      if (qE > m_qHi)
	m_qHi = qE;

      if (tS < m_tLo || m_tLo == -1)
	m_tLo = tS;
      if (qS < m_qLo || m_qLo == -1)
	m_qLo = qS;

      m_targetStart.push_back(tS);
      m_targetEnd.push_back(tE);
      m_queryStart.push_back(qS);
      m_queryEnd.push_back(qE);
    }
  }


private:
  svec<int> m_targetStart;
  svec<int> m_targetEnd;
  svec<int> m_queryStart;
  svec<int> m_queryEnd;
  svec<int> m_fast;

  int m_tLo, m_tHi, m_qLo, m_qHi;
 

};






//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{

  for (int xx=0; xx<argc; xx++)
    cout << argv[xx] << " ";


  commandArg<string> mStringCmmd("-master","name of the submit host");
  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<string> oStringCmmd("-o","output file (binary)", "");
  //commandArg<bool>  aBoolCmmd("-rc","reverse complement ONLY", false);
  commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
  commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
  commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);

  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.9999);
  commandArg<double> cutoffCmmd("-cutoff","signal selection cutoff", 1.8);
  commandArg<double> cutoffFastCmmd("-cutoff_fast","signal selection cutoff (fast)", 2.9);

  commandArg<string> pairsCmmd("-pairs","valid block pairs to process", "");

  commandLineParser P(argc,argv);
  P.SetDescription("Compares two sequences via cross-correlation");
  P.registerArg(mStringCmmd);
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(lIntCmmd);

  P.registerArg(qChunkCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(probCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(cutoffFastCmmd);
  P.registerArg(pairsCmmd);

  P.parse();

  string master = P.GetStringValueFor(mStringCmmd);
  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);

  double topCutoff = P.GetDoubleValueFor(cutoffCmmd);
  double topCutoffFast = P.GetDoubleValueFor(cutoffFastCmmd);
  string thePairs = P.GetStringValueFor(pairsCmmd);

  int i, j, k;

  //if (output == "") {
  //  cout << "You MUST specify a valid output file (option '-o')!" << endl;
  //  return -1;
  //}


  char out_dir[4096];
  strcpy(out_dir, output.c_str());
  for (i = strlen(out_dir)-1; i>=0; i--) {
    if (out_dir[i] == '/') {
      out_dir[i+1] = 0;
      break;
    }
  }


  cout << "Reporting to master host " << master << endl;

  cout << "Command line:" << endl;
  for (i=0; i<argc; i++) {
    cout << argv[i] << " ";
  }
  cout << endl << endl;

  vecDNAVector target, query;
  vecDNAVector targetRaw, queryRaw;

  svec<string> targetNames, queryNames;

  Load(queryRaw, queryNames, sQuery);
  Load(targetRaw, targetNames, sTarget);


  cout << "Query sequence:  " << sQuery << endl;
  cout << "Target sequence: " << sTarget << endl;

  cout << "Done loading." << endl;

  ValidChunkPairs validPairs;
  if (thePairs != "")
    validPairs.LoadFromString(thePairs);
  

  ChunkManager cmTarget(targetChunk, targetChunk / 4);
  ChunkManager cmQuery(queryChunk, 0);

  svec<SeqChunk> targetInfo, queryInfo;

  cmTarget.ChunkIt(target, targetInfo, targetRaw, targetNames, 0, 0);
  cmQuery.ChunkIt(query, queryInfo, queryRaw, queryNames, 0, 0);
  

  MultiMatches chained;      

  queryRaw.clear();
  targetRaw.clear();
  

  cout << "Done chunking." << endl;
  
 
  HomologyByXCorr hx;

  hx.SetMinProb(minProb);
  cout << "Keeping alignments more like real than " << minProb << endl;
  hx.SetTopCutoff(topCutoff);

  MultiMatches multi;
  multi.SetNames(queryNames, targetNames);



  double targetTotal = 0;
  for (i=0; i<cmTarget.GetCount(); i++) {
    multi.SetTargetSize(i, cmTarget.GetSize(i));
    targetTotal += (double)cmTarget.GetSize(i);
  }
  for (i=0; i<cmQuery.GetCount(); i++)
    multi.SetQuerySize(i, cmQuery.GetSize(i));


  hx.SetMultiMatches(&multi);

  hx.SetMinimumAlignLen(minLen);

  hx.SetTargetSize(targetTotal);

  CCSignal * pQuerySig;
  CCSignal * pTargetSig;

  pQuerySig = new CCSignal;
  pTargetSig = new CCSignal;


  // 30 sec timeout, then we'll exit
  int timeout = 2500;

  //------------------------- Main loop ---------------------------------------
  while (1==1) {

    if (output == "") {


      string die = out_dir;
      die += "/slaves.directive";
      cout << "Checking " << die << endl;
      FILE * pDie = fopen(die.c_str(), "r");
      if (pDie != NULL) {
	fclose(pDie);
	cout << "Exiting." << endl;
	return 0;
      }



      cout << "Listening to master..." << endl;
      
      bool bOK = false;
      char msg[2048];
      strcpy(msg, "");
      
      for (j=0; j<timeout; j++) {

	SCommReceiver * pRec = GetReceiver(master.c_str());

	if (pRec->Get(msg, sizeof(msg))) {
	  bOK = true;
	  break;
	}
	delete pRec;
	
	//if (pRec->IsFatal()) {
	//  cout << "COMMUNICATION ERROR: I don't feel so well, I'll better die now!!" << endl;
	//  return 0;
	//}
	cout << "Checking for directive..." << endl;
	FILE * pDie = fopen(die.c_str(), "r");
	if (pDie != NULL) {
	  fclose(pDie);
	  cout << "Exiting." << endl;
	  return 0;
	}
	

	//if (j % 10 == 0)
	cout << "Waiting..." << endl;
	sleep(1);
      }

      cout << "Got message: " << msg << endl;

 
      if (!bOK) {
	cout << "WARNING: exiting, master did not respond (timeout=" << timeout << "s)!" << endl;
	cout << "COMMUNICATION ERROR: I don't feel so well, I'll better die now!!" << endl;
	return 0;
	break;
      }
      
      if (strcmp(msg, "exit") == 0) {
	cout << "Received directive EXIT. Will exit." << endl;
	break;
      }
    

      StringParser sp;
      sp.SetLine(msg);
      
      output = sp.AsString(0);
      thePairs = sp.AsString(1);

      multi.ClearMatches();

    }
    cout << "Will process " << thePairs << " output to " << output << endl;

    validPairs.LoadFromString(thePairs);
    
    
    for (j=0; j<(int)target.size(); j++) {
      
      if (target[j].isize() == 0) {
	//cout << "Target size=0 " << j << endl;
	continue;
      }
      
      if (!validPairs.IsValidTarget(j)) {
	//cout << "Invalid target " << j << endl;
	continue;
      }
      
      CCSignal & targetSignal = *pTargetSig;
      CCSignal & querySignal = *pQuerySig;  
      targetSignal.SetSequence(target[j], targetChunk * 2);
      
      const SeqChunk & tChunk = targetInfo[j];
      
      cout << "Target chunk: " << j << endl;
      
      for (i=0; i<(int)query.size(); i++) {
	if (query[i].isize() == 0)
	  continue;
	
	bool bFast = false;
	if (!validPairs.IsValid(j, i, bFast))
	  continue;
	
	//cout << "Query: " << i << endl;

	if (bFast) {
	  hx.SetTopCutoff(topCutoffFast);
	} else {
	  hx.SetTopCutoff(topCutoff);
	}

	
	const SeqChunk & qChunk = queryInfo[i];
	
	if (i > 0 && i % 1000 == 0) {
	  cout << "Processed " << i << " query chunks (of " << query.size() << ")" << endl;
	}
	
	
	int offset = qChunk.GetStart();
	hx.SetStuff(qChunk.GetID(),
		    tChunk.GetID(),
		    tChunk.GetStart(), 
		    offset, 
		    (int)query[i].size(), 
		    cmQuery.GetSize(qChunk.GetID()));
	
	
	querySignal.SetSequence(query[i], targetChunk * 2);           
	hx.Align(querySignal, targetSignal, query[i], target[j], i);
	
	DNAVector rcQuery = query[i];
	
	rcQuery.ReverseComplement();
	querySignal.SetSequence(rcQuery, targetSignal.GetFullSize());
	
	hx.Align(querySignal, targetSignal, rcQuery, target[j], i, true);
	
      }
      
      
    }
    
    cout << "Saving alignments..." << endl;
    multi.Write(output);
    cout << "done for now" << endl;
    multi.ClearMatches();
    output = "";
    thePairs = "";

  }


  delete pQuerySig;
  delete pTargetSig;

  return 0;
}

