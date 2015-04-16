#ifndef FORCE_DEBUG
#define NDEBUG
#endif




#include "analysis/DNAVector.h"


#include <string>

#include "base/CommandLineParser.h"
#include "analysis/CrossCorr.h"
#include "analysis/SequenceMatch.h"
#include "analysis/SeqChunk.h"
#include "analysis/AlignProbability.h"
#include "analysis/MatchDynProg.h"
#include "util/mutil.h"
#include "aligns/KmerAlignCore.h"



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



int SetGuideChunks(vecDNAVector & target, 
		   vecDNAVector & query, 
		   svec<SeqChunk> & targetInfo, 
		   svec<SeqChunk> & queryInfo, 
		   vecDNAVector & targetRaw, 
		   vecDNAVector & queryRaw, 
		   const svec<string> & targetNames, 
		   const svec<string> & queryNames, 
		   const MultiMatches & chained,
		   int size)
{
  int i, j;
  int maxLen = 100000;
  int minLen = 20;

  for (i=1; i<chained.GetMatchCount(); i++) {
    const SingleMatch & m = chained.GetMatch(i);
    const SingleMatch & n = chained.GetMatch(i-1);

    if (m.GetTargetID() != n.GetTargetID())
      continue;
    if (m.GetQueryID() != n.GetQueryID())
      continue;
    if (m.IsRC() != n.IsRC())
      continue;


    int forceOri = 1;
    if (m.IsRC())
      forceOri = -1;

    int len = n.GetLength();
    int lap = 32;          

    int startT = n.GetStartTarget() + len - lap;
    int startQ = n.GetStartQuery() + len - lap;

    if (startT < 0)
      startT = 0;
    if (startQ < 0)
      startQ = 0;

    int endT = m.GetStartTarget() + lap;
    int endQ = m.GetStartQuery() + lap;


    int targetSize = endT - startT;
    int querySize = endQ - startQ;
    

    if (targetSize > maxLen || querySize > maxLen)
      continue;
    if (targetSize < minLen || querySize < minLen)
      continue;

    int max = targetSize;
    if (querySize > max)
      max = querySize;

    int pieces = 1 + max / size;
    int tChunk = targetSize / pieces;
    int qChunk = querySize / pieces;

    //cout << "Chunking into " << pieces << " pieces (sizeT=" << tChunk << " sizeQ=" << qChunk << ")" << endl;
    //cout << " -> between target " << startT << "-" << endT << " query " << startQ << "-" << endQ << endl;


    int tID = m.GetTargetID();
    int qID = m.GetQueryID();
    const DNAVector & bT = targetRaw[tID];
    const DNAVector & bQ = queryRaw[qID];

    if (m.IsRC()) {
      startQ = (int)bQ.size() - endQ;
      endQ = startQ + querySize;
    }

    int tIter = startT;
    int qIter = startQ;

    while (qIter < endQ && tIter < endT) {
      DNAVector bQuery, bTarget;

      
      if (qIter + qChunk >= (int)bQ.size())
	break;
      if (tIter + tChunk >= (int)bT.size())
	break;
      
      if (qIter < 0) {
	cout << "Minor error (q)" << endl;
	qIter = 0;
      }
      if (tIter < 0) {
	cout << "Minor error (t)" << endl;
	tIter = 0;
      }

      bQuery.SetToSubOf(bQ, qIter, qChunk);
      bTarget.SetToSubOf(bT, tIter, tChunk);
      
      target.push_back(bTarget);
      query.push_back(bQuery);

      SeqChunk infoQ, infoT;

      infoQ.Set(queryNames[qID], qIter, qID);
      infoT.Set(targetNames[tID], tIter, tID);

      infoQ.ForceOrientation(forceOri);
      infoT.ForceOrientation(forceOri);

      targetInfo.push_back(infoT);
      queryInfo.push_back(infoQ);
      

      qIter += tChunk / 2;
      tIter += qChunk / 2;

    }     


  }
  return 0;
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

  bool IsValid(int iTarget, int iQuery) {
    if (m_targetStart.isize() == 0)
      return true;

    if (iTarget < m_tLo || iTarget > m_tHi)
      return false;
    if (iQuery < m_qLo || iQuery > m_qHi)
      return false;

    int i;
    for (i=0; i<m_targetStart.isize(); i++) {
      if (iTarget >= m_targetStart[i] && iTarget <= m_targetEnd[i] &&
	  iQuery >= m_queryStart[i] && iQuery <= m_queryEnd[i])
	return true;
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

  int m_tLo, m_tHi, m_qLo, m_qHi;
 

};


class PreFilter
{
public:
  PreFilter() {
    m_k = 12;
    m_trans.SetSize(m_k);
    m_table.resize(m_trans.GetBoundValue(), 0);
    m_core.SetTranslator(&m_trans);
    m_core.SetNumTables(2);
    m_core.SetMax12Mer(5);
  }

  void Setup(const DNAVector & d) {
    int i;
    for (i=0; i<=d.isize()-m_k; i++) {
      DNAVector tmp;
      tmp.SetToSubOf(d, i, m_k);
      int n = m_trans.BasesToNumber(tmp, 0);
      if (n > 0)
	m_table[n]++;
      tmp.ReverseComplement();
      n = m_trans.BasesToNumber(tmp, 0);
      if (n > 0)
	m_table[n]++;
    }
    vecDNAVector both;
    both.push_back(d);
    DNAVector rc = d;
    rc.ReverseComplement();
    both.push_back(rc);
    m_core.AddData(both);
    m_core.SortAll();
  }

  bool Valid(const DNAVector & d) {
    //return true;

    int i;
    /*
    for (i=0; i<=d.isize()-m_k*2; i++) {
      svec<KmerAlignCoreRecord> matches;
      m_core.GetMatches(matches, d, i);
      if (matches.isize() == 1)
	return true;
    }
    return false;
    */
      
    int count = 0;
    bool bLast = false;
    for (i=0; i<=d.isize()-m_k; i++) {
      int n = m_trans.BasesToNumber(d, i);
      if (m_table[n] > 0) {
	if (bLast)
	  count++;
	if (count >= 2*m_k-1)
	  return true;
	bLast = true;
      } else {
	bLast = false;
	count = 0;
      }
    }
    
    return false;
      
  }
private:
  TranslateBasesToNumberExact m_trans;
  KmerAlignCore m_core;

  int m_k;
  svec<short> m_table;
};

//================================================================
//================================================================
//================================================================
//================================================================
int main( int argc, char** argv )
{

#ifdef ALLOW_PROFILE
  RunTime();
#endif


  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<int> cqIntCmmd("-cq","query sequence id", -1);
  commandArg<int> ctIntCmmd("-ct","target sequence id", -1);
  commandArg<int> oStringCmmd("-o","output file (binary)");
  //commandArg<bool>  aBoolCmmd("-rc","reverse complement ONLY", false);
  commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
  commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
  commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);

  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.9999);
  commandArg<int> nBlocksCmmd("-nblocks","total number of blocks to process (when parallel)", 0);
  commandArg<int> blockCmmd("-block","block # to process (when parallel)", 0);

  commandArg<int> nBlocksQueryCmmd("-nblocks_query","total number of QUERY blocks to process (when parallel)", 0);
  commandArg<int> blockQueryCmmd("-block_query","QUERY block # to process (when parallel)", 0);

  commandArg<bool>  sameCmmd("-same_only","only align sequences that have the same name.", false);
  commandArg<string> gStringCmmd("-guide","XCorr input file (binary) to guide alignments", "");

  commandArg<bool> protCmmd("-proteins","translate seqeunces to proteins when aligning.", false);
  commandArg<double> cutoffCmmd("-cutoff","signal selection cutoff", 1.8);

  commandArg<bool> chainCmmd("-chain","chain matches before writing them out", false);
  commandArg<string> selectCmmd("-select","selection file", "");
  commandArg<string> pairsCmmd("-pairs","valid block pairs to process", "");
  commandArg<double> lineCmmd("-line","pick one block and compare to the query", -1.0);

  commandLineParser P(argc,argv);
  P.SetDescription("Compares two sequences via cross-correlation");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(cqIntCmmd);
  P.registerArg(ctIntCmmd);
  P.registerArg(lIntCmmd);
  //P.registerArg(aBoolCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(qChunkCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(nBlocksCmmd);
  P.registerArg(blockCmmd);
  P.registerArg(nBlocksQueryCmmd);
  P.registerArg(blockQueryCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(probCmmd);
  P.registerArg(gStringCmmd);
  P.registerArg(protCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(chainCmmd);
  P.registerArg(selectCmmd);
  P.registerArg(pairsCmmd);
  P.registerArg(lineCmmd);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string guideFile = P.GetStringValueFor(gStringCmmd);
  int minLen = P.GetIntValueFor(lIntCmmd);
  int contigQ = P.GetIntValueFor(cqIntCmmd);
  int contigT = P.GetIntValueFor(ctIntCmmd);
  //bool rc2 = P.GetBoolValueFor(aBoolCmmd);
  int targetChunk = P.GetIntValueFor(tChunkCmmd);
  int queryChunk = P.GetIntValueFor(qChunkCmmd);
  int nBlocks = P.GetIntValueFor(nBlocksCmmd);
  int myBlock = P.GetIntValueFor(blockCmmd);
  int nBlocksQuery = P.GetIntValueFor(nBlocksQueryCmmd);
  int myBlockQuery = P.GetIntValueFor(blockQueryCmmd);
  bool bSameOnly = P.GetBoolValueFor(sameCmmd);
  double minProb = P.GetDoubleValueFor(probCmmd);

  bool bProtein = P.GetBoolValueFor(protCmmd);
  double topCutoff = P.GetDoubleValueFor(cutoffCmmd);
  bool bChain = P.GetBoolValueFor(chainCmmd);
  string selection = P.GetStringValueFor(selectCmmd);
  string thePairs = P.GetStringValueFor(pairsCmmd);
  double oneLine = P.GetDoubleValueFor(lineCmmd);

  int i, j, k;

  if (output == "") {
    cout << "You MUST specify a valid output file (option '-o')!" << endl;
    return -1;
  }
  cout << "Command line:" << endl;
  for (i=0; i<argc; i++) {
    cout << argv[i] << " ";
  }
  cout << endl << endl;


  if (guideFile != "" && nBlocks != 0) {
    cout << "Guided mode, chunking not available!!" << endl;
    return -1;
  }

  vecDNAVector target, query;
  vecDNAVector targetRaw, queryRaw;

  svec<string> targetNames, queryNames;

  Load(queryRaw, queryNames, sQuery);
  Load(targetRaw, targetNames, sTarget);


  cout << "Query sequence:  " << sQuery << endl;
  cout << "Target sequence: " << sTarget << endl;

  cout << "Done loading." << endl;


  if (contigT != -1) {
    cout << "Erasing contigs..." << endl;
  }
     

  ValidChunkPairs validPairs;
  if (thePairs != "")
    validPairs.LoadFromString(thePairs);
  

  ChunkManager cmTarget(targetChunk, targetChunk / 2);
  ChunkManager cmQuery(queryChunk, 0);

  svec<SeqChunk> targetInfo, queryInfo;

  if (selection != "") {
    cout << "Initializing from selection file..." << endl;
    SeqChunkSelect sel;
    sel.Read(selection);
    cmTarget.ChunkItSelect(target, targetInfo, targetRaw, targetNames, sel.Target(), nBlocks, myBlock);
    cmQuery.ChunkItSelect(query, queryInfo, queryRaw, queryNames, sel.Query(), nBlocksQuery, myBlockQuery);
  } else {
    cmTarget.ChunkIt(target, targetInfo, targetRaw, targetNames, nBlocks, myBlock, oneLine);
    cmQuery.ChunkIt(query, queryInfo, queryRaw, queryNames, nBlocksQuery, myBlockQuery);
  }

  bool bOneOnOne = false;
  MultiMatches chained;      

  if (guideFile != "") {
    bOneOnOne = true;
    cout << "Reading " << guideFile << endl;
    chained.Read(guideFile);
    cout << "Done loading, recomputing chunks." << endl;
    target.clear();
    query.clear();
    targetInfo.clear();
    queryInfo.clear();

    SetGuideChunks(target, 
		   query, 
		   targetInfo, 
		   queryInfo, 
		   targetRaw, 
		   queryRaw, 
		   targetNames, 
		   queryNames, 
		   chained,
		   targetChunk);
  }

  queryRaw.clear();
  targetRaw.clear();
  

  cout << "Done chunking." << endl;
  
  //if (rc2) { 
  //cout << "rc'ing the query sequence." << endl;
  //query[contigQ].ReverseComplement();
  //}

  HomologyByXCorr hx;

  hx.SetMinProb(minProb);
  cout << "Keeping alignments more like real than " << minProb << endl;
  hx.SetTopCutoff(topCutoff);

  MultiMatches multi;
  multi.SetNames(queryNames, targetNames);

  int skipped = 0;
  int accepted = 0;

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

  if (guideFile != "") {
    cout << "Using target size (guided) " << targetChunk << endl;
    hx.SetTargetSize(targetChunk);
  }

  CCSignal * pQuerySig;
  CCSignal * pTargetSig;

  if (!bProtein) {
    pQuerySig = new CCSignal;
    pTargetSig = new CCSignal;
  } else {
    pQuerySig = new CCSignalWithCodons;
    pTargetSig = new CCSignalWithCodons;
  }


  PreFilter filter;

  for (j=0; j<(int)target.size(); j++) {
    
    if (target[j].isize() == 0)
      continue;


    if (!validPairs.IsValidTarget(j))
      continue;

    if (oneLine >= 0.) {
      cout << "Num targets: " << target.size() << endl;
      filter.Setup(target[j]);
    }

    CCSignal & targetSignal = *pTargetSig;
    CCSignal & querySignal = *pQuerySig;  
    targetSignal.SetSequence(target[j], targetChunk * 2);

    const SeqChunk & tChunk = targetInfo[j];

    cout << "Target chunk: " << j << endl;

    for (i=0; i<(int)query.size(); i++) {

      if (i > 0 && i % 1000000 == 0) {
	cout << "Processed " << i << " query chunks (of " << query.size() << ")" << " skip=" << skipped << " accepted=" << accepted << endl;
      }
      if (query[i].isize() == 0)
	continue;
      if (oneLine >= 0.) {
	if (!filter.Valid(query[i])) {
	  //cout << "Skip." << endl;
	  skipped++;
	  continue;
	} else {
	  accepted++;
	  //cout << "Passed" << endl;
	}
      }

      if (!validPairs.IsValid(j, i))
	continue;
      

      const SeqChunk & qChunk = queryInfo[i];

 
      
      if (bSameOnly) {
	if (targetNames[tChunk.GetID()] != queryNames[qChunk.GetID()])
	  continue;
      }

      if (bOneOnOne) {
	int r = i - j;
	if (r < -3 || r > 3)
	  continue;
      }
      


      int offset = qChunk.GetStart();
      hx.SetStuff(qChunk.GetID(),
		  tChunk.GetID(),
		  tChunk.GetStart(), 
		  offset, 
		  (int)query[i].size(), 
		  cmQuery.GetSize(qChunk.GetID()));


      bool bDoForward = true;
      bool bDoBackwards = true;
      if (qChunk.GetForceOrientation() == -1)
	bDoForward = false;
      if (qChunk.GetForceOrientation() == 1)
	bDoBackwards = false;

      if (bDoForward) {
	querySignal.SetSequence(query[i], targetChunk * 2);           
	hx.Align(querySignal, targetSignal, query[i], target[j], i);
      }
      
      if (bDoBackwards) {
	DNAVector rcQuery = query[i];
      
	rcQuery.ReverseComplement();
	querySignal.SetSequence(rcQuery, targetSignal.GetFullSize());
     
	hx.Align(querySignal, targetSignal, rcQuery, target[j], i, true);
      }
      // Let's put it back...
      //query[i].ReverseComplement();
      
    }
  

    //hx.CoverageStats();
    //hx.ClearCoverage();
  }

  if (guideFile != "") {
    cout << "Merging w/ guide..." << endl;
    multi.MergeRead(guideFile);
    multi.Sort();
  }


  if (bChain) {
    MultiMatches out;
    multi.Sort();
    RunMatchDynProg(out, multi);
    cout << "Saving alignments..." << endl;
    out.Write(output);
    cout << "all done!" << endl;
    
  } else {
    cout << "Saving alignments..." << endl;
    multi.Write(output);
    cout << "all done!" << endl;
  }


  delete pQuerySig;
  delete pTargetSig;

  return 0;
}

