//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif


#include "analysis/DNAVector.h"


#include <string.h>

#include "base/CommandLineParser.h"
#include "analysis/CrossCorr.h"
#include "analysis/SequenceMatch.h"
#include "analysis/XCorrDynProg.h"


//--------------------------------------------------------------------------------

void Load(vecDNAVector & out, svec<string> & names, const string & file, bool addRC) 
{
  /*
  FullNameParser fwParser;
  FastaSequenceFilestream fsf(file, &fwParser);

  veccompseq markerBases;
 
  int estSize = fsf.estimatedSize();
  longlong dataSize = fsf.estimatedData();
  names.reserve(estSize);
  markerBases.Reserve(dataSize, estSize);

  //cout << "Reading fasta..." << endl;
  fsf.parse(names, markerBases);
  //cout << "Building bases table from markers..." << endl;
      
  out.resize(markerBases.size());
  int i;
  for (i=0; i<(int)markerBases.size(); i++) {
    out[i] = markerBases[i].asBasevector();
    }*/

  out.Read(file, names);


  if (addRC) {
    int n = out.size();
    out.resize(2 * n);
    for (int i=0; i<n; i++) {
      out[i+n] = out[i];
      out[i+n].ReverseComplement();
    }
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
  commandArg<string> iStringCmmd("-i","XCorr input file (binary)");
  commandArg<string> oStringCmmd("-o","output file (summary)");
  commandArg<string> qltStringCmmd("-qlt","QLT-style output file", "");
  commandArg<int> lToffCmmd("-t_offset","offset of target sequence", 0);
  commandArg<int> lQoffCmmd("-q_offset","offset of query sequence", 0);
  commandArg<bool> bSwitch("-reverse","reverse query and target id's", false);
  commandArg<double> fMin("-min_score","minimum score to accept", 0.6);
  commandArg<bool> cmdSelf("-self","suppresses self-matches", false);


  commandLineParser P(argc,argv);
  P.SetDescription("Merges matches found by HomologyByXCorr");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(qltStringCmmd);
  P.registerArg(iStringCmmd);
  P.registerArg(lQoffCmmd);
  P.registerArg(lToffCmmd);
  P.registerArg(bSwitch);
  P.registerArg(fMin);
  P.registerArg(cmdSelf);

  P.parse();

  string sQuery = P.GetStringValueFor(aStringCmmd);
  string sTarget = P.GetStringValueFor(bStringCmmd);
  string input = P.GetStringValueFor(iStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  string qlt = P.GetStringValueFor(qltStringCmmd);
  int t_offset = P.GetIntValueFor(lToffCmmd);
  int q_offset = P.GetIntValueFor(lQoffCmmd);
  bool bReverse = P.GetBoolValueFor(bSwitch);
  bool bSelf = P.GetBoolValueFor(cmdSelf);
  double minAccept = P.GetDoubleValueFor(fMin);


  int i, j, k;

  vecDNAVector target, query;

  svec<string> targetNames, queryNames;

  Load(query, queryNames, sQuery, true);
  Load(target, targetNames, sTarget, false);

  int rcStart = (int)query.size() / 2;

  MultiMatches multi;
  MultiMatches raw;

  raw.Read(input);
  cout << "Done loading." << endl;


  cout << "Filtering." << endl;

  multi.SetCounts(raw.GetTargetCount(), raw.GetQueryCount());
  for (i=0; i<raw.GetTargetCount(); i++) {
    multi.SetTargetSize(i, raw.GetTargetSize(i));
    multi.SetTargetName(i, raw.GetTargetName(i));
  }
  for (i=0; i<raw.GetQueryCount(); i++) {
    multi.SetQuerySize(i, raw.GetQuerySize(i));
    multi.SetQueryName(i, raw.GetQueryName(i));
  }
  cout << "Matches before filtering: " << raw.GetMatchCount() << endl;
  for (i=0; i<raw.GetMatchCount(); i++) {
    const SingleMatch & m = raw.GetMatch(i);    
    if (m.GetProbability() > minAccept) {
      multi.AddMatch(m);
      //const SingleMatch & m2 = multi.GetMatch(multi.GetMatchCount()-1);
      //m2.Print();
    }
				  
  }
  cout << "       after filtering:  " << multi.GetMatchCount() << endl;

  cout << "Sorting..." << endl;
  multi.Sort();

  int total = 0;
  int lastEnd = 0;




  FILE * pOut = fopen(output.c_str(), "w");

  FILE * pQLT = NULL;
  if (qlt != "")
    pQLT = fopen(qlt.c_str(), "w");

  vecqualvector covered;
  covered.resize(target.size());
  for (i=0; i<target.size(); i++) {
    qualvector & q = covered[i];
    q.resize(target[i].size());
    for (j=0; j<(int)q.size(); j++)
      q[j] = 0;
  }

  cout << "Number of matches: " << multi.GetMatchCount() << endl;
  //  cout << "IMPORTANT NOTE: query coordinates for rc matches are on the rc of the query sequence!" << endl;

  for (i=0; i<multi.GetMatchCount(); i++) {
    const SingleMatch & m = multi.GetMatch(i);
    int start = m.GetStartTarget();
    if (m.GetTargetID() < 0) {
      cout << "ERROR at match " << i << endl;
      m.Print();
      continue;
    }
    qualvector & q = covered[m.GetTargetID()];

    //cout << m.GetTargetID() << " len=" << q.isize() << endl;

    if (start + m.GetLength() >= q.isize()) {
      cout << "ERROR: startT=" <<  start << " startQ=" << m.GetStartQuery() << " len=" << m.GetLength() << endl;
      continue;
    }

    if (m.GetLength() <= 0)
      cout << "ERROR!!" << endl;
    for (k=start; k<start + m.GetLength(); k++) {
      if (q[k] < 100)
	q[k]++;
    }
  }


  int allBases = 0;
  int coveredBases = 0;
  int depthCov = 0;

  for (i=0; i<covered.isize(); i++) {
    qualvector & q = covered[i];
    allBases += (int)q.size();
    for (j=0; j<(int)q.size(); j++) {
      depthCov += q[j];
      if (q[j] > 0)
	coveredBases++;
    }
  }


  cout << "Covered by at least one alignment: " << coveredBases << " (" << 100 * (double)coveredBases/(double)allBases << " %)" << endl;
  


  for (i=0; i<multi.GetMatchCount(); i++) {
    SingleMatch m = multi.GetMatch(i);
    int start = m.GetStartTarget();
    int startQ = m.GetStartQuery();


    //qualvector & q = covered[m.GetTargetID()];
    //for (k=start; k<m.GetStartTarget() + m.GetLength(); k++)
    //q[k]++;



    for (j=i+1; j<multi.GetMatchCount(); j++) {
      const SingleMatch & m2 = multi.GetMatch(j);
      //cout << "ms=" << m.GetStartTarget() << " ms2=" << m2.GetStartTarget() << " len=" <<  m.GetLength() << endl;
      
      
      if (m2.GetQueryID() == m.GetQueryID() && m2.IsRC() == m.IsRC() &&
	  m2.GetTargetID() == m.GetTargetID() &&
	  m2.GetStartTarget() < m.GetStartTarget() + m.GetLength() && 
	  m2.GetStartQuery() < m.GetStartQuery() + m.GetLength() &&
	  m2.GetStartTarget() >= m.GetStartTarget() &&
	  m2.GetStartQuery() >= m.GetStartQuery()) {
      
	//if (1 == 1) {
	int oldEndTarget =  m.GetStartTarget() + m.GetLength();
	m = m2;

	if (m.GetStartTarget() + m.GetLength() < oldEndTarget)
	  m.SetLength(oldEndTarget-m.GetStartTarget());

	continue;
      }
      break;
    }


    //const SingleMatch & mx1 = multi.GetMatch(i);
    //cout << "BEFORE MERGE!" << endl;
    //SeqMatch s;
    //s.Set(mx1.GetStartTarget(), mx1.GetStartQuery(), mx1.GetLength(), 0.);
    
    //if (mx1.IsRC()) {
    //  PrintMatch(query[mx1.GetQueryID()+rcStart], target[mx1.GetTargetID()], s);
    //} else {
    //  PrintMatch(query[mx1.GetQueryID()], target[mx1.GetTargetID()], s);
    //}

    const string & targetName = multi.GetTargetName(m.GetTargetID());
    const string & queryName = multi.GetQueryName(m.GetQueryID());
    int endQ = m.GetStartQuery() + m.GetLength();
    
    if (j > i+1) {
      cout << "Found overlapping matches, merging " << j-i << " blocks (" << i << " -> " << j << ")"  << endl;

      //cout << "Start=" << start << " m.StartTarget=" << m.GetStartTarget() << " m.Len=" << m.GetLength() << endl;
    }


    MatchDynProg dp(start, m.GetStartTarget() + m.GetLength() - start);
    int tID = m.GetTargetID();
    int qID = m.GetQueryID();
    int qIDOrig = qID;
    if (m.IsRC())
      qID += rcStart;
    for (k=i; k<j; k++) {
      dp.AddMatch(target[tID], query[qID], multi.GetMatch(k));
      //cout << "k=" << k << " start=" << multi.GetMatch(k).GetStartTarget() << " len=" << multi.GetMatch(k).GetLength() << endl;
                 
      //const SingleMatch & mx = multi.GetMatch(k);
      //cout << "BEFORE MERGE!" << endl;
      //SeqMatch s;
      //s.Set(mx.GetStartTarget(), mx.GetStartQuery(), mx.GetLength(), 0.);
      
      //if (mx.IsRC()) {
      //PrintMatch(query[mx.GetQueryID()+rcStart], target[mx.GetTargetID()], s);
      //} else {
      //PrintMatch(query[mx.GetQueryID()], target[mx.GetTargetID()], s);
      //}
      
      
    }
    svec<int> index;
    int len = dp.Merge(index);

    int qLen = (int)query[qID].size();
    
    if (m.IsRC()) {
      int tmpQ = startQ;
      startQ = qLen - endQ;
      endQ = qLen - tmpQ;
    }

    double ident = 0.; 
    if (!bSelf || targetName != queryName || m.IsRC() || t_offset + start != q_offset + startQ) {

      if (bReverse) {
	cout << "Query " << targetName << " ["<< t_offset + start << "-" << t_offset + m.GetStartTarget() + m.GetLength() << "] vs target " << queryName;
	cout << " [" << q_offset + startQ << "-" << q_offset + startQ << "]";
      } else {
	cout << "Query " << queryName << " ["<< q_offset + startQ << "-" << q_offset + endQ << "] vs target " << targetName;
	cout << " [" << t_offset + start << "-" << t_offset + m.GetStartTarget() + m.GetLength() << "]";
      }
      //cout << "Start: " << start << " end: " <<  m.GetStartTarget() + m.GetLength() << endl;
      
      if (m.IsRC())
	cout << " -";
      else 
	cout << " +";
      
      cout << " length " << m.GetStartTarget() + m.GetLength() - start << " check " << len << endl;
      
      ident = dp.PrettyPrint(index, start, target[tID], query[qID], len);

      total += m.GetStartTarget() + m.GetLength() - start;
    

      //if (m.GetStartTarget() < lastEnd) {
      //total -= lastEnd - m.GetStartTarget();
      //cout << "Subtracting for overlap: " << lastEnd - m.GetStartTarget() << endl;
      //}
      
      
      fprintf(pOut, "%s\t%d\t%d\t%s\t", targetName.c_str(), start, m.GetStartTarget() + m.GetLength(), queryName.c_str());
      fprintf(pOut, "%d\t%d\t%f\t", startQ, endQ, ident);
      if (m.IsRC()) {
	fprintf(pOut, "-\n");
      } else {
	fprintf(pOut, "+\n");
      }
    }

    if (pQLT != NULL) {
      fprintf(pQLT, "\n");
      int rc = 0;
      if (m.IsRC())
	rc = 1;

      int tlen = multi.GetTargetSize(tID);

      fprintf(pQLT, "QUERY\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t1\t0\t%d\t0\n", 
	      qIDOrig, startQ, endQ, endQ-startQ, rc, tID,
	      start, m.GetStartTarget() + m.GetLength(), tlen, endQ-startQ);
      
      fprintf(pQLT, "%d", qIDOrig);
      if (m.IsRC())
	fprintf(pQLT, "rc vs");
      else
	fprintf(pQLT, "fw vs");

      int mis = (int)((1. - ident) * (endQ-startQ)); 
      
      fprintf(pQLT, " %d, %d mismatches/0 indels (of %d), from %d-%d to %d-%d (of %d)\n", 
	      tID, mis, endQ-startQ, startQ, endQ, start, m.GetStartTarget() + m.GetLength(), tlen);
    }


    lastEnd = m.GetStartTarget() + m.GetLength();
    //cout << "Start: " << start << " end: " <<  m.GetStartTarget() + m.GetLength() << endl;
    
    i = j-1;
  }

  fclose(pOut);

  if (pQLT != NULL)
    fclose(pQLT);

  cout << endl;
  cout << "Total bases aligned:               " << total << endl;
  cout << "Total bases in target sequence:    " << allBases << endl;
  cout << "Covered by at least one alignment: " << coveredBases << " (" << 100 * (double)coveredBases/(double)allBases << " %)" << endl;
  cout << "All bases in aligns: " << depthCov << endl;

  return 0;
}

