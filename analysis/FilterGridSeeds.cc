#include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerAlignCore.h"
#include "analysis/SequenceMatch.h"

int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-t","target fasta file");
  commandArg<string> bStringCmmd("-q","query fasta file");
  commandArg<string> oStringCmmd("-o","output file");
  commandArg<int> distCmmd("-distance","distance between seeds", 1);
  commandArg<int> numCmmd("-w","width of the filter", 2);
  commandArg<bool> selfCmmd("-self","self-alignments", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Tool to test the VecDNAVector.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(distCmmd);
  P.registerArg(numCmmd);
  P.registerArg(selfCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);
  string output = P.GetStringValueFor(oStringCmmd);
  int distance = P.GetIntValueFor(distCmmd);
  int num12 = P.GetIntValueFor(numCmmd);
  bool bSelf = P.GetBoolValueFor(selfCmmd);

  if (bSelf) {
    cout << "ERROR: option -self is not implemented." << endl;
    return -1;
  }
 
  vecDNAVector target, query;
  
  target.Read(aString);
  query.Read(bString);
 

  int i, j, l;
  MultiMatches all;

  svec<string> queryNames, targetNames;

  for (i=0; i<target.isize(); i++) 
    targetNames.push_back(target.NameClean(i));
  for (i=0; i<query.isize(); i++) 
    queryNames.push_back(query.NameClean(i));
  all.SetNames(queryNames, targetNames);


  for (i=0; i<target.isize(); i++) {
    all.SetTargetSize(i, target[i].isize());
  }
  for (i=0; i<query.isize(); i++)
    all.SetQuerySize(i, query[i].isize());



  KmerAlignCore core;
  core.SetNumTables(num12);
  TranslateBasesToNumberExact trans;
  core.SetTranslator(&trans);
  
  cout << "Adding target." << endl;
  core.AddData(query);
  core.SortAll();
  cout << "done" << endl;
  int k = 12*num12;

  core.SetMax12Mer(150);

  for (i=0; i<target.isize(); i++) {
    const DNAVector & d = target[i];
    cout << target.Name(i) << endl;
    for (j=0; j<=d.isize()-k; j+=distance) {
      //cout << "len=" << d.isize() << " j=" << j << endl;
      DNAVector sub;
      sub.SetToSubOf(d, j, k);
      svec<KmerAlignCoreRecord> matches;
      //cout << "look" << endl;
      core.GetMatches(matches, sub);
      int found = 0;
      //cout << "=====" << endl;
      if (matches.isize() == 1) {
	for (l=0; l<matches.isize(); l++) {
	  //cout << target.NameClean(i) << "\t" << j << "\t";
	  //cout << query.NameClean(matches[l].GetContig()) << "\t";
	  //cout << matches[l].GetPosition() << "\t+" << endl;
	  bool bExtend = false;
	  if (all.GetMatchCount() > 0) {
	    SingleMatch & last = all.GetMatchDirect(all.GetMatchCount()-1);
	    if (last.GetQueryID() == matches[l].GetContig() && last.GetTargetID() == i &&
		last.GetStartTarget() == j-1 && last.GetStartQuery() == matches[l].GetPosition() - 1) {	     
	      last.SetPos(last.GetStartQuery(), last.GetStartTarget(), last.GetLength()+1, false);
	      bExtend = true;
	      //cout << "Merged (fw)" << endl;
	    }
	  }
	  if (!bExtend) {
	    SingleMatch m;
	    m.SetQueryTargetID(matches[l].GetContig(), i, query[matches[l].GetContig()].isize());
	    m.SetPos(matches[l].GetPosition(), j, k, false);
	    m.SetProbability(1.);
	    m.SetIdentity(1.);
	    m.AddMatches(k);
	    all.AddMatch(m);
	    found++;
	  }
	  
	}
      }
      matches.clear();
      sub.ReverseComplement();
      core.GetMatches(matches, sub);
      //cout << "Matches: " << matches.isize();
      if (matches.isize() == 1) {
	for (l=0; l<matches.isize(); l++) {
	  //cout << target.NameClean(i) << "\t" << j << "\t";
	  //cout << query.NameClean(matches[l].GetContig()) << "\t";
	  //cout << matches[l].GetPosition() << "\t-" << endl;
	  bool bExtend = false;
	  if (all.GetMatchCount() > 0) {
	    SingleMatch & last = all.GetMatchDirect(all.GetMatchCount()-1);

	    int newQPos = query[matches[l].GetContig()].isize() - k - matches[l].GetPosition();

	    //cout << "Checking rc, " << last.GetStartQuery() << " " <<  newQPos << endl;

	    if (last.GetQueryID() == matches[l].GetContig() && last.GetTargetID() == i &&
		last.GetStartTarget() == j-1 && last.GetStartQuery() == newQPos-1) {	     
	      last.SetPos(last.GetStartQuery(), last.GetStartTarget(), last.GetLength()+1, true);
	      //cout << "Merged (rc)" << endl;
	      bExtend = true;
	    }
	  }

	  if (!bExtend) {
	    SingleMatch m;
	    m.SetQueryTargetID(matches[l].GetContig(), i, query[matches[l].GetContig()].isize());
	    m.SetPos(query[matches[l].GetContig()].isize() - k - matches[l].GetPosition(), j, k, true);
	    //m.SetPos(matches[l].GetPosition(), j, k, true);
	    m.SetProbability(1.);
	    m.SetIdentity(1.);
	    m.AddMatches(k);
	    all.AddMatch(m);
	    found++;
	  }
	}
      }
      //cout << "===" << endl;
      //if (found == 0) {
      //cout << query.NameClean(i) << "\t" << j << "\tNO MATCH" << endl;	
      //}
    }
  }
  
  all.Write(output);
  
  return 0;

}
  
