#include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerDynProg.h"
//#include "src/SequenceMatch.h"

int main(int argc,char** argv)
{

  
  commandArg<string> targetCmmd("-t","target fasta file");
  commandArg<string> queryCmmd("-q","query fasta file");
  commandArg<string> outCmmd("-o","output file");
  commandArg<bool> quietCmmd("-quiet","minimal output to the screen", false);
  commandArg<int> numCmmd("-k","number of 12-mers used for seeding", 2);
  commandArg<int> laCmmd("-la","number of lookahead k-mers", 12);
  commandArg<int> skipCmmd("-la_skip","skip k-mers in lookahead", 2);
  commandArg<int> wordCmmd("-w","size of one k-mer", 12);

  commandArg<int> blockCmmd("-block","process sequence block", 0);
  commandArg<int> n_blockCmmd("-n_block","total blocks (processes)", 0);
  commandLineParser P(argc,argv);

  commandArg<bool> selfCmmd("-self","ignore matches (query=target)", false);
  commandArg<bool> repCmmd("-noreps","do not use soft-masekd repeats", false);

  P.SetDescription("Welcome to Papaya, alignments with seeds!");
  P.registerArg(targetCmmd);
  P.registerArg(queryCmmd);
  P.registerArg(outCmmd);
  P.registerArg(quietCmmd);
  P.registerArg(numCmmd);
  P.registerArg(laCmmd);
  P.registerArg(skipCmmd);
  P.registerArg(selfCmmd);
  P.registerArg(repCmmd);
  P.registerArg(wordCmmd);

  P.registerArg(blockCmmd);
  P.registerArg(n_blockCmmd);

  P.parse();

  string targetName = P.GetStringValueFor(targetCmmd);
  string queryName = P.GetStringValueFor(queryCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  bool bQuiet = P.GetBoolValueFor(quietCmmd);
  int num_kmers = P.GetIntValueFor(numCmmd);
  int look_ahead = P.GetIntValueFor(laCmmd);
  int look_ahead_skip = P.GetIntValueFor(skipCmmd);

  int myBlock = P.GetIntValueFor(blockCmmd);
  int nBlocks = P.GetIntValueFor(n_blockCmmd);
  int wordSize = P.GetIntValueFor(wordCmmd);

  bool bSelf = P.GetBoolValueFor(selfCmmd);
  bool bNoReps = P.GetBoolValueFor(repCmmd);
  
  svec<string> queryNames, targetNames;

  vecDNAVector query;
  vecDNAVector target;
  
  if (bNoReps) {
    cout << "Reading query (masked)..." << endl;
    query.Read(queryName, false, false, false);
    cout << "Reading target (masked)..." << endl;
    target.Read(targetName, false, false, false);
  } else {
    cout << "Reading query..." << endl;
    query.Read(queryName, queryNames);
    cout << "Reading target..." << endl;
    target.Read(targetName, targetNames);
  }
  cout << "done!" << endl;
    
  //test.ReverseComplement();


  int i, j;


  KmerSuperAligner sa;

  sa.SetWordSize(wordSize);
  sa.SetNumKmers(num_kmers);
  sa.SetLookAhead(0);
  sa.SetNewLookahead(look_ahead, look_ahead_skip);

  sa.SetRefBases(target);

  FILE * pOut = fopen(outName.c_str(), "w");

  int from = 0;
  int to = query.isize();

  if (nBlocks > 0) {
    int blockSize = query.isize() / nBlocks;
    from = myBlock * blockSize;
    to = (myBlock + 1) * blockSize;
    if (to >= query.isize())
      to = query.isize();
    if (query.isize() - to < blockSize)
      to = query.isize();
    cout << "Query sequences: " << query.isize() << " will do seqs " << from << " through " << to-1 << endl; 
  }


  for (i=from; i<to; i++) {
    
    vecDNAVector tmpBases;
    // Slightly stupid for now
    tmpBases.push_back(query[i]);
    svec<int> contigStarts;
    svec<int> contigDevs;
    contigStarts.push_back(0);
    contigDevs.push_back(0);

    SuperAlign result;
    int skip = -1;

    if (bSelf)
      skip = i;
    sa.Align(result, tmpBases, contigStarts, contigDevs, skip);


    //cout << "Matches: " << result.GetMatchCount() << endl;
    for (j=0; j<result.GetMatchCount(); j++) {
      const SuperMatch & m = result.GetMatch(j);
      int contig = m.GetContig();
      int startQuery = m.GetFirstBase();
      int endQuery =  m.GetLastBase();


      int targetID = m.GetRefID();
      bool targetRC = m.GetRC();
      int startTarget = m.GetRefStart();
      int endTarget = m.GetRefEnd();
      int len = m.GetLastBase() - m.GetFirstBase() + 1;
      double ident = 1.;
      //cout << "c=" << contig << " t=" << targetID << endl;

      fprintf(pOut, "%s\t%d\t%d\t%s\t", target.NameClean(targetID), startTarget, endTarget, query.NameClean(i));
      fprintf(pOut, "%d\t%d\t%f\t", startQuery, endQuery, ident);
      //fprintf(pOut, "%s\t%d\t%d\t%s\t", query.NameClean(i), startQuery, endQuery, target.NameClean(targetID));
      //fprintf(pOut, "%d\t%d\t%f\t", startTarget, endTarget, ident);
      if (m.GetRC()) {
	fprintf(pOut, "-");
      } else {
	fprintf(pOut, "+");
      }
      fprintf(pOut, "\t%d\n", query[i].isize());
    }

  }

  fclose(pOut);
  cout << "Closing." << endl;
  string doneName = outName + ".done";
  FILE * pDone = fopen(doneName.c_str(), "w");
  fprintf(pDone, "done\n");
  fclose(pDone);
  cout << "All done!" << endl;

  return 0;

}
  
