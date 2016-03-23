
#include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "src/DPAlign.h"
#include "analysis/SatsumaAlign.h"
#include <fstream>


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-q","query sequence");
  commandArg<string> bStringCmmd("-t","target sequence");
  commandArg<string> sStringCmmd("-s","papaya seeds file");
  commandArg<string> outFile("-o","papaya output");
  commandArg<string> outFileSummary("-summary","papaya output (summary file)", "");

  commandArg<double> identCmmd("-min","min identity", 0.999);
  commandArg<int> maxCmmd("-max_len","maximum alignment length", 8000);
  commandArg<int> maxGapCmmd("-max_gap","maximum gap", 25);
  commandArg<int> minCmmd("-min_len","minimum alignment length", 16);
  commandArg<string> tStringCmmd("-target","this target only", "");
  commandArg<bool> verboseCmmd("-verbose","pretty-print alignments", false);
  

  commandLineParser P(argc,argv);
  P.SetDescription("Aligns two fasta files.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(sStringCmmd);
  P.registerArg(tStringCmmd);
  P.registerArg(identCmmd);
  P.registerArg(maxCmmd);
  P.registerArg(maxGapCmmd);
  P.registerArg(outFile);
  P.registerArg(outFileSummary);
  P.registerArg(minCmmd);
  P.registerArg(verboseCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);

  string tString = P.GetStringValueFor(tStringCmmd);
  string sString = P.GetStringValueFor(sStringCmmd);
  double minIdent = P.GetDoubleValueFor(identCmmd);
  int maxLen = P.GetIntValueFor(maxCmmd);
  string out_file = P.GetStringValueFor(outFile);
  string out_file_summary = P.GetStringValueFor(outFileSummary);
  int minLen = P.GetIntValueFor(minCmmd);
  bool bVerbose = P.GetBoolValueFor(verboseCmmd);
  
  int maxGap = P.GetIntValueFor(maxGapCmmd);

  cout << "Reading satsuma file..." << endl;
  SAlignVec aligns;
  aligns.Read(sString);


  vecDNAVector query, target;
  
  cout << "Reading query " << aString << endl;
  query.Read(aString);
  cout << "Reading target " << bString << endl;
  target.Read(bString);

 //  int maxGap = 25;
 
  int i, j;

  //int off = 10;

  bool bLast = false;
  int printOff = 0;

  ofstream out_strm(out_file.c_str(), ios::out | ios::binary);
  //double minIdent = 0.99999;
  
  FILE * pSummary = NULL;
  if (out_file_summary != "")
    pSummary = fopen(out_file_summary.c_str(), "w");

  int alignCount = 0;
  for (i=0; i<aligns.isize(); i++) {
 
    const SAlign & one = aligns[i];

    for (j=i+1; j<aligns.isize(); j++) {
      const SAlign & two = aligns[j];

      //cout << "Considering target " << one.Target() << endl;    
      if (two.Target() != one.Target()) {
	break;
      }

      if (two.Query() != one.Query()) {
	break;
      }

      if (two.Direction() != one.Direction()) {
	break;
      }

      if (two.StartQuery() > aligns[j-1].EndQuery() + maxGap) {
	break;
      }

      if (one.Direction() == 1) {
	if (two.StartTarget() > aligns[j-1].EndTarget() + maxGap) {
	  break;
	}
	if (two.StartTarget() < aligns[j-1].StartTarget()) {
	  break;
	}
      } else {
	if (two.EndTarget() < aligns[j-1].StartTarget() - maxGap) {
	  break;
	}
	if (two.StartTarget() > aligns[j-1].StartTarget()) {
	  break;
	}
      }

    }
    const SAlign & two = aligns[j-1];

    // Align
    int startTarget = one.StartTarget();
    int endTarget = two.EndTarget();

    int startQuery = one.StartQuery();
    int endQuery = two.EndQuery();

    if (one.Direction() == 1) {
    } else {
      startTarget = two.StartTarget();
      endTarget = one.EndTarget();
    }

    i = j-1;

    if (endQuery - startQuery < minLen)
      continue;
    
 
    DNAVector t, q;
   
    string ugly = ">" + one.Target();
    string moreugly = ">" + one.Query();
    DNAVector & t_full = target(ugly);
    DNAVector & q_full = query(moreugly);


    t.SetToSubOf(t_full, startTarget, endTarget-startTarget+1);
    q.SetToSubOf(q_full, startQuery, endQuery-startQuery+1);

    if (one.Direction() == -1) {
      q.ReverseComplement();
    }

    if (!bLast && bVerbose) {
      
      cout << one.Target() << " [" << startTarget << "-" << endTarget << "] vs " << one.Query();
      cout << " [" << startQuery << "-" << endQuery << "] ";
      
      if (one.Direction() == -1)
	cout << "-" << endl;
      else
	cout << "+" << endl;
    }

    FullAlignment align;
    align.SetNameA(one.Target());
    align.SetNameB(one.Query());

    DPAligner aligner;
    aligner.SetMaxSkip(maxGap);
   
    aligner.UseRewardFunc(false);
    aligner.SetRect(t.isize(), q.isize());
    aligner.Align(align, t, q);

    if (one.Direction() == -1) {
      q.ReverseComplement();
      align.Flip(q.isize());
    }

    double ident = align.Identity(t, q);
    if (bVerbose) {
      align.PrettyPrint(t, q);
      cout << endl;
      cout << "Identity: " << ident << endl << endl;
    } else {
      if (i % 1000 == 0) {
	cout << "Blocks processed: " << i << " alignments found: " << alignCount << endl;
      }
    }
 
   
    if (ident >= minIdent) {
      //cout << startTarget << " -> " <<  endTarget << endl;
      for (int x=startTarget; x<= endTarget; x++) {
	
	t_full.SetQual(x, 1);
      }
    }

    if (pSummary != NULL) {
      fprintf(pSummary, "%s\t%d\t%d\t%s\t%d\t%d\t%f\t", 
	      one.Target().c_str(),
	      startTarget,
	      endTarget,
	      one.Query().c_str(),
	      startQuery,
	      endQuery,
	      ident);
      if (one.Direction() == -1)
	fprintf(pSummary, "-\n");
      else
	fprintf(pSummary, "+\n");
      fflush(pSummary);
    }




    alignCount++;
    align.SetOffsetA(startTarget);
    align.SetOffsetB(startQuery);
      
    align.WriteBinary(out_strm);

  }

  int numCov = 0;
  int numNuc = 0;
  for (i=0; i<target.isize(); i++) {
    const DNAVector & d = target[i];
    for (j=0; j<d.isize(); j++) {
      numNuc++;
      if (d.Qual(j) > 0)
	numCov++;
    }
  }

  cout << "Total nucleotides in target: " << numNuc << endl;
  cout << "Covered by at least 1 align at identity >= " << minIdent << ": " << numCov;
  cout << " (" << 100.*(double)numCov/(double)numNuc << "%)" << endl;

  out_strm.close();
  if (pSummary != NULL)
    fclose(pSummary);

  return 0;

}
  
