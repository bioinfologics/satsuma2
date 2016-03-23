#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Cola/NSGAaligner.h"
#include "src/Cola/NOIAligner.h"

class Stats
{
public:
  Stats() {
    m_match = 0;
    m_edits = 0;
  }

  void Add(const AlignmentCola& pAlign, int start, int stop) {
    AlignmentCola subAlign(pAlign);
    subAlign.keepSubalignment(start, stop);
    m_match += subAlign.getBaseMatched();
    m_edits += subAlign.calcEditCount();
  }

  void Print(const string & pre) {
    cout << "Stats: " << pre << " matches: " << m_match << " edits: " << m_edits << endl;
  }
private:
  int m_match;
  int m_edits;
};



void Fill(const DNAVector & t, svec<int> & same, const AlignmentCola& pAlign, int start, int stop)
{ 
  int i;
  
  for (i=start; i<stop; i++) {
    char c = pAlign.getQueryAlignCharForTarget(i);

    //if (i > 0 && pAlign->getQueryAlignIndexForTarget(i) - pAlign->getQueryAlignIndexForTarget(i-1) > 1
    //  && pAlign->getQueryAlignIndexForTarget(i) != -1
    //  && pAlign->getQueryAlignIndexForTarget(i-1) != -1) {
    if (i > start && i<stop-1) {
  if (pAlign.getQueryAlignIndexForTarget(i-1) == -1
      || pAlign.getQueryAlignIndexForTarget(i+1) == -1) {
    continue;
  }
      //insert[i-1] = 1;
    }
    if (c == t[i])
      same[i]++;     
  }
  
}

void PrintMotif(const string & pre, const DNAVector & d, const svec<int> & same, int count) 
{
  int len = 0;
  for (int i=0; i<same.isize(); i++) {
    int index = (10*(same[i]))/count;
    if (index < 9) {
      if (len > 5) {
  cout << pre << " Motif: ";
  for (int j=i-len; j<i; j++) {
    cout << d[j];
  }
  cout << endl;
      }
      len = 0;
      continue;
    }
    len++;
  }
}

int main(int argc,char** argv)
{


  commandArg<string> aStringCmmd("-i","list of multi-fasta files");
  
  //commandArg<string> bStringCmmd("-t","target sequence");
  //commandArg<bool> selfCmmd("-s","self-alignments", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Aligns sequences in a fasta file.");
  P.registerArg(aStringCmmd);
  //P.registerArg(bStringCmmd);
  //P.registerArg(selfCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  //string bString = P.GetStringValueFor(bStringCmmd);
  //bool bSelf = P.GetBoolValueFor(selfCmmd);
  
  svec<int> percentSWAff;
  percentSWAff.resize(11, 0);
  svec<int> percentCola;
  percentCola.resize(11, 0);
  svec<int> percentColaAff;
  percentColaAff.resize(11, 0);
  svec<int> percentSW;
  percentSW.resize(11, 0);

  FlatFileParser parser;

  int nESW = 0;
  int nECola = 0;
  int nESWAff = 0;
  int nEColaAff = 0;

  Stats sSW, sCola, sSWAff, sColaAff;

  
  parser.Open(aString);
  
  int nn = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    cout << "Processing: " << parser.AsString(0) << endl;

    nn++;
    if (nn % 10 == 0) {
      cout << "Prelim: " << endl;
      cout << "#\tSW\tCola\tSWAffi\tColaAff" << endl;
      for (int i=0; i<percentCola.isize(); i++) 
        cout << i << "\t" << percentSW[i] << "\t" << percentCola[i] << "\t" << percentSWAff[i] << "\t" << percentColaAff[i] << endl;
    }


    vecDNAVector seq;
    
    seq.Read(parser.AsString(0));

    if (seq.isize() < 20)
      continue;
    seq.RemoveGaps();
    
    if (seq[0].isize() > 1000 || seq[0].isize() < 100)
      continue;
    
    int i, j;
    
    
    const DNAVector & t = seq[0];
    svec<int> sameSW;
    svec<int> sameCola;
    svec<int> sameSWAff;
    svec<int> sameColaAff;
    sameSW.resize(t.isize(), 0);
    sameCola.resize(t.isize(), 0);
    sameSWAff.resize(t.isize(), 0);
    sameColaAff.resize(t.isize(), 0);
    svec<int> insert;
    insert.resize(t.isize(), 0);

    int count = 0;
    for (j=1; j<seq.isize(); j++) {
      //if (j > 1)
      //break;
      if (seq[j].isize() < 100)
        continue;

      cout << seq.Name(0) << " vs " << seq.Name(j) << endl;
      //        NSaligner aligner(target[i], query[j]);
      //        SWGAaligner aligner(target[i], query[j]);
      NSaligner alignerCola(t, seq[j]);
      //NSGAaligner aligner(t, seq[j], -12, -1, -2);
      NSGAaligner alignerColaAff(t, seq[j]);
      SWGAaligner alignerSWAff(t, seq[j]);
      SWGAaligner alignerSW(t, seq[j]);
      int start = 0;
      int stop = t.isize();
      
      cout << "Smith Waterman" << endl;
      alignerSW.align();
      cout << "Cola" << endl;
      alignerCola.align();
      cout << "Smith Waterman Affine" << endl;
      alignerSWAff.align();
      cout << "Cola Affine" << endl;
      alignerColaAff.align();

      //nESW += alignerSW.getAlignment()->getNumOfEdits();
      //nECola += alignerCola.getAlignment()->getNumOfEdits();
      //nESWAff += alignerSWAff.getAlignment()->getNumOfEdits();
      //nEColaAff += alignerColaAff.getAlignment()->getNumOfEdits();
  
      
      int s = alignerCola.getAlignment().getTargetOffset();
      int e = s + alignerCola.getAlignment().getTargetBaseAligned();
      //cout << 
      if (s > start)
        start = s;
      if (e < stop)
        stop = e;

      s = alignerSWAff.getAlignment().getTargetOffset();
      e = s + alignerSWAff.getAlignment().getTargetBaseAligned();
      if (s > start)
        start = s;
      if (e < stop)
        stop = e;
      s = alignerColaAff.getAlignment().getTargetOffset();
      e = s + alignerColaAff.getAlignment().getTargetBaseAligned();
      if (s > start)
        start = s;
      if (e < stop)
        stop = e;
      s = alignerSW.getAlignment().getTargetOffset();
      e = s + alignerSW.getAlignment().getTargetBaseAligned();
      if (s > start)
        start = s;
      if (e < stop)
        stop = e;

      if (stop - start < 150)
        continue;

      count++;
      Fill(t, sameCola, alignerCola.getAlignment(), start, stop);
            Fill(t, sameSWAff, alignerSWAff.getAlignment(), start, stop);
      Fill(t, sameColaAff, alignerColaAff.getAlignment(), start, stop);
      Fill(t, sameSW, alignerSW.getAlignment(), start, stop);

      sSW.Add(alignerSW.getAlignment(), start, stop);
      sSWAff.Add(alignerSWAff.getAlignment(), start, stop);
      sCola.Add(alignerCola.getAlignment(), start, stop);
      sColaAff.Add(alignerColaAff.getAlignment(), start, stop);
      
      
    }
    if (count < 20)
      continue;

    PrintMotif("SW", t, sameSW, count);
    PrintMotif("Cola", t, sameCola, count);
    PrintMotif("SWAff", t, sameSWAff, count);
    PrintMotif("ColaAff", t, sameColaAff, count);
    

    //cout << "Start=" << start << " Stop=" << stop << endl;
    for (i=0; i<t.isize(); i++) {
      //cout << "Conservation " << t[i] << "\t" << seq.isize() - 1 - diff[i];
      
      int index = (10*(sameCola[i]))/count;
      percentCola[index]++;
      index = (10*(sameSWAff[i]))/count;
      percentSWAff[index]++;
      index = (10*(sameColaAff[i]))/count;
      percentColaAff[index]++;
      index = (10*(sameSW[i]))/count;
      percentSW[index]++;
      
      //if (insert[i] > 0)
      //  cout << "\t---";
      //cout << endl;
    }
  }

  cout << "RESULT: " << endl;
  cout << "#\tSW\tCola\tSWAffi\tColaAff" << endl;
  for (int i=0; i<percentCola.isize(); i++) 
    cout << "histogram " << i << "\t" << percentSW[i] << "\t" << percentCola[i] << "\t" << percentSWAff[i] << "\t" << percentColaAff[i] << endl;

  //cout << "Number of edits:" << endl;
  //cout << "SW:\t" << nESW << endl;
  //cout << "Cola:\t" << nECola << endl;
  //cout << "SWAff:\t" << nESWAff << endl;
  //cout << "ColaAff:\t" << nEColaAff << endl;

  cout << endl;
  sSW.Print("Method SW");
  sCola.Print("Method Cola");
  sSWAff.Print("Method SWAff");
  sColaAff.Print("Method ColaAff");

  return 0;

}

