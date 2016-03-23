#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include<fstream>
#include "../src/Cola/NSaligner.h"
#include "../src/Cola/SWGAaligner.h"
#include "../src/Cola/NSGAaligner.h"
#include "../src/Cola/NOIAligner.h"
#include "RefineSatsuma.h"

//=====================================================================
void RefineSatsuma::init(const string& satsumaFile) { 
  FlatFileParser satsumaParser; /// parser object containing the satsuma alignment coords
  satsumaParser.Open(satsumaFile);
  blockElements.clear();
  AlignmentBlock curr;
  while (curr.parse(satsumaParser, false))
  {
    blockElements.push_back(curr);
  }
  //TODO fileparser is not returning when it should be - fix
  blockElements.resize(blockElements.size()-1);
}

void RefineSatsuma::alignAll(int mergeBoundary, int numOfThreads) const {
  int totSize   = blockElements.size();
  ThreadHandler th;
  if(numOfThreads>totSize) { numOfThreads = totSize; }
  for (int i=0; i<numOfThreads; i++) {
    char tmp[256];
    sprintf(tmp, "%d", i);
    string init = "init_";
    init += tmp;
    int from = i*(totSize/numOfThreads);
    int to   = (i+1)*(totSize/numOfThreads);
    if(i==(numOfThreads-1)) { to = totSize-1; } 
    th.AddThread(new AlignBlocksThread(*this, mergeBoundary, (realignFile+"."+tmp),
                                       (gapFile+"."+tmp), from, to, (i==numOfThreads-1), i));    
    th.Feed(i, init);
  }

  while (!th.AllDone()) {
    usleep(10000);
  }
}

void RefineSatsuma::alignSubset(int startIdx, int endIdx, int mergeBoundary, 
                                ofstream& realignFout, ofstream& gapFout, bool processLastBlock) const {
  AlignmentBlock last, curr, lastGap;
  for(int i=startIdx; i<=endIdx; i++) { 
    curr = blockElements[i];
    if (last.merge(curr, mergeBoundary)) {
      if(outputMode==0) {
        realignFout << "Merged blocks." << endl << last.getTargetChrom() << " : " 
                    << last.getTargetStart() << " - " << last.getTargetStop()-1;
      }
      lastGap = curr;
      if((i==endIdx)) {  //Last iteration : process here
        if(outputMode==0) { realignFout << "Re-aligning Satsuma block: " << endl; }
        alignBlock(last, realignFout); 
      }
      continue;
    }
    if(last.getTargetChrom() == "" || last.getTargetStart() < 0 || last.getQueryStart() < 0) {
      last = curr;
      lastGap = curr;
      if((i==endIdx)) { 
        if(outputMode==0) { realignFout << "Re-aligning Satsuma block: " << endl; }
        alignBlock(last, realignFout); 
      }
      continue;
    }
    if(outputMode==0) { realignFout << "Re-aligning Satsuma block: " << endl; }
    alignBlock(last, realignFout); 
    processGap(lastGap, curr, gapFout);
    last    = curr;
    lastGap = curr;
    if((i==endIdx) && (processLastBlock)) {  //Last iteration and last block should be processed
      if(outputMode==0) { realignFout << "Re-aligning Satsuma block: " << endl; }
      alignBlock(last, realignFout); 
    }
  }
}

void RefineSatsuma::alignBlock(const AlignmentBlock& aBlock, ofstream& realignFout) const {
  printRealignInfo(aBlock, realignFout);
  DNAVector q, t;
  q.SetToSubOf(query(aBlock.getQueryChrom()), aBlock.getQueryStart(),
     aBlock.getQueryStop()-aBlock.getQueryStart());
  t.SetToSubOf(target(aBlock.getTargetChrom()), aBlock.getTargetStart(),
     aBlock.getTargetStop()-aBlock.getTargetStart());
  if (aBlock.isReversed()) { q.ReverseComplement(); }
   q.SetName(aBlock.getQueryChrom());
   t.SetName(aBlock.getTargetChrom());

  realignFout << "##BLOCK" << endl;
  alignSeqs(t, aBlock.getTargetStart(), '+', q, aBlock.getQueryStart(), aBlock.getOrient(), realignFout);
}

void RefineSatsuma::alignSeqs(const DNAVector& t, int tOffset, char tStrand, 
                              const DNAVector& q, int qOffset, char qStrand, ostream& sout) const {
  Cola cola1 = Cola();
  AlignmentCola cAlign = cola1.createAlignment(q, t, aligner);
  cAlign.setSeqAuxInfo(qOffset, tOffset, qStrand, tStrand); //Note target and query mean the opposite in the underlying aligner
  cAlign.print(outputMode, 1.0, sout, screenWidth); 
}

void RefineSatsuma::processGap(const AlignmentBlock& lastGap,
              const AlignmentBlock& curr, ofstream& gapFout) const {
  if (!curr.isCompatible(lastGap)) { return; } 
  int tGapStart = lastGap.getTargetStop();
  int tGapStop  = curr.getTargetStart();
  int qGapStart = lastGap.getQueryStop();
  int qGapStop  = curr.getQueryStart();
  if (curr.isReversed()) {
      qGapStart = curr.getQueryStop();
      qGapStop  = lastGap.getQueryStart(); 
  }
  switch(outputMode) {
    case 0:  
    gapFout << "Aligning gap: " << lastGap.getTargetStop() 
      << " " << curr.getTargetStart() << ", " << qGapStart << " " << qGapStop << endl
      << lastGap.getTargetChrom() << " " << curr.getTargetChrom() 
      << " " << lastGap.getQueryChrom() << " " << curr.getQueryChrom() 
      << " " << lastGap.getOrient() << " " << curr.getOrient() 
      << tGapStop - tGapStart << " - "
      << qGapStop - qGapStart << endl;
      break;
    case 1:
    gapFout << lastGap.getTargetChrom() << "," << tGapStart << "," 
      << tGapStop << "," << tGapStop - tGapStart << "," 
      << curr.getQueryChrom() << "," << qGapStart << "," 
      << qGapStop << "," << qGapStop - qGapStart << ",";
  }

  // Check length of gap not to exceed certain limit
  if((tGapStop - tGapStart < MIN_GAP_SIZE || tGapStop - tGapStart > MAX_GAP_SIZE) 
      || (qGapStop - qGapStart < MIN_GAP_SIZE || qGapStop - qGapStart > MAX_GAP_SIZE)) {
    return; // Gap will not be processed
  }
  DNAVector q, t;
  t.SetToSubOf(target(lastGap.getTargetChrom()), tGapStart, tGapStop - tGapStart);
  t.SetName(lastGap.getTargetChrom());
  q.SetToSubOf(query(lastGap.getQueryChrom()), qGapStart, qGapStop - qGapStart);
  q.SetName(lastGap.getQueryChrom());
  if (lastGap.isReversed()) { q.ReverseComplement(); }

  gapFout << "##GAP" << endl;
  alignSeqs(t, tGapStart, '+', q, qGapStart, curr.getOrient(), gapFout);
}

void RefineSatsuma::printRealignInfo(const AlignmentBlock& aBlock, ofstream& realignFout) const {
  switch(outputMode) {
    case 0:
      realignFout << "Alignment (full): " << aBlock.getTargetChrom() << " : " 
                  << aBlock.getTargetStart() << " - " << aBlock.getTargetStop()-1
                  << " vs. " << aBlock.getQueryChrom() << " : " 
                  << aBlock.getQueryStart() << " - " << aBlock.getQueryStop()-1 
                  << " " << aBlock.getOrient() << endl
                  << aBlock.getTargetStop() - aBlock.getTargetStart() << " - "
                  << aBlock.getQueryStop() - aBlock.getQueryStart() << endl;
      break;
    case 1: 
    realignFout  << aBlock.getTargetChrom() << "," << aBlock.getTargetStart()
                 << "," << aBlock.getTargetStop() 
                 << "," << aBlock.getTargetStop() - aBlock.getTargetStart() 
                 << "," << aBlock.getQueryChrom()
                 << "," << aBlock.getQueryStart() << "," << aBlock.getQueryStop()
                 << "," << aBlock.getQueryStop() - aBlock.getQueryStart()
                 << "," << aBlock.getOrient() << ",";
  }
}
      
void RefineSatsuma::output(const string& str) const {
  cout<<str;
}

//======================================================
bool AlignBlocksThread::OnDo(const string & msg) {
  m_refinerUnit.alignSubset(m_fromIdx, m_toIdx, m_mergeBoundary, m_realignFout, m_gapFout, m_processLastBlock);
  return true;
}


//======================================================
