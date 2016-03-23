#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Cola/NOIAligner.h"

#define MIN_ALIGN_LEN 5 

//=====================================================================

const AlignmentSet& NOIAligner::align() {
  recurseAlign(0, targetSeq.isize(), 0, querySeq.isize());
  return alignments;
}

void NOIAligner::recurseAlign(int targetStart, int targetStop, int queryStart, int queryStop) {
  DNAVector t,q;
  t.SetToSubOf(targetSeq, targetStart, targetStop-targetStart);
  q.SetToSubOf(querySeq, queryStart, queryStop-queryStart);
  Cola cola1 = Cola();
  AlignmentCola algn = cola1.createAlignment(t, q, aligner); 
  // No significant alignment was found in region hence no recursion needed
  // TODO parameterise
  if(algn.calcPVal()>0.1) { return; } 

  alignments.add(algn); 

  // The part before the current alignment starts from target/queryStart upto target/queryStop1
  int targetStop1  = algn.getTargetOffset();
  int queryStop1   = algn.getQueryOffset();
  if((targetStop1-targetStart>MIN_ALIGN_LEN) && (queryStop1-queryStart>MIN_ALIGN_LEN)) {
    recurseAlign(targetStart, targetStop1, queryStart, queryStop1);
  }
  // The part after the current alignment starts from target/queryStart2 upto target/queryStop
  int targetStart2 = targetStart + algn.getTargetOffset() + algn.getTargetBaseAligned();
  int queryStart2  = queryStart  + algn.getQueryOffset()  + algn.getQueryBaseAligned();
  if((targetStop-targetStart2>MIN_ALIGN_LEN) && (queryStop-queryStart2>MIN_ALIGN_LEN)) {
    recurseAlign(targetStart2, targetStop, queryStart2, queryStop);
  }
}

