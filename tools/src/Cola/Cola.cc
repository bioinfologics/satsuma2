#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "Cola.h"
#include "NSGAaligner.h"

const AlignmentCola& Cola::createAlignment(const DNAVector& tSeq, const DNAVector& qSeq,
              AlignerParams params) {
  IAligner* aligner;
  switch(params.getType()) {
    case NSGA: 
      aligner = new NSGAaligner(tSeq, qSeq, params);
      break;
    case NS: 
      aligner = new NSaligner(tSeq, qSeq, params);
      break;
    case SWGA: 
      aligner = new SWGAaligner(tSeq, qSeq, params);
      break;
    case SW: 
      aligner = new SWGAaligner(tSeq, qSeq, params); 
      break;
    default:
      //TODO error message
      aligner = new NSGAaligner(tSeq, qSeq, params);
  }
  latestAlignment = aligner->align();
  delete aligner;
  return latestAlignment;
}
