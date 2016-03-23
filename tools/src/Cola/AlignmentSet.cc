#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Cola/AlignmentSet.h"

//=====================================================================

void AlignmentSet::print(int outputType, double pValLimit, ostream& sout,  int screenWidth) const {
  for( vector<AlignmentCola>::const_iterator iter = alignments.begin(); iter != alignments.end(); iter++) {
    (*iter).print(outputType, pValLimit, sout, screenWidth);
  }
}

