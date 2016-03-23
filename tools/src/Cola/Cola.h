#ifndef _COLA_H_
#define _COLA_H_

#include "AlignmentCola.h"
#include "AlignerParams.h"

//=====================================================================

/**
 * Factory class used for obtaining one of the aligner types:
 * NSGAaligner, NSaligner, SWGAaligner, and SWaligner
 */
class Cola
{
public:
  /** 
   * @return pointer to the requested aligner object 
   */
  const AlignmentCola& createAlignment(const DNAVector& tSeq, const DNAVector& qSeq,
              AlignerParams params); 

  const AlignmentCola getAlignment() { return latestAlignment; }

private:
  AlignmentCola latestAlignment;
};

#endif //_COLA_H_
