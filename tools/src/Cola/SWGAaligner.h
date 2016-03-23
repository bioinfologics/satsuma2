#ifndef _SWGAALIGNER_H_
#define _SWGAALIGNER_H_

#include "NSaligner.h"

class SWGAaligner: public NSaligner 
{
public:
  /** 
   * @param[in]  The target sequence
   * @param[in]  The query sequence
   * Note that the targetSeq and querySeq are not copied 
   * N.B. The depth is set as two here to use three elements for each cell
   * one for the vertical, once for horizontal, and one for the diagonal move
   */
  SWGAaligner(const DNAVector& tSeq, const DNAVector& qSeq, 
              const AlignerParams& p = AlignerParams(SWGA), int maxD = 2) 
    :NSaligner(tSeq, qSeq, p, maxD) {}

  ~SWGAaligner() {}

protected:
  /** 
   * Overrides the visitNode function in the parent class (i.e. NSaligner)
   * @param[in]  The node to be visited 
   * @param[in] Index of column that should be checkpointed in current iteration
   */
  virtual  void visitNode(EditGraphNode* currNode, int currCheckpointColIndex); 
};


#endif //_SWGAALIGNER_H_
