#ifndef _NSGAALIGNER_H_
#define _NSGAALIGNER_H_

#include "SWGAaligner.h"

class NSGAaligner: public SWGAaligner 
{
public:
  /** 
   * @param[in]  The target  sequence
   * @param[in]  The query sequence
   * Note that the targetSeq and querySeq are not copied 
   * N.B. The depth is set as a high number in the same fashion as NS
   */
  NSGAaligner(const DNAVector& tSeq, const DNAVector& qSeq, 
             const AlignerParams& p = AlignerParams(NSGA))
    :SWGAaligner(tSeq, qSeq, p, 10000) {} 

  ~NSGAaligner() {}

protected:
  /** 
   * Overrides the visitNode function in the parent class
   * @param[in]  The node to be visited 
   * @param[in] Index of column that should be checkpointed in current iteration
   */
  void visitNode(EditGraphNode* currNode, int currCheckpointColIndex); 
  
  /**
   * Used to visit nodes that are not handled by the vertical and horizontal moves
   * from the visitNode function inherited from SWGAligner
   * The nodes visited by this function include all the different diagonal move options
   */
  void visitNodeCola(EditGraphNode* currNode, int  currCheckpointColIndex); 
};


#endif //_NSGAALIGNER_H_
