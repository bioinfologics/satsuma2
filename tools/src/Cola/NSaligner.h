#ifndef _NSALIGNER_H_
#define _NSALIGNER_H_


#include <string>
#include <cmath>
#include "EditGraph.h"
#include "AlignmentCola.h"
#include "IAligner.h"
#include "AlignerParams.h"

//=====================================================================
/**
 * Cola NS -  Contiguous Optimal Local Aligner - Nonlinear Scoring
 * This alignment algorithm is based on the basic Smith-Waterman
 * solution for pairwise alignment.
 * The main feature of this algorithm is:
 * Contiguous matches (i.e. matches without interruption by a gap or mismatch)
 * are tracked and scored accordingly. 
 */
class NSaligner: public IAligner 
{
public:
  /** 
   * @param[in]  The target  sequence
   * @param[in]  The query sequence
   * Note that the targetSeq and querySeq are not copied 
   * Default for penalties can be accessed from AlignerParams.h 
   */
  NSaligner(const DNAVector& tSeq, const DNAVector& qSeq, 
            const AlignerParams& p = AlignerParams(NS), int maxDepth=10000)
    :editGraph(tSeq.size(), qSeq.size(), maxDepth, p.getBandWidth()), alignment(tSeq, qSeq, p), params(p) {}

  ~NSaligner() {}

  /** 
   * main function to call for performing the alignment.
   * @return Returns the Alignment object which contains the alignment strings and other data
   */
  virtual const AlignmentCola& align();

  /** Returns the alignment object which contains the alignment strings and info **/
  virtual const AlignmentCola& getAlignment() { return alignment; }

  /** Returns the target sequence used for the alignment */
  virtual const DNAVector& getTargetSeq() { return alignment.getTargetSeq(); } 

  /** Returns the query sequence used for the alignment */
  virtual const DNAVector& getQuerySeq() { return alignment.getQuerySeq(); } 

protected:
  /** 
   * Visit all the nodes in the given subsection of the editgraph in turn
   * and assign their relevant score and find the middle node on 
   * the optimal path checkpoint column and pass it to next recursion.
   * @param[in] starting point (row and column), ending point (start and end)
   * @return Returns the mean contiguity depth traversed - Note that this would include affine depths
   */
  double traverseGraph(int startRow, int startCol, int endRow,
       int endCol, int endDepth, const EditGraphDepth& checkpointCell);

  /** 
   * Takes a node and decide how it should be visited to set its  score and origin
   * @param[in]  The node to be visited 
   * @param[in] Index of column that should be checkpointed in current iteration
   */
  virtual void visitNode(EditGraphNode* currNode, int currCheckpointColIndex); 

  /** 
   * Takes a node with contiguity depth of zero and sets the node scores
   * @param[in]  The node to be visited 
   * @param[in] Index of column that should be checkpointed in current iteration
   */
  void visitNodeContigZero(EditGraphNode* currNode, int currCheckpointColIndex); 

  /** 
   * Takes a node with contiguity depth of greater than zero and sets the node scores
   * @param[in]  The node to be visited 
   * @param[in] Index of column that should be checkpointed in current iteration
   */
  void visitNodeContigNoneZero(EditGraphNode* currNode, int currCheckpointColIndex);


  EditGraph     editGraph; /// The edit graph used for calculating the path scores
  AlignmentCola alignment; /// Object containing the backtraced alignment
  AlignerParams params;    /// The generic object which includes the relevant penalties
};



#endif //_NSALIGNER_H_

