#ifndef _ALIGNMENT_COLA_H_
#define _ALIGNMENT_COLA_H_


#include "../Alignment.h"
#include "EditGraph.h"
#include "AlignerParams.h"

#define GAP_CHAR      '-'
#define MATCH_CHAR    '|'
#define MISMATCH_CHAR ' '


//===================================================================
/** 
 * Alignment class represents an alignment, which can be created from an aligned editGraph
 * The alignment can be represented by the one editgraphNode or by backtracing this and 
 * finding the aligned sequence coordinates.
 */
class AlignmentCola:public Alignment
{
public:

  // Ctor  - Can be used as default with no params or given as many of the params in order as available
  AlignmentCola(const DNAVector& tSeq = DNAVector(), const DNAVector& qSeq = DNAVector(),
                const AlignerParams& p = AlignerParams())
    :Alignment(tSeq,qSeq), pathNodes(), params(p)  {}

  virtual ~AlignmentCola() {} 

  /** Get the number of elements in the alignment (i.e. alignment length) */
  int getLength() const { return pathNodes.size(); }

  /** 
   * Add node to the optimal path nodes
   */
  void addNodeToPath(EditGraphNode* node);

  /**
   * Produce alignment from the path nodes
   * @parameter - choose whether beginning part of alignment that has negative score is traced
   */
  void traceAlignment(bool tracePrefix = false);

  /**
   * Confines the existing alignment to the start/end boundaries
   * Used for only keeping the alignment from start to end points
   * and discarding the rest. 
   * All fields of alignment/info will be updated accordingly
   */
  void keepSubalignment(int start, int end);
  
  /** 
   * Return the significance (P-Value) of the similarity or identity score:
   * This is done based on an empirical model of the expected score
   * as a function of the length of the alignment and assumption
   * of an underlying binomial distribution estimated by the findColaSigDist module
   */
  virtual double calcPVal() const;
  
  virtual void printFull(double pValLimit, ostream& sout,  int screenWidth, bool withInfo=true) const;

private:
  /** Use to update smith-waterman score after each iteration in finding the alignment path**/
  void updateSWScore();

  map<int,EditGraphNode> pathNodes;  /// Nodes on the optimal path mapped to increasing value: row+col 
  AlignerParams params;              /// Include the type of aligner and relevant params used for this alignment
};


#endif //_ALIGNMENT_COLA_H_

