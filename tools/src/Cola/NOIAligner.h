#ifndef _NOICOLA_H_
#define _NOICOLA_H_

#include "Cola.h"
#include "AlignmentSet.h"

//=====================================================================
/**
 * NOICOLA - None Overlapping Iterative - Contiguous Optimal Local Aligner:
 * The local alignment algorithm, Cola is ran on given regions to iteratively
 * find the best scoring none-overlapping alignments.
 */
class NOIAligner 
{
public:
  /** 
   * @param[in]  The target  sequence
   * @param[in]  The query sequence
   * Note that the targetSeq and querySeq are not copied 
   */
  NOIAligner(const DNAVector& tSeq, const DNAVector& qSeq)
    :targetSeq(tSeq), querySeq(qSeq), alignments(), aligner() {}

   ~NOIAligner() {}

  /** 
   *  Set the Type of aligner to use and also the relevant parameters 
   *  If this function is not called the class will use its defaults,
   *  which is the NSGA aligner and its default parameters.
   */
  void setAligner(AlignerParams algn) { aligner = algn; }

  /** 
   * main function to call for performing the alignment.
   * @return Returns the Alignment object which contains the alignment strings and other data
   */
  const AlignmentSet& align();

  /** Returns the alignment object which contains the alignments strings and info **/
  AlignmentSet* getAlignments() { return &alignments; }

  /** Returns the target sequence used for the alignment */
  const DNAVector& getTargetSeq() { return targetSeq; } 

  /** Returns the query sequence used for the alignment */
  const DNAVector& getQuerySeq() { return querySeq; } 

private:
  DNAVector targetSeq;        /// The target sequence that has been aligned to the query sequence
  DNAVector querySeq;         /// The query sequence that has been aligned to the target sequence
  AlignmentSet alignments;    /// Object containing the noneoverlapping alignments
  AlignerParams aligner;      /// Object containing aligner type and relevant parameters

  void recurseAlign(int targetStart, int targetStop, int queryStart, int queryStop); 
};



#endif //_NOICOLA_H_

