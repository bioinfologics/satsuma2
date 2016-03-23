#ifndef _ALIGNMENTSET_H_
#define _ALIGNMENTSET_H_

#include "AlignmentCola.h"

#define MISMATCH_PENALTY -7.0 

//=====================================================================
/**
 * NOICOLA - None Overlapping Iterative - Contiguous Optimal Local Aligner:
 * The local alignment algorithm, Cola is ran on given regions to iteratively
 * find the best scoring none-overlapping alignments.
 */
class AlignmentSet 
{
public:
  /** 
   * Constructor just initialises the alignments 
   */
  AlignmentSet():alignments() {}
  
  /** 
   * Add an alignment to the list of alignments held by the set
   * A const reference to the alignment is accepted and used to create copy to keep
   */
  void add(const AlignmentCola& algn) { alignments.push_back(AlignmentCola(algn));}

  /** 
   * Return the Nth alignment in the order that they were created
   */
  const AlignmentCola& getNthAlignmentCola(int rank) const { return alignments[rank-1]; } 

  /** 
   * Use to print the alignment and a summary of the details 
   * @param outputType: choose what to print 0- full alignment 1-info in CSV format
   */
  void print(int outputType, double pValLimit, ostream& sout,  int screenWidth) const; 

private:
  vector<AlignmentCola> alignments;   /// Vector containing alignments in the order they were found
};



#endif //_NOICOLA_H_

