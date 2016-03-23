#ifndef _IALIGNER_H_
#define _IALIGNER_H_

#include "AlignmentCola.h"

//=====================================================================
/**
 * Interface for Aligner objects (NSGAaligner, NSaligner, SWGAaligner, and SWaligner) 
 */
class IAligner 
{
public:
   virtual ~IAligner() {}
  /** 
   * main function to call for performing the alignment.
   * @return Returns the Alignment object which contains the alignment strings and other data
   */
  virtual const AlignmentCola& align() = 0;

  /** Returns the alignment object which contains the alignment strings and info **/
  virtual const AlignmentCola& getAlignment() = 0;

  /** Returns the target sequence used for the alignment */
  virtual const DNAVector& getTargetSeq() = 0; 

  /** Returns the query sequence used for the alignment */
  virtual const DNAVector& getQuerySeq() = 0; 
};



#endif //_IALIGNER_H_
