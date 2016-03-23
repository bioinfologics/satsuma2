#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_


#include <vector>
#include "../analysis/DNAVector.h"

#define GAP_CHAR      '-'
#define MATCH_CHAR    '|'
#define MISMATCH_CHAR ' '

//==================================================================
// Forward Declaration
class Alignment; 
class AlignmentCola;

//===================================================================
/** 
 * Alignment info contain all the different measures that have
 * been computed for a certain alignment.The alignment class 
 * contains an instance of this object and can be used for further src.
 */
class AlignmentInfo
{
  friend class Alignment;
  friend class AlignmentCola;

public:
  // Ctor  - Can be used as default with no params or given as many of the params in order as available
  AlignmentInfo(int tl=-1, int ql=-1, int to=-1, int qo=-1, 
                int tba=-1, int qba=-1, int tm=-1, int al=0, int sws= -1,
                double rs=-1, double rt=-1, double rtc=-1, double id=-1, double eval=-1)
    :tLen(tl), qLen(ql),
     tOffset(to), qOffset(qo), tBaseAligned(tba),
     qBaseAligned(qba), baseMatched(tm), alignmentLen(al), 
     smithWatermanScore(sws), rawScore(rs),
     runtime(rt), runtimeCoef(rtc), identity(id), eValue(eval) {}

  ~AlignmentInfo() {} 

  int getTargetLength() const       { return tLen;                            }
  int getQueryLength() const        { return qLen;                            }
  int getTargetOffset() const       { return tOffset;                         }
  int getQueryOffset() const        { return qOffset;                         }
  int getTargetStop() const         { return (tOffset+tBaseAligned-1);        }
  int getQueryStop() const          { return (qOffset+qBaseAligned-1);        }
  int getTargetBaseAligned() const  { return tBaseAligned;                    }
  int getQueryBaseAligned() const   { return qBaseAligned;                    }
  int getMaxBaseAligned() const     { return max(tBaseAligned, qBaseAligned); }
  int getMinBaseAligned() const     { return min(tBaseAligned, qBaseAligned); }
  int getAlignmentLen() const       { return alignmentLen;                    }
  int getBaseMatched() const        { return baseMatched;                     }
  int getSWScore() const            { return smithWatermanScore;              }
  double getRawScore()  const       { return rawScore;                        }
  double getRuntime() const         { return runtime;                         }
  double getRuntimeCoef()const      { return runtimeCoef;                     }
  double getEValue()const           { return eValue;                          }
  double getIdentity()              { return (identity==-1)?identity=calcIdentity():identity; }
  double getIdentity()const         { return (identity==-1)?calcIdentity():identity; }

  /** 
   * Use to obtain the identity score of the alignment.
   * The identitiy score is the ratio of the basepair matches 
   * over the entire length of the alignment
   */
  double calcIdentity() const;

  /** 
   * Use to reset the calculated alignment infos to zero.
   * This can also be used when mapping the alignment
   * reaches a null score and needs to set the offset values.
   */
  void resetCalculations(int tOffset, int qOffset);

protected:
  int tLen;               /// Length of the aligned target sequnce
  int qLen;               /// Length of the aligned query sequence
  int tOffset;            /// The position in the target sequence where the alignment starts 
  int qOffset;            /// The position in the query sequence where the alignment starts 
  int tBaseAligned;       /// The number of bases in the target that are part of the alignment
  int qBaseAligned;       /// The number of bases in the query that are part of the alignment
  int baseMatched;        /// The number of bases in the target/query that are matched
  int alignmentLen;       /// The total length of the alignment (including all match/mismach & insertions)
  int smithWatermanScore; /// The smithWaterman score of this alignment (+1:match, -1:mismatch and -1:gap)
  double rawScore;        /// The Raw score for the alignment (This score is native to the alignment method)
  double runtime;         /// Time take to find the alignment (in seconds)
  double runtimeCoef;     /// Factor by which runtime increases as opposed to basic SW
  double identity;        /// The number of matches over the entire length of the alignment
  double eValue;          /// The expected chance of obtaining such alignment by chance (For database search)
};

//===================================================================
/** 
 * Alignment class represents an alignment.
 * This class can also be used to hold a generic alignment regardless of method used for 
 * alignment (might need further extensions for certain alignment types).
 */
class Alignment
{
public:
  // Ctor  - Can be used as default with no params or given as many of the params in order as available
  Alignment(const DNAVector& tSeq = DNAVector(), const DNAVector& qSeq = DNAVector(), 
           const AlignmentInfo& inf = AlignmentInfo(),
            const string& tStr = string(), const string& qStr = string(), const string& mStr = string())
    :targetSeq(tSeq), querySeq(qSeq), info(inf), targetStr(tStr), 
     queryStr(qStr), matchesStr(mStr) ,
     targetIdxsInQuery(tSeq.size(), -1), queryIdxsInTarget(qSeq.size(), -1), 
     targetOrigOffset(0) , queryOrigOffset(0), targetSeqStrand('+'), querySeqStrand('+') {

    info.qLen = qSeq.size();
    info.tLen = tSeq.size();
  }

  virtual ~Alignment() {} 

  /** Get information available on the info object directly from the Alignment class */
  int getTargetLength() const       { return info.getTargetLength();      }
  int getQueryLength() const        { return info.getQueryLength();       }
  int getTargetOffset() const       { return info.getTargetOffset();      }
  int getQueryOffset() const        { return info.getQueryOffset();       }
  int getTargetBaseAligned() const  { return info.getTargetBaseAligned(); }
  int getQueryBaseAligned() const   { return info.getQueryBaseAligned();  }
  int getMaxBaseAligned() const     { return info.getMaxBaseAligned();    }
  int getMinBaseAligned() const     { return info.getMinBaseAligned();    }
  int getBaseMatched() const        { return info.getBaseMatched();       }
  int getAlignmentLen() const       { return info.getAlignmentLen();      }
  int getSWScore() const            { return info.getSWScore();           }
  double getRawScore() const        { return info.getRawScore();          }
  double calcIdentityScore() const  { return info.calcIdentity();         }
  double getIdentityScore()         { return info.getIdentity();          }
  double getIdentityScore() const   { return info.getIdentity();          }
  double getEValue()const           { return info.getEValue();            }
  double getRuntime() const         { return info.getRuntime();           }
  double getRuntimeCoef()const      { return info.getRuntimeCoef();       }

  /** Set information that is needed to be set external to the Alignment class */
  void setRuntime(double rt)        { info.runtime = rt; }
  void setRuntimeFactor(double rtc) { info.runtimeCoef = rtc; }

  /** Getter function for the auxillary information related to the query/target sequences */
  void getSeqAuxInfo(int& tOrigOffset, int& qOrigOffset, char& tSeqStrand, char& qSeqStrand); 

  /** Setter function for the auxillary information related to the query/target sequences */
  void setSeqAuxInfo(int targetOrigOffset, int queryOrigOffset, char targetSeqStrand, char querySeqStrand); 
  void setSeqAuxInfo(int tOrigOffset, int qOrigOffset, bool tSeqStrand, bool qSeqStrand); 

  /**
   * Prints either the full alignment or info in CSV format 
   * @param outputMode: choose what to print 0- full alignment 1-info in CSV format
   * Note that screenWidth is only used if choosing to print full alignment
   */
  virtual void print(int outputMode, double pValLimit, ostream& sout,  int screenWidth, bool withInfo=true) const;

  /**
   * Prints the alignment by assuming a minimum screen width
   * @param pValLimit is the upper bound for the  p-value of the alignment to be printed
   * @param sout the output stream for sending the output to
   * @param screenWidth is the width of each row of the alignment being printed
   */
  virtual string toString(int screenWidth, bool withInfo=true) const;
  virtual void printFull(double pValLimit, ostream& sout,  int screenWidth, bool withInfo=true) const;

  /**
   * Output the alignments in multiple alignment format
   * E.g. chr1 + pos-start ACT---AT--ATCG pos-end
   */
  void printMFAFormat(double pValLimit, ostream& sout, int screenWidth) const; 

  /**
   * Outputs results in NCBI Blast XML output format 
   */
  void printXMLFormat(double pValLimit, ostream& sout, int screenWidth) const; 

  /**
   * Output the complete set of information for an alignment in CSV format.
   * This includes the offsets, aligned bases, identity score and pValue
   */
  void printInfoCSV(ostream& sout) const; 


  /**
   * Use to get the corresponding pair in the query alignment for a given target base location
   * Note that the correspondence might be to a gap in which case GAP_CHAR is returned.
   */  
  virtual char getQueryAlignCharForTarget(int targetIndex) const; 

  /**
   * Use to get the corresponding index in the query alignment for a given target base location
   * Note that the correspondence might be to a gap in which case -1 is returned
   */  
  virtual int getQueryAlignIndexForTarget(int targetIndex) const; 

  /**
   * Use to get the corresponding pair in the target alignment for a given query base location
   * Note that the correspondence might be to a gap in which case GAP_CHAR is returned.
   */  
  virtual char getTargetAlignCharForQuery(int queryIndex) const; 

  /**
   * Use to get the corresponding index in the target alignment for a given query base location
   * Note that the correspondence might be to a gap in which case -1 is returned
   */  
  virtual int getTargetAlignIndexForQuery(int queryIndex) const; 

  /** Returns the target sequence used for the alignment */
  const DNAVector& getTargetSeq() const { return targetSeq; } 

  /** Returns the query sequence used for the alignment */
  const DNAVector& getQuerySeq() const { return querySeq; } 

  // Functions for obtaining stringified alignment output
  const string& getTargetString() const { return targetStr;  }
  const string& getQueryString()  const { return queryStr;   }
  const string& getMatchString()  const { return matchesStr; }

  /** Calculate and return the average length of ungapped regions in the alignment */
  double calcMeanContig() const;
  
  /** Calculate and return the edit count for the alignment */
  int calcEditCount() const;
  
  /** Calculate and return the modified Smith-Waterman score that limits gap penalties to 25 */  
  int calcModSWScore() const;

  /** 
   * Return the significance (P-Value) of the similarity or identity score:
   * This is done based on an empirical model of the expected score
   * as a function of the length of the alignment and assumption
   * of an underlying binomial distribution estimated by the findColaSigDist module
   */
  virtual double calcPVal() const;

  /**
   * Adapted from the method by which blast calculates it's raw score and E-value
   */
  virtual double calcEVal() const { return 0; } //TODO implement 

  /** 
   * Resetting containers before a new alignment is being set
   */
  void resetContainers();

protected:
  DNAVector targetSeq;               /// The target sequence that has been aligned to the query sequence
  DNAVector querySeq;                /// The query sequence that has been aligned to the target sequence
  AlignmentInfo info;                /// Object containing info such as score, offsets, etc.
  string targetStr;                  /// The aligned target substring
  string queryStr;                   /// The aligned query substring
  string matchesStr;                 /// String representation of alignment matches
  /// The following data is repetitive but kept as it's 
  /// relatively low cost and to save from recalculation on every use
  vector<int> targetIdxsInQuery;     /// The corresponding indexes of target sequence in the alignment
  vector<int> queryIdxsInTarget;     /// The corresponding indexes of query sequence in the alignment

  //Extra info about the target/query sequences
  int targetOrigOffset;              /// Original offset of the target sequence if it was extracted from a longer sequence
  int queryOrigOffset;               /// Original offset of the query sequence if it was extracted from a longer sequence
  char targetSeqStrand;              /// Target sequence strand
  char querySeqStrand;               /// Query sequence strand

};


#endif //_ALIGNMENT_H_

