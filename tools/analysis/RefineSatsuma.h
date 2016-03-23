#ifndef _REFINESATSUMA_H_
#define _REFINESATSUMA_H_

#include "../base/ThreadHandler.h"
#include "DNAVector.h"
#include "../src/AlignmentBlock.h"
#include "../src/Cola/Cola.h"

// Maximum and minimum gap size allowed for alignment 
#define MAX_GAP_SIZE 15000 
#define MIN_GAP_SIZE 1 

//==================================================================
/**
 * This class facilitates taking satsuma alignment coordinates and 
 * the sequences that were used in the satsuma alignment and to postprocess
 * these with a more optimal alignment tool such as Cola
 * The class is templated over the refining alignment class which is assumed
 * to have an "align" function associated with it. 
 **/
class RefineSatsuma 
{
public:
  RefineSatsuma(const vecDNAVector& t, const vecDNAVector& q, 
                 const string& satsumaFile, 
                 const string& reFile, const string& gFile, int oMode
               ):target(t), query(q), realignFile(reFile),
           gapFile(gFile), outputMode(oMode), screenWidth(200),
           aligner()
  {
    init(satsumaFile);
  }
  
  /** 
   *  Set the Type of aligner to use and also the relevant parameters 
   *  If this function is not called the class will use its defaults,
   *  which is the NSGA aligner and its default parameters.
   */
  void setAligner(AlignerParams algn) { aligner = algn; }

  /** Realign the aligned regions and also attempt to align gaps 
   *  @param outputMode: 
   */
  void alignAll(int mergeBoundary, int numOfThreads) const;
  void alignSubset(int startIdx, int endIdx, int mergeBound, ofstream& realignFout, 
                   ofstream& gapFout, bool processLastBlock) const; 

private:
  /**
   * load blocks from a given file for initialisation 
   */
   void init(const string& satsumaFile);
  /** 
   * Maps the start and end coordinates of an alignment block on to the 
   * orignating query/target sequences and realigns those regions
   */ 
  void alignBlock(const AlignmentBlock& aBlock, ofstream& realignFout) const; 

  /**
   * Processes a gap by aligning it if it is between the min-max limits
   */
  void processGap(const AlignmentBlock& lastGap, const AlignmentBlock& curr, ofstream& gapFout) const; 

  /** 
   * Create an instance of the alignerType class and align the given sequences 
   */
  void alignSeqs(const DNAVector& t, int tOffset, char tStrand, 
                 const DNAVector& q, int qOffset, char qStrand, ostream& sout=cout) const; 
  
  void printRealignInfo(const AlignmentBlock& aBlock, ofstream& realignFout) const; 
  void printGapInfo(const AlignmentBlock& lastGap, const AlignmentBlock& curr) const; 

  /** handle output of alignments and other info */
  // TODO handle output instead of scattered couts
  void output(const string& str) const;

  const vecDNAVector& target;           /// target genome contained in chromosomes
  const vecDNAVector& query;            /// query genome contained in chromosomes
  vector<AlignmentBlock> blockElements; /// vector containing the input blocks from the satsuma file
  string realignFile;                   /// The output file for the realignments
  string gapFile;                       /// The output file for the gap alignments
  int outputMode;                       /// 0-full alignments and info 1- only info needed for extracting stats 
  int screenWidth;                      /// Width to print alignments horizontally
  AlignerParams aligner;                /// Object containg aligner type and necessary parameters
};

class AlignBlocksThread : public IOneThread
{
public:
  AlignBlocksThread(const RefineSatsuma& refUnit,
                    int mBound, const string& realignFile, const string& gapFile, 
                    int from, int to, bool processLB, int tn): m_refinerUnit(refUnit), m_mergeBoundary(mBound),
                                                               m_fromIdx(from), m_toIdx(to), 
                                                               m_processLastBlock(processLB), m_threadIdx(tn) {
    m_realignFout.open(realignFile.c_str(), std::ofstream::app);
    m_gapFout.open(gapFile.c_str(), std::ofstream::app);
  }

  ~AlignBlocksThread() { 
    m_realignFout.close();
    m_gapFout.close();
  }

protected:
  virtual bool OnDie() { return true; }
  virtual bool OnDo(const string & msg); 
  virtual bool OnInitialize(const string & msg) {    return true; }

private:
  const RefineSatsuma&  m_refinerUnit;  
  int m_mergeBoundary;
  ofstream m_realignFout;                /// The outputstream for the realignments
  ofstream m_gapFout;                    /// The outputstream for the gap alignments
  int m_fromIdx;
  int m_toIdx;
  bool m_processLastBlock;               /// Flag indicating whether this is the process running the last block
  int m_threadIdx;
};


#endif //_REFINESATSUMA_H_
