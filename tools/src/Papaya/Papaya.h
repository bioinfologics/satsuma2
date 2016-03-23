#ifndef PAPAYA_H
#define PAPAYA_H

#include <string>
#include "base/ThreadHandler.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerDynProg.h"
#include "analysis/SatsumaAlign.h"

//Forward Declaration
class PapayaAligner;

class PapayaAlignerParams{

  friend class PapayaAligner;

public:
  PapayaAlignerParams(): m_maxGap(200), m_minLen(16), m_minIdent(0.999) {}

  PapayaAlignerParams(int maxGap, int minLen, float minIdent
                     ): m_maxGap(maxGap), m_minLen(minLen), m_minIdent(minIdent) {}

private:
  int   m_maxGap;        /// maximum gap (def: 200)
  int   m_minLen;        /// Minimum alignment Length (def: 16)
  float m_minIdent;      /// Minimum identity (def: 0.999)
};

class PapayaSeedingParams{

  friend class PapayaAligner;

public:
  PapayaSeedingParams(): m_numKmers(2), m_lookAhead(12), m_lookAheadSkip(2), m_wordSize(12),
                         m_self(true), m_bNoReps(false) {}

  PapayaSeedingParams(int numKmers, int lookAhead, int lookAheadSkip, int wordSize, 
                      bool self, bool bNoReps
                     ): m_numKmers(numKmers), m_lookAhead(lookAhead), m_lookAheadSkip(lookAheadSkip), 
                        m_wordSize(wordSize), m_self(self), m_bNoReps(bNoReps)  {}

private:
  int   m_numKmers;      /// number of 12-mers used for seeding (def: 2)
  int   m_lookAhead;     /// number of lookahead k-mers (def: 12) 
  int   m_lookAheadSkip; /// skip k-mers in lookahea (def: 2)
  int   m_wordSize;      /// size of one k-mer (def:12)
  bool  m_self;          /// ignore matches (query=target  (def: false)
  bool  m_bNoReps;       /// do not use soft-masekd repeats (def: false) 
};

class PapayaQueryInstance
{
public:
  PapayaQueryInstance(const PapayaAlignerParams& params, 
                      const string& qFile, 
                      const string& oFile): m_alignParams(params), m_queryFile(qFile), m_outFile(oFile) {}

  PapayaQueryInstance(): m_alignParams(), m_queryFile(), m_outFile() {}
  
  PapayaAlignerParams getParams()    const  { return m_alignParams; }
  string              getQueryFile() const  { return m_queryFile;   }
  string              getOutFile()   const  { return m_outFile;     }

private:
  PapayaAlignerParams m_alignParams;   /// Parameters to use for the alignment
  string m_queryFile;                  /// Query sequence file to align against the reference
  string m_outFile;                    /// File where output is recorded
};


//Forward Declaration
class PapayaAligner;
class PapayaAlignerQueue 
{
public:
  PapayaAlignerQueue(): m_numOfAvailThreads(0), m_availableThreads(), 
                        m_queryQueue(), m_threadHandler(), m_lock(), m_alignerUnit(NULL) {}
  
  PapayaAlignerQueue(PapayaAligner* aligner, int threadCnt=1): m_numOfAvailThreads(threadCnt), m_availableThreads(), 
                                     m_queryQueue(), m_threadHandler(), m_lock(), m_alignerUnit(NULL) {
    init(threadCnt, aligner);
  }

  void init(int numOfThreads, PapayaAligner* aligner) {
    m_alignerUnit = aligner;
    m_availableThreads.resize(numOfThreads);
    for(int i=0; i<numOfThreads; i++){ 
      m_availableThreads[i] = i; 
    }
    m_numOfAvailThreads = numOfThreads;
  }

  void receiveQuery(const PapayaAlignerParams params,
                    const string& queryFile, const string& outFile);
  void releaseThread(int threadQueueId);
  void checkQueue(); 

private:
  void processThread(int threadQueueId); 

  int m_numOfAvailThreads;                /// Number of available threads
  svec<int> m_availableThreads;
  svec<PapayaQueryInstance> m_queryQueue;
  ThreadHandler m_threadHandler;          /// Unit to handle execution of threads
  pthread_mutex_t m_lock;
  PapayaAligner*  m_alignerUnit;          /// Aligner unit which handles the alignment  (not this is not a const as threads need to be released)
};

//Forward Declaration
class PapayaAlignThread;

class PapayaAligner{
  friend class PapayaAlignThread;
public:
  PapayaAligner(): m_target(), m_KSAligner(), m_seedingParams(), m_queryQueue() {}

  void initialize(const string& targetFile, const PapayaSeedingParams& params, int numOfThreads);
  void processQuery(const string& queryFile, const string& outFile, const PapayaAlignerParams& alignParams); 
  void releaseThread(int threadQueueId);

protected:
  void align(const string& queryFile, const string& outFile, const PapayaAlignerParams& alignParams);

private:
  void align(const vecDNAVector& query, const string& outFile, const PapayaAlignerParams& alignParams); 

  vecDNAVector          m_target;         /// Reference used for initialization and seeding
  KmerSuperAligner      m_KSAligner;      /// Kmer super aligner instance initilized with the target reference
  PapayaSeedingParams   m_seedingParams;  /// Alignment-seeding parameters
  PapayaAlignerQueue    m_queryQueue;     /// Queue to handle the incoming alignment queries
};

class PapayaAlignThread : public IOneThread
{
public:
  PapayaAlignThread(PapayaAligner* alignerUnit, 
                    const PapayaAlignerParams params, const string& queryFile,
                    const string& outFile, int threadId): m_alignerUnit(alignerUnit), m_alignParams(params),
                                                          m_queryFile(queryFile), m_outFile(outFile), m_threadQueueId(threadId) {}

  PapayaAlignThread(const PapayaQueryInstance& pInst, PapayaAligner* alignerUnit, int threadId
                   ): m_alignerUnit(alignerUnit), m_alignParams(pInst.getParams()),
                      m_queryFile(pInst.getQueryFile()), 
                      m_outFile(pInst.getOutFile()), m_threadQueueId(threadId) {}

  void setThreadId(int threadId) { m_threadQueueId = threadId; }

protected:
  virtual bool OnDie() { return true; }
  virtual bool OnDo(const string & msg); 
  virtual bool OnInitialize(const string & msg) { return true; }

private:
  PapayaAligner*  m_alignerUnit;       /// Aligner unit which handles the alignment  (not this is not a const as threads need to be released)
  PapayaAlignerParams m_alignParams;   /// Parameters to use for the alignment
  string m_queryFile;                  /// Query sequence file to align against the reference
  string m_outFile;                    /// File where output is recorded
  int m_threadQueueId;                /// Id for this thread
};


#endif //PAPAYA_H 
