#include <string>

#include "analysis/DNAVector.h"
#include "aligns/KmerDynProg.h"
#include "src/DPAlign.h"
#include "analysis/SatsumaAlign.h"
#include "src/Cola/Cola.h"
#include "src/Papaya/Papaya.h"

 
void PapayaAligner::initialize(const string& targetFile, const PapayaSeedingParams& params, int numOfThreads) {
  m_seedingParams = params;
  svec<string> queryNames, targetNames;
  if (m_seedingParams.m_bNoReps) {
    cout << "Reading target (masked)..." << endl;
    m_target.Read(targetFile, false, false, false);
  } else {
    cout << "Reading target..." << endl;
    m_target.Read(targetFile, targetNames);
  }

  m_KSAligner.SetWordSize(m_seedingParams.m_wordSize);
  m_KSAligner.SetNumKmers(m_seedingParams.m_numKmers);
  m_KSAligner.SetLookAhead(m_seedingParams.m_lookAhead);
  m_KSAligner.SetNewLookahead(m_seedingParams.m_lookAhead, m_seedingParams.m_lookAheadSkip);
  m_KSAligner.SetRefBases(m_target);

  m_queryQueue.init(numOfThreads, this);

}

void PapayaAligner::processQuery(const string& queryFile, const string& outFile, const PapayaAlignerParams& alignParams) {
  m_queryQueue.receiveQuery(alignParams, queryFile, outFile);
} 
 
void PapayaAligner::align(const string& queryFile, const string& outFile, const PapayaAlignerParams& alignParams) {
  svec<string> queryNames;
  vecDNAVector query;
  if (m_seedingParams.m_bNoReps) {
    cout << "Reading query (masked)..." << endl;
    query.Read(queryFile, false, false, false);
  } else {
    cout << "Reading query..." << endl;
    query.Read(queryFile, queryNames);
  }
  align(query, outFile, alignParams);
}

void PapayaAligner::align(const vecDNAVector& query, const string& outFile, const PapayaAlignerParams& alignParams) {
  ofstream outStrm(outFile.c_str(), ios::out | ios::binary);
  for(int i=0; i<query.isize(); i++) {
    vecDNAVector tmpBases;
    tmpBases.push_back(query[i]);
    svec<int> contigStarts;
    svec<int> contigDevs;
    contigStarts.push_back(0);
    contigDevs.push_back(0);

    SuperAlign result;
    int skip = -1;
    if (m_seedingParams.m_self) { skip = i; }
    m_KSAligner.Align(result, tmpBases, contigStarts, contigDevs, skip);
    for (int j=0; j<result.GetMatchCount(); j++) {
      const SuperMatch & one = result.GetMatch(j);
      int k;
      for (k=j+1; j<result.GetMatchCount(); k++) {
        const SuperMatch & two= result.GetMatch(k);
        if (two.GetRefID() != one.GetRefID()) { 
          break;
        }
        if (two.GetContig() != one.GetContig()) {
          break;
        }
        if (two.GetRC() != one.GetRC()) {
          break;
        }
        if (two.GetFirstBase() > result.GetMatch(k-1).GetLastBase() + alignParams.m_maxGap) {
         break;
        }
        if (one.GetRC()==false) {
          if (two.GetRefStart() > result.GetMatch(k-1).GetRefEnd() + alignParams.m_maxGap) {
            break;
          }
          if (two.GetRefStart() < result.GetMatch(k-1).GetRefStart()) {
            break;
          }
        } else {
          if (two.GetRefEnd() < result.GetMatch(k-1).GetRefStart() - alignParams.m_maxGap) {
            break;
          }
          if (two.GetRefStart() > result.GetMatch(k-1).GetRefStart()) {
            break;
          }
        }
      }
      const SuperMatch & two= result.GetMatch(k-1);

      // Align
      int startTarget = one.GetRefStart();
      int endTarget = two.GetRefEnd();

      int startQuery = one.GetFirstBase();
      int endQuery = two.GetLastBase();

      if (one.GetRC() == false) {
      } else {
        startTarget = two.GetRefStart();
        endTarget = one.GetRefEnd();
      }

      j = k-1; //Update index by setting upto where merge blocks have been merged, ready for next iteration

      if (endQuery - startQuery < alignParams.m_minLen)
        continue;
  
      DNAVector t, q;
   
      string tName     = ">" + string(m_target.NameClean(one.GetRefID()));
      //string qName     = ">" + string(query.NameClean(one.GetContig()));
      string qName     = ">" + string(query.NameClean(i));
      const DNAVector & t_full = m_target(tName);
      const DNAVector & q_full = query(qName);

      t.SetToSubOf(t_full, startTarget, endTarget-startTarget+1);
      t.SetName(tName);
      q.SetToSubOf(q_full, startQuery, endQuery-startQuery+1);
      q.SetName(qName);

      if (one.GetRC()) {
        q.ReverseComplement();
      }

      Cola cola1 = Cola();
      AlignmentCola cAlign = cola1.createAlignment(t, q, AlignerParams());
      cAlign.setSeqAuxInfo(startTarget, startQuery, true, !one.GetRC()); 
      cAlign.print(2, alignParams.m_minIdent, outStrm, 150); //TODO Parameterise 
    }
  }

  //outStrm << endl << "<DONE/>" << endl;
  outStrm.close();
  cout << "Done aligning - Closing. " << outFile << endl;
  string doneName = outFile + ".done";
  FILE * pDone = fopen(doneName.c_str(), "w");
  fprintf(pDone, "done\n");
  fclose(pDone);
  cout << "All done!" << endl;
}


void PapayaAligner::releaseThread(int threadQueueId){
  m_queryQueue.releaseThread(threadQueueId);
}

//======================================================
bool PapayaAlignThread::OnDo(const string & msg) {
  cout << "Processing: " << msg << endl;
  m_alignerUnit->align(m_queryFile, m_outFile, m_alignParams);    
  m_alignerUnit->releaseThread(m_threadQueueId);
  return true;
}
//======================================================



//======================================================
void PapayaAlignerQueue::receiveQuery(const PapayaAlignerParams params,
                                      const string& queryFile, const string& outFile) {
  pthread_mutex_lock(&m_lock);
  m_queryQueue.push_back(PapayaQueryInstance(params, queryFile, outFile));
  pthread_mutex_unlock(&m_lock);
  checkQueue(); 
}

void PapayaAlignerQueue::releaseThread(int threadQueueId){
  pthread_mutex_lock(&m_lock);
  m_availableThreads[m_numOfAvailThreads] = threadQueueId;
  m_numOfAvailThreads++;
  pthread_mutex_unlock(&m_lock);
  checkQueue();
}

void PapayaAlignerQueue::checkQueue() {
  pthread_mutex_lock(&m_lock);
  while(m_numOfAvailThreads>0) {
    if(m_queryQueue.isize()>0) {
      processThread(m_availableThreads[m_numOfAvailThreads-1]);
      m_numOfAvailThreads--;
    } else {
      pthread_mutex_unlock(&m_lock);
      return; 
    }
  }
  pthread_mutex_unlock(&m_lock);
}
 
void PapayaAlignerQueue::processThread(int threadQueueId) {
   char tmp[256];
   sprintf(tmp, "Processing query file: %s Writing to outputFile: %s Processed by thread%d",
           m_queryQueue[0].getQueryFile().c_str(), m_queryQueue[0].getOutFile().c_str(), threadQueueId);
   string init = "init_";
   init += tmp;
   int handlerThreadId = m_threadHandler.AddThread(new PapayaAlignThread(m_queryQueue[0], m_alignerUnit, threadQueueId));
   m_queryQueue.erase(m_queryQueue.begin()); 
   m_threadHandler.Feed(handlerThreadId, init);
}
//======================================================
