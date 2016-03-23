



#ifndef _KMERDYNPROG_H_
#define _KMERDYNPROG_H_


#include "aligns/KmerAlignCore.h"


#define KSA_K 24
#define NULL_KMER_ID -1

// This is how the alignment comes back...
class SuperMatch
{
 public:
  SuperMatch() {
    m_contig = -1;
    m_firstBase = -1;
    m_lastBase = -1;
    m_refID = -1;
    m_rc = false;
    m_refStart = -1;
    m_refEnd = -1;
    m_score = 0;
  }

  void Set(int contig,
	   int firstBase,
	   int lastBase,
	   int refID,
	   bool rc,
	   int refStart,
	   int refEnd)
    {
      m_contig = contig;
      m_firstBase = firstBase;
      m_lastBase = lastBase;
      m_refID = refID;
      m_rc = rc;
      m_refStart = refStart;
      m_refEnd = refEnd;
    }

  void PrettyPrint(int realContigID = -1) const {
    cout << "c";
    if (realContigID != -1) {
      cout << realContigID;
    } else {
      cout << m_contig;
    }
    cout << " [ " << m_firstBase << " - " << m_lastBase << " ] ==> ";
    cout << m_refID << " [ " << m_refStart << " - " << m_refEnd << " ] ";
    if (m_rc) {
      cout << "rc";
    } else {
      cout << "fw";
    }
    cout << " len=" << m_lastBase - m_firstBase + 1 << endl;
  }

  void SetScore(int s) {m_score = s;}

  int GetContig() const {return  m_contig;}
  int GetFirstBase() const {return m_firstBase;}
  int GetLastBase() const {return m_lastBase;}
  int GetRefID() const {return m_refID;}
  bool GetRC() const {return m_rc;}
  int GetRefStart() const {return m_refStart;}
  int GetRefEnd() const {return m_refEnd;}
  int GetScore()  const {return m_score;}


  bool IsContiguous(const SuperMatch & two) const {
    if (m_contig != two.m_contig)
      return false;
    if (m_refID != two.m_refID)
      return false;
    if (m_rc != two.m_rc)
      return false;
    if (!m_rc) {
      if (m_refStart + 1 != two.m_refStart)
	return false;
    } else {
      if (m_refStart - 1 != two.m_refStart)
	return false;   
    }
    return true;
  }

 private:
  int m_contig;
  int m_firstBase;
  int m_lastBase;
  int m_refID;
  bool m_rc;
  int m_refStart;
  int m_refEnd;
  int m_score;
};



class SuperAlign
{
 public:
  SuperAlign() {}

  void AddMatch(const SuperMatch & m) {m_matches.push_back(m);}

  int GetMatchCount() const {return m_matches.isize();}
  const SuperMatch & GetMatch(int i) const {return m_matches[i];}
  const svec<SuperMatch> & GetMatches() const {return m_matches;}
private:
  svec<SuperMatch> m_matches;
};



// Dynamic programming stuff...


class KSAMatch
{
 public:
  KSAMatch() {
    m_refID = -1;
    m_refStart = 0;
    m_prevIndex = -1;
    //m_rc = false;
    m_score = 0x7FFFFFFF;
    m_feed = -1;
    m_penalty = 0;
  }

  void Set(int refID, 
	   int refStart, 
	   int refStop,
	   bool rc,
	   int penalty = 0) {
    //if (refStart < 0)
    //cout << "refStart=" << refStart << endl;
    if (penalty > 0)
      m_penalty = (m_penalty | 0x80000000);
    m_refStart = (refStart & 0x7FFFFFFF);
    m_refID = refID;
    if (rc)
       m_refStart = (m_refStart | 0x80000000);
  }

  int GetRefID() const {return m_refID;}
  int GetRefStart() const {return (m_refStart & 0x7FFFFFFF);}
  int GetRefStop() const {return GetRefStart() + KSA_K - 1;}
  bool IsRC() const {
    if ((m_refStart & 0x80000000) == 0)
      return false;
    else
      return true;
  }

  bool IsContiguous(const KSAMatch & m, int shift = 1) const {
    if (GetRefID() == NULL_KMER_ID &&  m.GetRefID() == NULL_KMER_ID)
      return true;

    if (GetRefID() != m.GetRefID())
      return false;
    if (IsRC() != m.IsRC())
      return false;
    if (IsRC()) {
      if (GetRefStart() - shift != m.GetRefStart())      
	return false;
    } else {
      if (GetRefStart() + shift != m.GetRefStart())
	return false;
    }
    return true;
  }

  void SetPrevIndex(int i) {m_prevIndex = i;}
  void SetPrevFrame(int i) {
    m_penalty = ((m_penalty & 0x80000000) + i);
  }

  void SetScore(int s) {m_score = s;}
  int GetPrevIndex() const {return m_prevIndex;}
  int GetPrevFrame() const {return (m_penalty &  0x7FFFFFFF);}
  int GetScore() const {
    //if (m_penalty > 0)
    //cout << "Adding penalty " << m_penalty << endl;
    if (m_score < 0x7FFFFFFF)
      return m_score + GetPenalty();
    else 
      return m_score;
  }
  int GetPenalty() const {
    if ((m_penalty & 0x80000000) == 0)
      return 0;
    else
      return 1;
  }
  

  int GetFeed() const {return m_feed;}
  void SetFeed(int f) {m_feed = f;}

  bool operator < (const KSAMatch & m) const {
    if (IsRC() != m.IsRC()) {
      if (IsRC())
	return false;
      else
	return true;
    }
    if (m_refID != m.m_refID)
      return (m_refID < m.m_refID);
    return (GetRefStart() < m.GetRefStart());

  }


 private:
  int m_refID;
  int m_feed;

  int m_refStart;
  int m_prevIndex;
  int m_score;

  int m_penalty;
  //bool m_rc;
};


//=======================================================================================



class KSAItem
{
 public:
  KSAItem() {
    m_fromOnSuper = -1;
    m_fromOnContig = -1;
    m_contigIndex = -1;
    m_dev = 0;
  }

  void Set(int fromOnSuper,
	   int fromOnContig,
	   int contigIndex,
	   int dev)
  {
    m_fromOnSuper = fromOnSuper;
    m_fromOnContig = fromOnContig;
    m_contigIndex = contigIndex;
    m_dev = dev;
  }


  int GetFromOnSuper() const {return m_fromOnSuper;}
  int GetFromOnContig() const {return m_fromOnContig;}
  int GetContigIndex() const {return m_contigIndex;}
  int GetCount() const {return m_matches.isize();}
  
  int GetDev() const {return m_dev;}

  const KSAMatch & Get(int i) const {return m_matches[i];}
  KSAMatch & GetIt(int i) {return m_matches[i];}

  bool Update(int i, int score, int prevIndex, int prevFrame) {
    if (score <= m_matches[i].GetScore()) {
      m_matches[i].SetPrevIndex(prevIndex);
      m_matches[i].SetPrevFrame(prevFrame);
      m_matches[i].SetScore(score);
      return true;
    }
    return false;
  }

  void AddMatch(const KSAMatch & m) {
    m_matches.push_back(m);
  }

  void AddMatch(int refID, 
		int refStart, 
		int refStop,
		bool rc) {
    KSAMatch m;
    m.Set(refID, refStart, refStop, rc);
    m_matches.push_back(m);
  }

  int GetSize() const {return m_matches.isize();}
  void SetSize(int i) {
    m_matches.resize(i);
  }

  void AddMatch(int refID, 
		int refStart, 
		int refStop,
		bool rc,
		int index,
		int penalty) {
    //if (penalty > 0)
    //cout << "Set penalty " << penalty << endl;
    m_matches[index].Set(refID, refStart, refStop, rc, penalty);
  }


  void SortMatches() {
    Sort(m_matches);
  }

 private:
  int m_fromOnSuper;
  int m_fromOnContig;
  int m_contigIndex;
  int m_dev;
  svec<KSAMatch> m_matches;
};


    



class MatchScore
{
 public:
  MatchScore( ) {
    m_nullScore = 1;
    m_globScore = 50000;
  }

  void SetLog( ostream *log ) { m_log = log; }
  int GetScore(const KSAMatch & from, const KSAMatch & to, int sep, int dev);
  int GetGlobalDistanceLimit() const {return m_globScore;}

 private:
  ostream *m_log;
  int m_nullScore;
  int m_globScore;

};


// Finally, the aligner itself...
class KmerSuperAligner
{
 public:
  KmerSuperAligner( ostream *m_log = 0 );
  ~KmerSuperAligner();


  void SetRefBases(const vecDNAVector & b);

  // Allow for repeat masking
  void SetRefBases(const vecDNAVector & b, const vecNumVector & tags, int min = 4);

  void Align(SuperAlign & result, 
	     const vecDNAVector & contigBases,
	     svec<int> & contigStartInSuper,
	     svec<int> & gapDevInSuper,
	     int skipContig = -1);

  void SetKmerShift(int s) {
    m_kmerShift = s;
  }

  void SetMaxFreq(int i) {m_maxFreq = i;}

  void SetUseCombmers(bool b) {
    m_bUseCombs = b;
  }
  void SetUseProteins(bool b) {
    m_bUseProteins = b;
  }

  void SetNumKmers(int i) {
    m_numKmers = i;
  }

  void SetNewLookahead(int frames, int skip) {
    m_lookAheadFrames = frames;
    m_lookAheadStep = skip;
  }

  void SetLookAhead(int i) {m_core.SetLookAhead(i);}
  void SetLookAheadMaxFreq(int i) {
    m_core.SetLAMaxFreq(i);
    m_lookAheadMaxFreq = i;
  }


  void SetWordSize(int i) {
    m_numberTrans.SetSize(i);
    m_combTrans.SetSize(i);
    m_proteinTrans.SetSize(i);
   }

 private:
  // Returns the index to the best score...
  int MatchCloseExtensions(int fromIndex, int toIndex, int addScore);
  int ExtendAll(int fromIndex, int bestIndex, int toIndex);

  KmerAlignCore<KmerAlignCoreRecordWithScore> m_core;
  TranslateBasesToNumberExact m_numberTrans;
  TranslateBasesToNumberDoubleComb m_combTrans;
  TranslateBasesToNumberProtein m_proteinTrans;
  int m_numKmers;
  bool m_bUseCombs;
  bool m_bUseProteins;

  bool m_bUseCore;

  svec<KSAItem> m_dynArray;
  MatchScore m_matchScore;
  int m_lastBestScore;
  int m_upCount;

  int m_kmerShift;

  int m_maxFreq;

  ostream *m_log;

  int m_lookAheadFrames;
  int m_lookAheadStep;
  int m_lookAheadMaxFreq;


};


#endif //_KMERSUPERALIGNER_H_

