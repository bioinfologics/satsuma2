#ifndef SEQCHUNK_H_
#define SEQCHUNK_H_

//#include "base/SVector.h"
#include "analysis/DNAVector.h"



//================================================================

class SeqChunkSelect
{
 public:
  SeqChunkSelect() {}

  void Read(const string & file);
  void Write(const string & file);
  void MergeRead(const string & file);


  int GetQuerySize() const {return m_query.isize();}
  int GetTargetSize() const {return m_target.isize();}
  int GetQuery(int i) const {return m_query[i];}
  int GetTarget(int i) const {return m_target[i];}

  const svec<int> & Query() const {return m_query;}
  const svec<int> & Target() const {return m_target;}

  void Set(int querySize, int targetSize) {
    m_query.resize(querySize, 0);
    m_target.resize(targetSize, 0);
  }
  void SetQuery(int i, int v) {m_query[i] = v;}
  void SetTarget(int i, int v) {m_target[i] = v;}

 private:
  svec<int> m_query;
  svec<int> m_target;
};




class SeqChunk
{
public:
  SeqChunk() {
    m_id = -1;
    m_ori = 0;
  }

  void Set(const string & name, 
	   int start,
	   int id) {
    m_name = name;
    m_start = start;
    m_id = id;
  }

  const string & GetName() const {return m_name;}
  int GetStart() const {return m_start;}
  int GetID() const {return m_id;}

  void ForceOrientation(int i) {m_ori = i;}
  int GetForceOrientation() const {return m_ori;}

  void Print() {
    cout << "Sequence " << m_name << " offset " << m_start;
  }

private:
  string m_name;
  int m_start;
  int m_id;
  int m_ori;
};

class ChunkManager
{
public:
  ChunkManager(int size, int overlap) {   
    m_size = size;
    m_overlap = overlap;
  }

  int GetCount() const {return m_lengths.isize();}
  int GetSize(int i) const {return m_lengths[i];}
  
  void ChunkIt(vecDNAVector & out, 
	       svec<SeqChunk> & chunks, 
	       const vecDNAVector & in, 
	       const svec<string> & names,
	       int nBlocks = 0,
	       int myBlock = 0,
	       double line = -1.0) {
    
    svec<int> selection;
    ChunkItSelect(out, chunks, in, names, selection, nBlocks, myBlock, line);
     
  }
  
  void ChunkItSelect(vecDNAVector & out, 
		     svec<SeqChunk> & chunks, 
		     const vecDNAVector & in, 
		     const svec<string> & names,
		     const svec<int> & selection,
		     int nBlocks = 0,
		     int myBlock = 0,
		     double line = -1.0);
     


private:

  svec<int> m_lengths;
  int m_size;
  int m_overlap;


};










#endif //SEQCHUNK_H_

