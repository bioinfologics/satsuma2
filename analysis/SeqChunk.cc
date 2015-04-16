#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "analysis/SeqChunk.h"
#include "util/mutil.h"


void SeqChunkSelect::Read(const string & file)
{
  m_target.clear();
  m_query.clear();
  MergeRead(file);
}

void SeqChunkSelect::Write(const string & file)
{
  CMWriteFileStream f;
  f.Open(file.c_str());
  int len = 0;

  int i = 1;
  f.Write(i);

  len = m_query.isize();
  f.Write(len);

  for (i=0; i<len; i++) {
    f.Write(m_query[i]);
  }

  len = m_target.isize();
  f.Write(len);

  for (i=0; i<len; i++) {
    f.Write(m_target[i]);
  }
}


void SeqChunkSelect::MergeRead(const string & file)
{
  CMReadFileStream f;
  f.Open(file.c_str());
  int len = 0;

  int i;
  f.Read(i);
  cout << "SeqChunkVersion " << i << endl;

  f.Read(len);
  m_query.resize(len, 0);

  for (i=0; i<len; i++) {
    int v;
    f.Read(v);
    m_query[i] += v;
  }

  f.Read(len);
  m_target.resize(len, 0);

  for (i=0; i<len; i++) {
    int v;
    f.Read(v);
    m_target[i] += v;
  }
}



void ChunkManager::ChunkItSelect(vecDNAVector & out, 
				 svec<SeqChunk> & chunks, 
				 const vecDNAVector & in, 
				 const svec<string> & names,
				 const svec<int> & selection,
				 int nBlocks,
				 int myBlock,
				 double line) 
{
  
  m_lengths.resize(in.size(), 0);
  
  int i, j;
  int n = 0;
  for (i=0; i<(int)in.size(); i++) {
    if (in[i].size() < 6)
      continue;
    int l = (int)in[i].size();
    n += 1 + l / (m_size - m_overlap);
  }
  
  
  chunks.resize(n);
  out.resize(n);
  
  
  cout << "select=" << selection.isize() << "\tchunks=" << n << endl;
  if (selection.isize() > 0 && selection.isize() != n) {
    cout << "*************************************************" << endl;      
    cout << "ERROR: chunk size and selection do not match up!!" << endl;      
  }
  
  int lastBlock = n;
  int firstBlock = 0;
  if (nBlocks > 0) {
    int bSize = (n + nBlocks - 1)/ nBlocks;
    firstBlock = myBlock * bSize;
    lastBlock = (myBlock + 1) * bSize;
    if (firstBlock == lastBlock) {
      cout << "Nonsensical block assignment!" << endl;
      lastBlock = n;
      firstBlock = 0;	
    }      
    cout << "Will process target chunks " << firstBlock << " through " << lastBlock - 1 << endl;
  }
  
  
  if (line >= -0.1) {
    int theBlock = (int)(line * (double)n + 0.5);
    cout << "#### Assignment=" << line << ", pick block " << theBlock << endl;
    firstBlock = theBlock;
    lastBlock = theBlock + 1;
  } 

  
  int k = 0;
  for (i=0; i<(int)in.size(); i++) {
    m_lengths[i] = (int)in[i].size();
    int l = (int)in[i].size();
    if (l < 6)
      continue;
    int nChunks = 1 + l / (m_size - m_overlap);
    for (j=0; j<nChunks; j++) {
      int start = j * (m_size - m_overlap);
      if (start < 0)
	start = 0;
      int end = (j + 1) * (m_size - m_overlap) + m_overlap;
      if (end >= l)
	end = l;
      
      //cout << "Selection i=" << i << " k=" << k << " sel=" << selection[k] << endl;
      
      if (selection.isize() > 0) {
	if (selection[k] > 0) {
	  out[k].SetToSubOf(in[i], start, end-start);
	}
      } else {
	if (k >= firstBlock && k < lastBlock) {
	  out[k].SetToSubOf(in[i], start, end-start);
	  int n = CountNs(out[k]);
	  if (n >= out[k].isize())
	    out[k].clear();
	}
      }
      
      
      
      chunks[k].Set(names[i], start, i);
      k++;
      
    }
  }
  cout << "chunks: " << n << endl;
}

