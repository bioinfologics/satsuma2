#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <fstream>
#include <cmath>
#include "AlignmentCola.h"

//=====================================================================
void AlignmentCola::traceAlignment(bool tracePrefix) {
  if(pathNodes.begin()==pathNodes.end()) { return; } // Return if path is empty
  
  int currRow, currCol,prevRow, prevCol;
  
  // Start from the first node on the path and build up the alignment
  // by tracing to the end. To do this we need to compare each cell
  // with its predecessor and choose the score based on the transition
  map<int, EditGraphNode>::iterator iter = pathNodes.begin();
  currRow = iter->second.getRow()-1;
  currCol = iter->second.getCol()-1;

  // Make sure containers are reset
  info.resetCalculations(iter->second.getCol(), iter->second.getRow());
  resetContainers();

  for( ; iter != pathNodes.end(); iter++) {
    if(!tracePrefix) {
      //Reset everything if node with 0 score - This removes the 
      // prefix section of the alignment which would have a negative score.
      if(iter->second.getScore()<=0) {
        resetContainers();
        // Set target and query offsets and reset other fields
        if(++iter==pathNodes.end()) { break; }
        info.resetCalculations(iter->second.getCol(), iter->second.getRow());
        iter--;
        // Set the currRow/col to be used in next basepair alignment
        currRow = iter->second.getRow() - 1;
        currCol = iter->second.getCol() - 1;
        continue;  
      }
    }
    prevRow  = currRow;
    prevCol  = currCol; 
    currRow  = iter->second.getRow();
    currCol  = iter->second.getCol();
    
    char queryLetter  = querySeq[currRow];
    char targetLetter = targetSeq[currCol];

    if( (currRow != prevRow && currCol != prevCol) ){ // Diagonal move
      targetStr += targetLetter;
      targetIdxsInQuery[currCol] = currRow;
      queryStr += queryLetter;
      queryIdxsInTarget[currRow] = currCol;
      if( queryLetter == targetLetter) {
        matchesStr += MATCH_CHAR;
        info.baseMatched++;
      } else {
        matchesStr += MISMATCH_CHAR;  
      }
      info.tBaseAligned++;
      info.qBaseAligned++;
    } else if( currRow != prevRow) {                  // Vertical move
      targetStr += GAP_CHAR;
      queryStr += queryLetter;
      matchesStr += MISMATCH_CHAR;  
      info.qBaseAligned++;
    } else {                                          // Horizontal move
      queryStr   += GAP_CHAR;
      targetStr  += targetLetter;
      matchesStr += MISMATCH_CHAR;  
      info.tBaseAligned++;
    }
    updateSWScore(); // update the smith-waterman score
  }

  // Set score value from the last node on the path (get from iterator used above)
  info.rawScore = (--iter)->second.getScore();  
  info.alignmentLen = matchesStr.size();
} 

void AlignmentCola::updateSWScore() {
  // If the latest added bases match, add 1 to the score and otherwise subtract 1
  if(queryStr[queryStr.size()-1] == targetStr[targetStr.size()-1]) { info.smithWatermanScore++;} 
  else { info.smithWatermanScore--; }
}  

void AlignmentCola::addNodeToPath(EditGraphNode* node) {
  pathNodes[(node->getRow() + node->getCol())] = *node;
} 

void AlignmentCola::keepSubalignment(int start, int end) {
  if(start>=getLength() || (start-end)>=getLength()) { cerr<<"Error cutting down alignment"<<endl; }
  
  int counter = 0;
  map<int, EditGraphNode>::iterator iter = pathNodes.begin();
  while( iter != pathNodes.end()) {
    if(counter<start || counter>end) {  
      pathNodes.erase(iter++);
    } else { iter++; }
    counter++;  
  }
  // Retrace alignment from updated path nodes
  traceAlignment();
  
  
}

double AlignmentCola::calcPVal() const {
// TODO Currently this only considers defaults
//  and also uses conditionals which is not nice  

  //NSGA defaults
  double lambda = 0.11;
  double mu     = 5.5;

  switch(params.getType()) {
    case NS: 
      lambda = 0.13;
      mu     = 6.4;
      break;
    case SWGA: 
      lambda = 0.72;
      mu     = 9.25;
      break;
    case SW: 
      lambda = 0.06;
      mu     = 44.6;
      break;
    case NSGA:
    case UNUSED:
      break;
  }
  double x      = calcModSWScore();
  double p      = 1 - exp(-exp(-lambda*(x-mu)));
  return p;
}

void AlignmentCola::printFull(double pValLimit, ostream& sout, int screenWidth, bool withInfo) const{
  if(calcPVal()>pValLimit) { return; } // Significance is less than required
  if(withInfo) {
    sout << "**********************************************"              << endl
         << "Aligner(1:NSGA 2:NS 3:SWGA 4:SW): " << params.getType()      << endl
         << "Gap Open Penalty:                 " << params.getGapOpenP()  << endl
         << "Gap Extension Penalty:            " << params.getGapExtP()   << endl
         << "Mismatch Penalty:                 " << params.getMismatchP() << endl
         << "**********************************************"              << endl;
  }
  Alignment::printFull(pValLimit, sout, screenWidth);
}
 
