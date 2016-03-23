#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "EditGraph.h"

//=====================================================================
void EditGraphNode::setCPACords(EditGraphNode* parent, int midCol) { 
  if(col == midCol) {
    CPArow   = row;
    CPAdepth = depth;
  } else if(col > midCol) {
    CPArow   = parent->getCPARow(); 
    CPAdepth = parent->getCPADepth(); 
  }
}

//=====================================================================
void EditGraphDepth::init(int row, int col) {
  for( int depth=0; depth<getSize(); depth++) {
    EditGraphNode* currNode = getNode(depth);
    // If a contiguity cell score is MINUS_INF higher levels will be inf for sure
    if(depth>3 &&  currNode->getScore()==MINUS_INF) { break; }
    else { currNode->setScore(MINUS_INF); } 
  }
}


//=====================================================================
void EditGraph::initCol(int col, int startRow, int endRow) {
  // Only need to reset nodes that fall within the bandwidth boundaries for banded alignment
  int start = max(startRow-1, col-bandWidth-1);
  int end   = min(endRow, col+bandWidth+1);
  for(int row=start; row<=end; row++) {
    getCell(row, col)->setBestNode(0); //Reset the bestNode index
    getCell(row, col)->init(row, col);
  }
}

void EditGraph::updateBest(EditGraphNode* bNode) {
  if(!bNode) { 
    bestScoredNode.setScore(MINUS_INF);
  } else if( bNode->getScore() > bestScoredNode.getScore()) {
    bestScoredNode = *bNode;
  }
}

void EditGraph::checkPoint(int col, int maxStartRow, int minEndRow) {
  // Only need to save nodes that fall within the bandwidth boundaries for banded alignment
  EditGraphColumn* colDat = getColumn(col);
  for(int row=maxStartRow; row<=minEndRow; row++) {
    checkpointCol.setCell(row, colDat->getCell(row));
  }
}
