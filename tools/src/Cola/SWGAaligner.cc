#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "SWGAaligner.h"

void SWGAaligner::visitNode(EditGraphNode* currNode, int  currCheckpointColIndex) {
  int i = currNode->getRow();
  int j = currNode->getCol();
  int k = currNode->getDepth(); 
  double s, score1, score2;
  switch(k) {
  case 0:  //Moving vertically
    score1 = editGraph.getNode(i-1, j, k)->getScore() + params.getGapExtP();
    score2 = editGraph.getBestScoreAtRowCol(i-1, j) + params.getGapOpenP();  
    if(score1>score2) {
      currNode->setScore(score1);
      currNode->setCPACords(editGraph.getNode(i-1, j, k),
        currCheckpointColIndex);
    } else {
      currNode->setScore(score2);
      currNode->setCPACords(editGraph.getBestNodeAtRowCol(i-1, j),
        currCheckpointColIndex);
    }
    break;
  case 1:  //Moving horizontally
    score1 = editGraph.getNode(i, j-1, k)->getScore() + params.getGapExtP();
    score2 = editGraph.getBestScoreAtRowCol(i, j-1) + params.getGapOpenP();  
    if(score1>score2) {
      currNode->setScore(score1);
      currNode->setCPACords(editGraph.getNode(i, j-1, k),
        currCheckpointColIndex);
    } else {
      currNode->setScore(score2);
      currNode->setCPACords(editGraph.getBestNodeAtRowCol(i, j-1),
        currCheckpointColIndex);
    }
    break;
  case 2: //Moving diagonally
    s = editGraph.getBestScoreAtRowCol(i-1,j-1);
    if( getTargetSeq()[j] == getQuerySeq()[i] ) {
      // The scoring is uniform for SW as opposed to NS
      if(i*j==0 && s==MINUS_INF) { s = 0; } // special case for first row/column
      currNode->setScore(s + 1);
    } else {
      currNode->setScore(s + params.getMismatchP());
    }
    currNode->setCPACords(editGraph.getBestNodeAtRowCol(i-1, j-1),
      currCheckpointColIndex);
    break;
  default: 
    currNode->setScore(MINUS_INF);
  }
} 


