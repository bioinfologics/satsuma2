#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "NSGAaligner.h"

void NSGAaligner::visitNode(EditGraphNode* currNode, int  currCheckpointColIndex) {
  if (currNode->getDepth() < 2) { // Can reuse affine zero-contiguity function
    SWGAaligner::visitNode(currNode, currCheckpointColIndex);
  } else {
    visitNodeCola(currNode, currCheckpointColIndex);
  }
} 

void NSGAaligner::visitNodeCola(EditGraphNode* currNode, int  currCheckpointColIndex) {
  // Choose the best vertical, horizontal (k=0,1) and diagonal move with zero depth of contig  
  int i = currNode->getRow();
  int j = currNode->getCol();
  int k = currNode->getDepth(); 
  if(k==2) {
    double score1, score2 = MINUS_INF;
    score1 = editGraph.getBestScoreAtRowCol(i, j);
    // Moving from the diagonal neighbour
    score2 = editGraph.getBestScoreAtRowCol(i-1, j-1) + params.getMismatchP();  
    if ( score2 >= score1 ) {
      currNode->setScore(score2);
      currNode->setCPACords(editGraph.getBestNodeAtRowCol(i-1, j-1),
          currCheckpointColIndex);
    } else {
      currNode->setScore(score1);
      // Pass in -1 for currCheckpoint Col to ensure that the parent nodes CPA gets set
      currNode->setCPACords(editGraph.getBestNodeAtRowCol(i,j), -1);
    }
  } else {
    // Moving from the diagonal neighbour only if there's  a  match.
    if( getTargetSeq()[j] == getQuerySeq()[i] ) {
      double s = editGraph.getNode(i-1,j-1,k-1)->getScore();
      if(i*j==0 && k==3 && s==MINUS_INF) { s = 0; } // special case for first row/column
      if (s != MINUS_INF) { 
        int depth = k - 2;
        currNode->setScore(s - pow(((depth-1)-0)/1.0, 3) + pow((depth-0)/1.0, 3));  
//        currNode->setScore(s - pow(1.9,((depth-1)-0)/1) + pow(1.9,(depth-0)/1));  
//        currNode->setScore(s + depth);  
        currNode->setCPACords(editGraph.getNode(i-1, j-1, k-1), currCheckpointColIndex);
      }
    } else {
      currNode->setScore(MINUS_INF);// There is no match:
    } 
  }
}
